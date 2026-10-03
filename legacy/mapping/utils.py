import hashlib
import re
import numpy
import json
import pandas as pd
import os
import requests

nodes = {}
edges = {}
W = 10
p_file = '../structure/name_id.json'
with open(p_file, 'r') as f:
    p_ids = json.load(f)
    if 404 not in p_ids:
        p_ids[404] = []
    p_ids[404] = set(p_ids[404])

class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """

    def default(self, obj):
        if isinstance(obj, (numpy.int_, numpy.intc, numpy.intp, numpy.int8,
                            numpy.int16, numpy.int32, numpy.int64, numpy.uint8,
                            numpy.uint16, numpy.uint32, numpy.uint64)):
            return int(obj)
        elif isinstance(obj, (numpy.float_, numpy.float16, numpy.float32,
                              numpy.float64)):
            return float(obj)
        elif isinstance(obj, (numpy.ndarray,)):
            return obj.tolist()
        elif isinstance(obj, (numpy.bool_,)):
            return bool(obj)
        return json.JSONEncoder.default(self, obj)


def get_hex(s):
    return hashlib.md5(str(s).encode(encoding='UTF-8')).hexdigest()


def get_hash(s):
    return int(get_hex(s), 16)


def rm_null(dic):
    return {k: v for k,v in dic.items() if v and not pd.isnull(v)}


def get_dic(obj):
    return rm_null(obj.__dict__)


def generate_str(s):
    s = s.strip()
    for p in ['wet ', ' method', 'mild ', 'aqueous ']:
        s = s.replace(p, '')

    yield s
    yield s.replace(' ', '')
    yield s.upper()             # 完全大写
    yield s[0].upper() + s[1:]  # 首字母大写
    yield ' '.join([x[0].upper() + x[1:] for x in s.split()])   # 每个单词首字母大写

    for pattern in [r'\u22c5\d+?H2O', r'\$\d+?H2O', r'\W', r'\(.+\)', '\d:\d']:
        for p in re.findall(pattern, s):
            yield s.replace(p, '').strip()

    for p in [' at ']:
        if s.count(p):
            yield s.split(p)[0].strip()

    for p in [':', '-']:
        if s.count(p):
            yield s.split(p)[-1].strip()
            
    



def get_precursor_id(name):
    """返回根据Name检索到的列表
    """
    if name in p_ids[404]:
        return ['P2' + get_hex(name)[-W:]]

    # P1，在字典中
    for s in generate_str(name):
        if s in p_ids:
            # print(s, ids[s][0])
            return ['P1' + str(x) for x in p_ids[s]]

    # P1，可检索
    for s in generate_str(name):
        try:
            r = requests.get(f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/concepts/name/JSON?name={s}')
            if r:       # 能够查到结果
                try:
                    res = r.json()['ConceptsAndCIDs']['CID']
                    p_ids[s] = res
                    p_ids[name] = res
                    return ['P1' + str(x) for x in res]
                except:
                    pass
        except:
            pass

    # if cid not exist, use hash(precursor name)
    print('name not found:', name)
    p_ids[404].add(name)
    return ['P2' + get_hex(name)[-W:]]


def save_ids():
    p_ids[404] = list(p_ids[404])
    with open(p_file, 'w') as f:
        json.dump(p_ids, f, ensure_ascii=False)    


class NameNode:
    def __init__(self, name, type):
        self.id = type[0] + '0' + get_hex(name)[-W:]
        if self.id in nodes:
            self.__dict__.update(nodes[self.id].__dict__)

        self.type = type
        if type == 'Paper':
            self.doi = name
        else:
            self.name = name
        
        nodes[self.id] = self

    def __str__(self):
        return f"{self.type} -- id: {self.id}"

    def __hash__(self):
        return get_hash(self.id)

    def __eq__(self, other):
        return self.id == other.id


class IdNode:
    def __init__(self, id, name, type, **kwargs):
        self.id = id
        if self.id in nodes:
            self.__dict__.update(nodes[self.id].__dict__)
        
        self.__dict__.update(rm_null(kwargs))
        self.type = type
        self.name = name
        nodes[self.id] = self

    def __hash__(self):
        return get_hash(self.id)

    def __eq__(self, other):
        return self.id == other.id


class MofNode:
    def __init__(self, row):
        for field in ['calculated_density', 'ccdc_number', 'chemical_name', \
                'formula', 'has_disorder','pressure', 'temperature']:
                self.__dict__[field] = row[field]

        self.type = "MOF"
        self.id = 'M0' + row['identifier']
        assert self.id not in nodes
        nodes[self.id] = self

    def __str__(self):
        return "MofNode -- id: %s, ccdc_number: %s, chemical_name: %s, formula: %s, has_disorder: %s, " \
               "pressure: %s, temperature: %s" % (
                   self.id, self.ccdc_number, self.chemical_name, self.formula, self.has_disorder,
                   self.pressure, self.temperature)

    def set(self, **kwargs):
        for k, v in kwargs.items():
            if v:
                self.__dict__[k] = v

    def __hash__(self):
        return get_hash(self.id)

    def __eq__(self, other):
        return self.id == other.id


class BasicEdge:
    '''
    source和target直接只能存在一种该type的边
    '''
    def __init__(self, source, target, type):
        self.source = source
        self.target = target
        self.type = type

    def update(self):
        edges[self.__id__()] = self

    def __id__(self):
        return ','.join([self.source, self.target, self.type])

    def __str__(self):
        return "%s -- source: %s, target: %s" % (
            self.type, self.source, self.target)

    def __hash__(self):
        return get_hash(self.__id__())

    def __eq__(self, other):
        return self.source == other.source and self.target == other.target and self.type == other.type


class SingleEdge(BasicEdge):
    '''
    source和target直接只能存在一种该type的边，融合多个边的属性，相同id只保存一个
    '''
    def __init__(self, source, target, type, **kwargs):
        super().__init__(source, target, type)
        if self.__id__() in edges:
            self.__dict__.update(edges[self.__id__()].__dict__)

        self.__dict__.update(rm_null(kwargs))


class MultipleEdge(SingleEdge):
    '''
    source和target直接只能存在一种该type的边，多余的会增加其count
    **共用全局edges变量**
    '''
    def __init__(self, source, target, type, **kwargs):
        super().__init__(source, target, type, **kwargs)
        if self.__id__() in edges:
            if 'count' in self.__dict__:
                self.count += 1
            else:
                self.count = 2


class PrecursorEdge(BasicEdge):
    '''
    为避免类型3的多重边，对hasSolvent/hasMetel/hasLinker进行去重
    多重属性以\t分隔
    要求(source, target, type, origin)元组唯一
    '''
    def __init__(self, source, target, type, origin, **kwargs):
        super().__init__(source, target, type)
        self.origin = origin
        if self.__id__() in edges:
            self.__dict__.update(edges[self.__id__()].__dict__)
    
        for k in kwargs:
            self.__setattr(k, kwargs[k])

    def __setattr(self, k, v):
        if not v:
            return
        value = set(self.__dict__.get(k, '').split('\t'))
        for v1 in v.split('\t'):
            value.add(v1)
        value.discard('')
        self.__dict__[k] = '\t'.join(sorted(value))

    def __id__(self):
        return ','.join([self.source, self.target, self.type, self.origin])

    def __eq__(self, other):
        return self.source == other.source and self.target == other.target \
            and self.type == other.type and self.origin == other.origin


def get_map():
    with open('data/map.txt', 'r') as f:
        lines = f.read()

    dic = {}
    for line in lines.split('\n'):
        a,b  = line.split()
        dic[b] = int(a) 

    print(dic)


def load_dir(path):
    dic = {}
    if not os.path.exists(path):
        return dic
    for file in os.listdir(path):
        with open(os.path.join(path ,file), 'r') as f:
            dic = {
                **dic,
                **json.load(f)
            }
    return dic


element = {'D': 0, 'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99}
relement = {v: k for k, v in element.items()}

def get_element(t):
    return str(element[t]).rjust(2, '0')


def plot_dic(dic, path='data/statistic.png', size=(15, 6), label=True, rotation=30, limit=10, hspace=0.4):
    '''
    limit 代表 此值以下不绘制
    '''
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(len(dic), figsize=(size[0], size[1] * len(dic)))
    plt.xticks(rotation=270)
    fig.subplots_adjust(hspace=hspace, wspace=0.2)
    
    for ix, t in enumerate(list(dic.keys())):
        d = {k: v for k, v in dic[t].items() if v >= limit}

        axis = ax if len(dic) == 1 else ax[ix]
        axis.title.set_text(f'{t}: {len(d)}/{len(dic[t])} (limit:{limit})')
        rects = axis.bar(list(d.keys()), list(d.values()))
        
        if label:
            for rect in rects:  #rects 是柱子的集合
                height = rect.get_height()
                plt.text(rect.get_x() + rect.get_width() / 2, height, str(height), size=10, ha='center', va='bottom')

        for tick in axis.get_xticklabels():
            tick.set_rotation(rotation)
    
    plt.savefig(path)


def sort_values(dic):
    return {x[0]: x[1]  for x in sorted(dic.items(), key=lambda x: x[1], reverse=True)}


def get_all_dic(graph_index=0):
    '''
    获取节点集、边集，index代表哪个图
        
    '''
    graph_index_mapping={
        0: ['synthesis_mapping', 'structure_mapping'],
        1: ['synthesis_mapping'],
        2: ['structure_mapping'],
        3: ['subgraph']
    }

    dic = {'nodes':[], 'edges': []}
    for target in graph_index_mapping[graph_index]:
        with open(f'data/{target}.json', 'r') as f:
            tmp = json.load(f)
            dic['nodes'] += tmp['nodes']
            dic['edges'] += tmp['edges']

    print('节点总数：', len(dic['nodes']))
    print('边总数：', len(dic['edges']))
    return dic


if __name__ == "__main__":
    pass