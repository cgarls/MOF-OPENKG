from utils import *
from collections import defaultdict


def plot_sep(type_dic, prefix=''):
    type_dic['nodes'] = sort_values(type_dic['nodes'])
    type_dic['edges'] = sort_values(type_dic['edges'])
    print(type_dic)
    plot_dic({'nodes': type_dic['nodes']}, f'data/draws/{prefix}nodes.png')
    plot_dic({'edges': type_dic['edges']}, f'data/draws/{prefix}edges.png')


def statistic(graph_index=0):
    dic = get_all_dic(graph_index)
    type_dic = {
        'nodes': defaultdict(int),
        'edges': defaultdict(int)
    }

    for t in ['nodes', 'edges']:
        for n in dic[t]:
            type_dic[t][n['type']] += 1
    plot_sep(type_dic)


def detail_statistic(graph_index=0):
    dic = get_all_dic(graph_index)
    type_dic = {
        'nodes': defaultdict(int),
        'edges': defaultdict(int)
    }

    node_dic = {}
    for n in dic['nodes']:
        node_dic[n['id']] = n

    for t in ['nodes', 'edges']:
        for n in dic[t]:
            if n['type'] == 'Precursor':
                for tt in ['Solvent', 'Linker', 'Metal']:
                    if n.get(tt):
                        type_dic[t][tt + '-' + n['id'][1]] += 1
            elif n['type'] == 'Kernel':
                type_dic[t]['Kernel-' + n['id'][1]] += 1
            elif n['type'] == 'hasKernel':
                type_dic[t]['hasKernel-' + n['target'][1]] += 1
            elif n['type'] == 'subKernel':
                type_dic[t]['subKernel-cen' if n.get('center') else 'subKernel'] += 1
            elif n['type'] in ['hasSolvent', 'hasLinker', 'hasMetal']:
                target_node = node_dic[n['target']]
                type_dic[t][n['type'] + '-' + target_node['id'][1]] += 1
            else:
                type_dic[t][n['type']] += 1
    plot_sep(type_dic, 'detail_')


def kernel_type_statistic():
    with open('data/structure_mapping.json', 'r') as f:
        dic = json.load(f)

    type_dic = {
        'Kernel': defaultdict(int),
        'center': defaultdict(int),
        'subtree': defaultdict(int),
        'hasKernel': defaultdict(int)
    }

    for n in dic['nodes']:
        if n['type'] == 'Kernel':
            type_dic['Kernel'][n['level']] += 1
    
    for e in dic['edges']:
        if e['type'] == 'center':
            type_dic['center'][e['target'][:2]] += 1
        if e['type'] == 'subtree':
            type_dic['subtree'][e['target'][:2]] += 1
        if e['type'] == 'hasKernel':
            type_dic['hasKernel'][e['target'][:2]] += 1
    
    print(type_dic)
    plot_dic(type_dic, f'data/draws/kernel_type.png', size=(6, 4), label=False)
    

def precursor_type_statistic():
    with open('data/precursor_mapping.json', 'r')  as f:
        dic = json.load(f)

    type_dic = {
        'Linker': defaultdict(int),
        'Solvent': defaultdict(int),
        'Metal': defaultdict(int)
    }
    for t in type_dic:
        type_dic['has' + t] = defaultdict(int)

    for n in dic['nodes']:
        if n['type'] in type_dic.keys():
            for k in type_dic:
                if n.get(k):
                    type_dic[k][n['id'][:2]] += 1

    for n in dic['edges']:
        if n['type'] in type_dic.keys():
            type_dic[n['type']][n['id'][:2]] += 1
    
    print(type_dic)
    plot_dic(type_dic, f'data/draws/precursor_type.png', size=(6, 4), label=False)


def distribution(detail=True):
    dic = get_all_dic()
    type_dic = {}

    def add_item(type, id):
        if type not in type_dic:
            type_dic[type] = defaultdict(int)
        type_dic[type][id] += 1
    
    def plot_all():
        dic = {}
        for k, v in type_dic.items():
            dic[k] = sort_values(v)
        prefix = 'detail_' if detail else ''
        plot_dic(dic, f'data/draws/{prefix}distribution.png', rotation=90, hspace=1.2, label=False, limit=100, size=(15, 5))

    node_dic = {}
    for n in dic['nodes']:
        node_dic[n['id']] = n

    if not detail:
        for n in dic['edges']:
            add_item(n['type'], n['target'])
        plot_all()
        return

    for n in dic['edges']:
        if n['type'] == 'hasKernel':
            add_item('hasKernel-' + n['target'][1], n['target'])
        elif n['type'] == 'subKernel':
            type = 'subKernel-cen' if n.get('center') else 'subKernel'
            add_item(type, n['target'])
        elif n['type'] in ['hasSolvent', 'hasLinker', 'hasMetal']:
            target_node = node_dic[n['target']]
            add_item(n['type'] + '-' + target_node['id'][1], n['target'])
    plot_all()


def multiple_edge_statistic():
    dic = get_all_dic()
    
    edge_dic = defaultdict(list)
    for e in dic['edges']:
        eid = e['source'] + e['target']
        edge_dic[eid].append(e)
    with open('data/multiple_edge.json', 'w') as f:
        json.dump({k: v for k, v in edge_dic.items() if len(v) > 1}, f, ensure_ascii=False, indent=2)

    type_dic = defaultdict(int)

    for k, v in edge_dic.items():
        if len(v) > 1:
            t = v[0]['type']
            for item in v:
                if item['type'] != t:
                    print(v)
                    type_dic['multipleType'] += 1
            type_dic[t] += 1

    plot_dic({'multiple_edge': type_dic}, 'data/draws/multiple_edge.png')


if __name__ == "__main__":
    graph_index = 3
    statistic(graph_index)
    detail_statistic(graph_index)
    
    # kernel_type_statistic()
    # precursor_type_statistic()
    
    # distribution(False)
    # distribution(True)

    # multiple_edge_statistic()