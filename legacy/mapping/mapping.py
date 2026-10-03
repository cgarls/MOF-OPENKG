import json
import re
from utils import *
import csv
import pandas as pd
from collections import defaultdict
from tqdm import tqdm

# mof_file = '../data/KAIST.json'
mof_file = './label/label.json'
csv_input_path = "/mnt/data1/csd/data.csv"
threads = 50

def csd_mapping(df = pd.read_csv(open(csv_input_path, 'r', encoding="utf-8-sig"))):
    citation_pattern = re.compile(
        r"Citation\(authors='(.*)', journal='Journal\((.*)\)', volume='(.*)', year=(\d*), first_page='(.*)', doi='(.*)'\)")

    df.reset_index(inplace=True)
    for i in tqdm(range(len(df))):
        row = df.loc[i]
        mof_node = MofNode(row)
        
        if type(row['color']) == str:
            color = NameNode(row['color'], 'Color')
            SingleEdge(mof_node.id, color.id, 'hasColor').update()
        
        if type(row['habit']) == str:
            habit = NameNode(row['habit'], 'Habit')
            SingleEdge(mof_node.id, habit.id, 'hasHabit').update()

        # if type(row['solvent']) == str:
        #     # print(row['solvent'])
        #     # CSD数据库中出现或者：chloroform/hexane，需要拆开分析
        #     for solvent in row['solvent'].split('/'):
        #         for cid in get_precursor_id(solvent):
        #             IdNode(cid, solvent, 'Solvent')
        #             PrecursorEdge(source=mof_node.id, target=cid, name=solvent, type='hasSolvent', origin='CSD').update()

        # 首次出版doi，有doi字段按doi，doi字段为空在publication中找
        first_publication_doi = ''
        if row['doi'] != '':
            first_publication_doi = row['doi']
        else:
            first_publication = re.search(citation_pattern, row['publication'])
            if first_publication is not None:
                first_publication_doi = first_publication.group(6)

        for publication in re.findall(citation_pattern, row['publications']):
            authors_group_list = publication[0].split(", ")
            journal = publication[1]
            volume = publication[2]
            year = publication[3]
            first_page = publication[4]
            doi = publication[5]

            paper_id = 'P0' + get_hex(doi)[-W:]
            SingleEdge(mof_node.id, paper_id, 'hasPaper', is_first_publish=doi != '' and doi == first_publication_doi).update()
            
            if paper_id not in nodes:
                paper_node = NameNode(doi, 'Paper')
                for author in authors_group_list:
                    author_node = NameNode(author, 'Author')
                    SingleEdge(paper_node.id, author_node.id, 'hasAuthor').update()

                journal_node = NameNode(journal, 'Journal')
                SingleEdge(paper_node.id, journal_node.id, 'hasJournal', \
                        volume=volume, year=year, first_page=first_page).update()


def set_dic(dic, **kwargs):
    try: 
        dic.set(**kwargs)
    except:
        dic.update(kwargs)


def precursor_mapping(mofs=json.load(open(mof_file))):
    mapping = {
        'M_precursor': 'Metal',
        'S_precursor': 'Solvent',
        'O_precursor': 'Linker'
    }

    for mof in tqdm(mofs):
        # print(mof)
        origin = mof.get('source', 'LABEL')     # 默认是标注
        # origin = mof.get('source', 'KASIT')     # 默认是标注
        mof_id = 'M0' + mof.get('identifier')
        set_dic(nodes[mof_id], Yield=mof.get('Yield', ''),
            melting_point=mof.get('melting_point', ''),
            decomposition_temperature=mof.get('decomposition_temperature', ''),
            boiling_point=mof.get('boiling_point', ''),
            time=mof.get('time', ''),
            temperature=mof.get('temperature', ''))
        
        if mof.get('method', ''):
            method_id = 'M1' + get_hex(mof.get('method', ''))[-W:]
            SingleEdge(mof_id, method_id, 'hasMethod').update()
            IdNode(method_id, mof.get('method', ''), 'Method')

        # todo cid not exist
        for k, v in mapping.items():
            for p in mof.get(k, []):
                for cid in  get_precursor_id(p.get('name')):
                    IdNode(cid, p.get('name'), 'Precursor', formula=p.get('formula'), smiles=p.get('smiles'), **{v: True})
                    PrecursorEdge(source=mof_id, target=cid, name=p.get('name'), 
                        composition=p.get('composition'), type='has' + v, 
                        origin=origin, confidence=p.get('confidence')).update()

        for p in mof.get('operation', []):
            name = 'undefined' if not p.get('name') else p.get('name')
            op_node = NameNode(name=name, type='Operation')
            SingleEdge(source=mof_id, target=op_node.id, type='hasOperation', 
                Temperature=p.get('Temperature', ''), Time=p.get('Time', ''), 
                Pressure=p.get('Pressure', '')).update()


def structure_mapping():
    all_bond = {
        'M0': load_dir('../structure/data/csd_bond'),
        'P1': load_dir('../structure/data/precursor_bond')
    }

    for el, number in element.items():
        IdNode('K1' + str(number).rjust(2, '0'), el, 'Kernel', level=1)

    for key, v in all_bond.items():
        for mof, bonds in v.items():
            # print(mof)
            neighbor = defaultdict(list)

            # 添加 Bond 与 hasBond 关系，Bond与 1-kernel关系
            for bond in bonds:
                try:
                    a1, a2 = re.findall(r'[A-Z][a-z]*\d', bond)
                except:
                    print('failed to convert:', mof, bond, a1, a2)
                    continue
                a, b = bond.split(',')
                neighbor[a].append(b)
                neighbor[b].append(a)
                b1, b2 = element[a1[:-1]], element[a2[:-1]]
                c1, c2 = min(b1, b2), max(b1, b2)
                d1, d2 = str(c1).rjust(2, '0'), str(c2).rjust(2, '0')
                bond_id = 'B0' + d1 + d2
                MultipleEdge(key + mof.split('.')[0], bond_id, 'hasBond').update()
                if bond_id not in edges:
                    IdNode(bond_id, relement[c1] + ',' + relement[c2], 'Bond')
                    SingleEdge(bond_id, 'K1' + d1, 'hasKernel').update()
                    SingleEdge(bond_id, 'K1' + d2, 'hasKernel').update()

            # 对每个结点初始编码
            # print(neighbor)
            encoding = defaultdict(dict)
            encoding[1] = {k: 'K1' + get_element(re.findall(r'[A-Z][a-z]*', k)[0])  for k in neighbor}
            for k in range(2, 4):
                # print(f'generating {k}-kernel...')
                for s, ts in neighbor.items():
                    new_encode = get_hex(encoding[k-1][s]) + \
                        str(sum([get_hash(encoding[k-1][t]) for t in ts]))
                    
                    node_id = 'K' + str(k) + get_hex(new_encode)[-W:]
                    encoding[k][s] = node_id

                    MultipleEdge(key + mof.split('.')[0], node_id, 'hasKernel').update()
                    
                    if node_id not in nodes:
                        IdNode(node_id, '', 'Kernel', level=k)
                        MultipleEdge(node_id, encoding[k-1][s], 'subKernel', center=True).update()
                        for t in ts:
                            MultipleEdge(node_id, encoding[k-1][t], 'subKernel').update()

def wrap_csd_mapping():
    # 多进程    
    import math
    from multiprocessing import Process

    df = pd.read_csv(open(csv_input_path, 'r', encoding="utf-8-sig"))
    pool = []  # 进程池
    step = math.ceil(len(df)/threads)
    for ix, start in enumerate(range(0, len(df), step)):
        p = Process(target=csd_mapping, args=(df.loc[start: start + step],))
        pool.append(p)
    for p in pool:
        p.start()
    for p in pool:
        p.join()


def wrap_precursor_mapping():
    # with open('data/csd_mapping.json', 'r') as f:
    #     res = json.load(f)
    # globals()['nodes'] = {x['id']: x for x in res['nodes']}
    # globals()['edges'] = {','.join([x['source'], x['target'], x['type']]) : x for x in res['edges']}

    # 多进程    
    import math
    from multiprocessing import Process

    mofs = json.load(open(mof_file))
    pool = []  # 进程池
    step = math.ceil(len(mofs)/threads)
    for ix, start in enumerate(range(0, len(mofs), step)):
        p = Process(target=precursor_mapping, args=(mofs[start: start + step],))
        pool.append(p)
    for p in pool:
        p.start()
    for p in pool:
        p.join()


    # precursor_mapping()



def mapping(typ='synthesis'):

    if typ == 'csd':
        oup = 'data/csd_mapping.json'
        csd_mapping()

    elif typ == 'synthesis':
        oup = 'data/synthesis_mapping_label.json'
        # oup = 'data/synthesis_mapping.json'
        print('===================csd================')
        # wrap_csd_mapping()
        csd_mapping()

        print('nodes', len(nodes))
        print('edges', len(edges))

        print('===================precursor================')
        # wrap_precursor_mapping()
        precursor_mapping()

    elif typ == 'structure':
        oup = 'data/structure_mapping.json'
        structure_mapping()

    elif typ == 'all':
        oup = 'data/mapping.json'
        precursor_mapping()
        csd_mapping()
        structure_mapping()

    res = {
        "nodes": [get_dic(node) for node in list(nodes.values())], 
        "edges": [get_dic(edge) for edge in list(edges.values())]
    }

    with open(oup, "w") as f:
        json.dump(res, f, ensure_ascii=False, cls=NumpyEncoder, indent=2)

        
def get_mofs():
    csd_mof = os.listdir('/mnt/data1/csd/cif/all')
    csd_mof = set([x.split('.')[0] for x in csd_mof])
    print('csd mof', len(csd_mof))


    with open('../data/KAIST.json') as mf:
        mofs = json.load(mf)

    kaist_mof = set()
    for mof in mofs:
        name = mof.get('identifier')
        if name in csd_mof:
            kaist_mof.add(name)
    
    with open('data/KAIST_MOF.txt', 'w') as f:
        f.write('\n'.join(kaist_mof))
    print('kaist mof:', len(kaist_mof))

    qmof_mof = set()  
    df = pd.read_csv('/mnt/data1/qmof/qmof_database/qmof.csv')
    for name in df['name']:
        name = name.replace('_FSR', '')
        if name.startswith('core_'):
            name = name.split('_')[1]
        if name in csd_mof:
            qmof_mof.add(name)
    print('qmof mof:', len(qmof_mof))
    
    core_mof = set()
    df = pd.read_excel('/mnt/data1/core/je9b00835_si_003.xlsx')
    for name in df['filename']:
        name = name.replace('_clean', '')
        if name in csd_mof:
            core_mof.add(name)
    for name in df['Matched_CSD_of_CoRE']:
        if name in csd_mof:
            core_mof.add(name)
    print('core mof:', len(core_mof))

    from matplotlib import pyplot as plt
    from matplotlib_venn import venn3, venn2

    lp_mof = kaist_mof.union(qmof_mof).union(core_mof)
    pf_mof = csd_mof.difference(qmof_mof)
    # venn3(subsets=[kaist_mof, qmof_mof, core_mof], set_labels=('KAIST', 'QMOF', 'CoRE'), set_colors=('r', 'g', 'b'))
    # plt.savefig('data/venn.png')
    # venn3(subsets=[lp_mof, pf_mof, csd_mof], set_labels=('Synthesis', 'Performance', 'CSD'), set_colors=('r', 'g', 'b'))
    venn2(subsets=[lp_mof, pf_mof], set_labels=('Synthesis', 'Performance'), set_colors=('r', 'b'))
    plt.savefig('data/venn2.png')

    all_mof = kaist_mof.union(qmof_mof).union(core_mof)
    with open('data/ALL_MOF.txt', 'w') as f:
        f.write('\n'.join(all_mof))
    print('all mof:', len(all_mof))


    save_ids()



if __name__ == "__main__":
    mapping()
    # get_mofs()
