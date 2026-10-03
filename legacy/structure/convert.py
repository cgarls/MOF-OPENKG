import enum
from tokenize import Double
from pymatgen.core.structure import Structure
from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
from pymatgen.io.ase import AseAtomsAdaptor
import os
from ase import io
import ase
from scipy.stats import rankdata
import numpy as np
from torch_geometric.utils import dense_to_sparse
import pandas as pd
from pymatgen.analysis import structure_matcher
import getpass
from pymatgen.core import Lattice, Structure, Molecule
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from moltocif import mol2cif
import json
import grakel


tmp_path = f'/tmp/{getpass.getuser()}/structure/'
os.makedirs(tmp_path, exist_ok=True)


processing_args = {
    "graph_max_radius": 8.0,
    "graph_max_neighbors": 12,
}


def get_matrix_trimmed(ase_crystal):
    distance_matrix = ase_crystal.get_all_distances(mic=True)

    ##Create sparse graph from distance matrix
    return threshold_sort(
        distance_matrix,
        processing_args["graph_max_radius"],
        processing_args["graph_max_neighbors"],
        adj=False,
    )

def get_matrix_voronoi(ase_crystal):
    Converter = AseAtomsAdaptor()
    pymatgen_crystal = Converter.get_structure(ase_crystal)
    # double check if cutoff distance does anything
    Voronoi = VoronoiConnectivity(
        pymatgen_crystal, cutoff=processing_args["graph_max_radius"]
    )
    connections = Voronoi.max_connectivity

    return threshold_sort(
        connections,
        9999,
        processing_args["graph_max_neighbors"],
        reverse=True,
        adj=False,
    )

def threshold_sort(matrix, threshold, neighbors, reverse=False, adj=False):
    mask = matrix > threshold
    distance_matrix_trimmed = np.ma.array(matrix, mask=mask)
    if reverse == False:
        distance_matrix_trimmed = rankdata(
            distance_matrix_trimmed, method="ordinal", axis=1
        )
    elif reverse == True:
        distance_matrix_trimmed = rankdata(
            distance_matrix_trimmed * -1, method="ordinal", axis=1
        )
    distance_matrix_trimmed = np.nan_to_num(
        np.where(mask, np.nan, distance_matrix_trimmed)
    )
    distance_matrix_trimmed[distance_matrix_trimmed > neighbors + 1] = 0

    if adj == False:
        distance_matrix_trimmed = np.where(
            distance_matrix_trimmed == 0, distance_matrix_trimmed, matrix
        )
        return distance_matrix_trimmed
    elif adj == True:
        adj_list = np.zeros((matrix.shape[0], neighbors + 1))
        adj_attr = np.zeros((matrix.shape[0], neighbors + 1))
        for i in range(0, matrix.shape[0]):
            temp = np.where(distance_matrix_trimmed[i] != 0)[0]
            adj_list[i, :] = np.pad(
                temp,
                pad_width=(0, neighbors + 1 - len(temp)),
                mode="constant",
                constant_values=0,
            )
            adj_attr[i, :] = matrix[i, adj_list[i, :].astype(int)]
        distance_matrix_trimmed = np.where(
            distance_matrix_trimmed == 0, distance_matrix_trimmed, matrix
        )
        return distance_matrix_trimmed, adj_list, adj_attr


def dense_to_sparse(adj):
    r"""Converts a dense adjacency matrix to a sparse adjacency matrix defined
    by edge indices and edge attributes.

    Args:
        adj (Tensor): The dense adjacency matrix.
     :rtype: (:class:`LongTensor`, :class:`Tensor`)
    """
    # assert adj.shpae() >= 2 and adj.shape() <= 3
    assert adj.shape[-1] == adj.shape[-2]
    index = adj.nonzero()
    edge_attr = adj[index]

    return index, edge_attr


def get_unified_structure(path, typ='trimmed'):
    # os.path.join(
    #         data_path, structure_id + "." + processing_args["data_format"]
    #     )
    ase_crystal = ase.io.read(path)
    if typ == 'trimmed':
        matrix = get_matrix_trimmed(ase_crystal)
    else:
        matrix = get_matrix_voronoi(ase_crystal)

    # print(ase_crystal.atomic_numbers)

    index, edge_attr = dense_to_sparse(matrix)
    print(index, edge_attr)

    print(index[0].shape, index[1].shape)
    print(edge_attr.shape)

    return edge_attr.shape[0]


def deduplicate(folder, new_folder):
    mofs = [] #initialize list to store Pymatgen structures
    entries = os.listdir(folder) #get all CIFs
    entries.sort() #alphabetical sort

    #for every CIF, store Pymatgen Structure in list
    for entry in entries:

        if '_pre.cif' not in entry:
            continue
        
        #read CIF
        mof_temp = Structure.from_file(os.path.join(folder,entry),primitive=False)

        #tag Pymatgen structure with its name
        mof_temp.name = entry
        mofs.append(mof_temp)

    #Initialize StructureMatcher
    sm = structure_matcher.StructureMatcher(primitive_cell=True)

    #Group structures
    groups = sm.group_structures(mofs)
    #Write out set of only unique CIFs
    if not os.path.exists(new_folder):
        os.mkdir(new_folder)
    for group in groups:
        mof_temp = group[0]
        mof_temp.to(filename=os.path.join(new_folder,mof_temp.name))

    return str(len(groups))+' unique out of '+str(len(entries))+' total'


def test1():
    p1 = '/mnt/data2/sy/Desktop/mofs/mofs-kg/structure/data/VUVSES.cif'
    p2 = '/mnt/data2/sy/Desktop/mofs/mofs-kg/structure/data/VUVSES_trim.cif'
    
    df = pd.DataFrame(columns=['original', 'trim', 'original_voronoi', 'trim_voronoi'])
    for radius in [2, 4, 8]:
        for neighbors in [8, 12]:
            processing_args['graph_max_radius'] = radius
            processing_args['graph_max_neighbors'] = neighbors
            df.loc[f'radius={radius}, neighbors={neighbors}'] = {
                'original': get_unified_structure(p1),
                'trim':  get_unified_structure(p2),
                'original_voronoi': get_unified_structure(p1, typ="voronoi"),
                'trim_voronoi': get_unified_structure(p2, typ="voronoi")
            }

    print(df)
    df.to_csv('data/test.csv') 


def test2():
    f = 'VUVSES_trim'
    path = f'/mnt/data2/sy/Desktop/mofs/mofs-kg/structure/data/{f}.cif'
    pre_path = tmp_path + path.split('/')[-1]

    if not os.path.exists(pre_path):
        mof = io.read(path)
        io.write(pre_path, mof)
    
    structure = Structure.from_file(pre_path)
    molecule = Molecule.from_sites(structure.sites)

    df = pd.DataFrame(columns=['len(bonds)'])
    for  tol in range(0, 50):
        tol /= 100
        bonds = molecule.get_covalent_bonds(tol=tol)
        df.loc[tol, 'len(bonds)'] = len(bonds)
        print(tol, len(bonds))
    
    df.plot()
    plt.savefig('data/bonds.png')


def is_H(s):
    if s[0] != 'H':
        return False
    if len(s) == 1:
        return True
    try: 
        int(s[1:])
        return True
    except:
        return False

def is_int(s):
    try: 
        int(s)
        return True
    except:
        return False

def is_outer_atom(s):
    if s[0] == '_' or s[-1] == '_':
        return False
    if 'A' <= s[0] <= 'Z':
        return is_H(s) or s[-1] in ['A', 'B', 'E'] or s.count('_') and is_int(s.split('_')[-1])
    return False


def get_relations(path, need_bonding_number=True):
    if not path.endswith('.cif') and not path.endswith('.sdf') and not path.endswith('.mol'):
        raise Exception('Can only deal with cif / sdf / mol files!')
    
    # SDF结构一定存在键信息，先转换为cif
    if path.endswith('.sdf') or path.endswith('.mol'):
        tmp = tmp_path + path.split('/')[-1].replace('.sdf', '.cif').replace('.mol', '.cif')
        path = mol2cif(path, tmp)

    with open(path, 'r') as f:
        content = f.read()
        lines = content.split('\n')

    # 已经存在键信息
    # 在len(bonds) < len(atoms)-1 的情况下，图不连通，需要重建边，仅在原有的基础上加！！！
    ret = {}
    if content.count('_geom_bond_atom_site'):
        i = 0
        headers = []
        while not lines[i].count('_geom_bond'):
            i += 1
        while len(lines[i]) and lines[i][0] == '_':
            headers.append(lines[i])
            i += 1
        atom1 = headers.index('_geom_bond_atom_site_label_1')
        atom2 = headers.index('_geom_bond_atom_site_label_2')
        count = None
        if '_ccdc_geom_bond_type' in headers:
            count = headers.index('_ccdc_geom_bond_type')
        while len(lines[i]) and lines[i][0] != '_' and lines[i] != 'loop_':
            lis = lines[i].split()
            if not is_outer_atom(lis[atom1]) and not is_outer_atom(lis[atom2]):
                # print(f'{lis[atom1]},{lis[atom2]}')
                if count:
                    ret[f'{lis[atom1]},{lis[atom2]}'] = {'S': 1,
                                            'D': 2,
                                            'T': 3}[lis[count]]
                else:
                    ret[f'{lis[atom1]},{lis[atom2]}'] = 0
            i += 1
    
    if not need_bonding_number:
        return list(ret.keys())
    
    # 不存在 _geom_bond，预处理cif，去除对称关系和成键信息
    i = 0
    cif_lines = []
    while i < len(lines):
        line = lines[i]
        if line.count('_symmetry'):
            i += 1
            continue
        if line.startswith('_chemical_name_systematic') or line.startswith('_chemical_properties_biological'):
            if lines[i+1][0] != '_':
                i += 2
                continue
        if line.startswith('loop_'):
            if lines[i + 1].startswith('_symmetry') or lines[i + 1].startswith('_geom_bond'):
                i += 1
                while len(lines[i]) and lines[i][0] == '_':
                    i += 1
                while len(lines[i]) and lines[i][0] != '_':
                    i += 1
                continue
        # print(line)
        lis = line.split()
        if len(lis) and is_outer_atom(lis[0]) or len(lis) > 1 and is_outer_atom(lis[1]):
            i += 1
            continue
        cif_lines.append(line)
        i += 1
    
    # print(cif_lines)
    # print(ret)
    pre_path = tmp_path + path.split('/')[-1]
    with open(pre_path, 'w') as f:
        f.write('\n'.join(cif_lines))

    print(path, pre_path)
    print(f'preprocessed: {path} -> {pre_path}')
    try:
        # 读取
        mof = io.read(pre_path)
        io.write(pre_path, mof)
        structure = Structure.from_file(pre_path)
    except Exception as e:
        print(f'Cannot open target structure!')
        return ret

    molecule = Molecule.from_sites(structure.sites)
    bonds = molecule.get_covalent_bonds(tol=0) 
    # print(bonds)

    coords_dic = {}
    atom_dic = defaultdict(int)

    for bond in bonds:
        lis = []
        for site_string in ['site1', 'site2']:
            site = bond.__dict__[site_string]
            coords = tuple(round(x, 2) for x in site.coords)
            if coords in coords_dic:
                s = coords_dic[coords]
            else:
                atom_dic[site.species_string] += 1
                s = site.species_string + str(atom_dic[site.species_string])
                coords_dic[coords] = s
            lis.append(s)
        # print(lis, f'{lis[0]},{lis[1]}')
        if f'{lis[0]},{lis[1]}' not in ret:
            print('additional bond!', lis, bond.count)
        if ret.get(f'{lis[0]},{lis[1]}', 0) == 0:
            ret[f'{lis[0]},{lis[1]}'] = bond.count 

    return ret

            


def get_bond_length():
    def is_alpha(s):
        if not s:
            return False
        for c in s.strip():
            if  'a' <= c <= 'z' or 'A' <= c <= 'Z':
                pass
            else:
                return False
        return True

    def digital(s):
        try:
            return int(s.strip())
        except:
            return 0

    def get_length(s):
        try:
            return digital(s.split()[0])
        except:
            return 0
    
    with open('data/bond_length.txt') as f:
        lines = f.readlines()
    
    dic = {}
    for i in range(len(lines)):
        line = lines[i].strip()
        if is_alpha(line):
            dic[line] = {k: get_length(lines[i+k]) for k in range(1, 4) if get_length(lines[i+k])}

    print(dic)
    return dic


def dump_relations(dir, out_dir, need_bonding_number=True):
    already_dic = {}
    index = 0
    if out_dir.count('.'):
        if os.path.exists(out_dir):
            with open(out_dir, 'r') as f:
                already_dic = json.load(f)
    else:
        os.makedirs(out_dir, exist_ok=True)
        for file in os.listdir(out_dir):
            index = max(index, int(file.split('.')[0]))
            with open(os.path.join(out_dir ,file), 'r') as f:
                already_dic = {
                    **already_dic,
                    **json.load(f)
                }
        index += 1
        
    print(f'already dump {len(already_dic)} structures!')

    new_dic = {}
    files = []
    if type(dir) == list:
        path = dir
    else:
        path = [dir]
    for p in path:
        files += [os.path.join(p, x) for x in os.listdir(p)]
        print(len(files))
    
    # exit(0)

    for f in files:
        name = f.split('/')[-1]
        if name in already_dic:
            print(name, end=' ')
            continue
        print('-'*10, name, '-'*10)
        t = get_relations(f, need_bonding_number=need_bonding_number)
        new_dic[name] = t
        if len(new_dic) % 5000 == 0:
            oup = out_dir if out_dir.count('.') else os.path.join(out_dir, f'{index}.json')
            with open(oup, 'w') as f:
                json.dump(new_dic, f)
            index += 1
            new_dic = {}

    if len(new_dic):
        oup = out_dir if out_dir.count('.') else os.path.join(out_dir, f'{index}.json')
        with open(oup, 'w') as f:
            json.dump(new_dic, f)


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


def statistic_relations(path):
    if type(path) == str:
        path = [path]

    fig, ax = plt.subplots(3, len(path), figsize=(15, 10))
    fig.subplots_adjust(hspace=0.4, wspace=0.2)
    
    for i, path_i in enumerate(path):
        dic = load_dir(path_i)

        print('total structure:', len(dic))
        count = {
            'edge': defaultdict(int),
            'bonding': defaultdict(int),
            'node': defaultdict(int)
        }

        for edges in dic.values():
            count['edge'][len(edges)] += 1
            count['bonding'][sum(edges.values())] += 1
            nodes = set()
            for e in edges:
                n1, n2 = e.split(',')
                nodes.add(n1)
                nodes.add(n2)
            count['node'][len(nodes)] += 1

        for ix, t in enumerate(count.keys()):
            ax[ix][i].title.set_text(f'barchart {t}: {len(count[t])} type')
            ax[ix][i].bar(list(count[t].keys()), list(count[t].values()))

    plt.savefig('data/statistic_relations.png')


def to_kernel(path):
    with open(path, 'r') as f:
        dic = json.load(f)

    grakel.WeisfeilerLehman()


def to_graph():
    path = ['/mnt/data1/csd/cif/all', '/mnt/data2/sy/Desktop/mofs/mofs-kg/structure/sdf']
    dump_relations(path[0], 'data/csd_bond.json', False)
    dump_relations(path[1], 'data/precursor_bond.json', False)


if __name__ == "__main__":
    

    # unified = '/mnt/data1/csd/cif/unified'
    # os.makedirs(unified, exist_ok=True)

    # compare()

    # deduplicate('data', 'out')
    to_graph()
    # statistic_relations(['data/csd_bond', 'data/precursor_bond'])
    
