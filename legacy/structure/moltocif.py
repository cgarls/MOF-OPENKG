from rdkit import Chem
from rdkit.Chem import AllChem
import re
import os
import datetime


rate = 40
atoms = []
bonds = []

class Atom:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if 'mol' in self.__dict__:
            x, y, z, self.t = self.mol.split()[:4]
            self.x = round(float(x)/rate, 5)
            self.y = round(float(y)/rate, 5)
            self.z = round(float(z)/rate, 5)
        atoms.append(self)

    @property
    def name(self):
        return f"{self.t}{self.id}"

    def cif(self):
        return "{:7}{:6}{:9}{:9}{:9}  0.00000  Uiso   1.00".format(self.name, self.t, self.x, self.y, self.z)

    def distance(self, atom):
        return round(((self.x - atom.x) ** 2 + (self.y - atom.y) ** 2 + (self.z - atom.z) ** 2) ** 0.5 * rate, 5)


class Bond:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if 'mol' in self.__dict__:

            first = self.mol.split()[0]
            if int(first) >= len(atoms):
                i = 1
                while i < len(first):
                    # 最短前缀作为第一个原子序号
                    if int(first[:i]) <= len(atoms) and int(first[i:]) <= len(atoms):
                        a = first[:i]
                        b = first[i:]
                        break
                    i += 1
                t = self.mol.split()[1]
            else:
                a, b, t = self.mol.split()[:3]
            self.a = int(a) - 1
            self.b = int(b) - 1
            self.t = {'1': 'S',
                '2': 'D',
                '3': 'T'}[t]
        bonds.append(self)

    def cif(self):
        # print(len(atoms), self.a, self.b)
        a = atoms[self.a]
        b = atoms[self.b]
        return "{:7}{:7}{:9}  .   {}".format(a.name, b.name, a.distance(b), self.t)


def mol2cif(path, cif_path=None):
    if cif_path is None:
        cif_path = path.replace('.mol', '.cif').replace('.sdf', '.cif')
    print(f'converting {path} -> {cif_path}')

    if path.endswith('.mol'):
        tmp_path = path.replace('.mol', '.txt')
        mol = Chem.MolFromMolFile(path)
        m2 = Chem.AddHs(mol) # 加氢原子
        AllChem.EmbedMolecule(m2) # 2D->3D化
        AllChem.MMFFOptimizeMolecule(m2)   # 使用MMFF94最小化RDKit生成的构象
        #m3 = Chem.RemoveHs(m2) # 删除氢原子
        #print(Chem.MolToMolBlock(m3))
        m2.SetProp('_Name', 'cyclobutane')
        Chem.MolToMolFile(m2, tmp_path)
    else:
        tmp_path = path

    with open(tmp_path,'r',encoding='utf-8') as f:
        content = f.read()
        lines = content.split('\n')
    if tmp_path != path:
        os.remove(tmp_path)

    globals()['atoms'] = []
    globals()['bonds'] = []
    # print(len(atoms), len(bonds))
    for line in lines[2:]:
        if line.split().count('0') > 8:
            Atom(mol=line, id=len(atoms)+1)
        elif line.startswith('M  END'):
            break
        elif re.findall(r'[A-Za-z]', line) or len(line.strip())==0:
            pass
        else:
            Bond(mol=line)
    
    # print(atoms, bonds)

    s = f'''data_cif_{os.path.split(path)[-1]}
_audit_creation_date              {datetime.date.today()}
_audit_creation_method            'Materials Studio'
_symmetry_space_group_name_H-M    'P1'
_symmetry_Int_Tables_number       1
_symmetry_cell_setting            triclinic
loop_
_symmetry_equiv_pos_as_xyz
  x,y,z
_cell_length_a                    {rate}.0000
_cell_length_b                    {rate}.0000
_cell_length_c                    {rate}.0000
_cell_angle_alpha                 90.0000
_cell_angle_beta                  90.0000
_cell_angle_gamma                 90.0000
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_occupancy
''' + '\n'.join([x.cif() for x in atoms]) + '''
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_geom_bond_site_symmetry_2
_ccdc_geom_bond_type
''' + '\n'.join([x.cif() for x in bonds]) + '\n'
    with open(cif_path, 'w', encoding='utf-8') as f:
        f.write(s)


    # 在len(bonds) < len(atoms)-1 的情况下，图不连通，需要重建边，仅在原有的基础上加！！！

    return cif_path

if __name__ == "__main__":
    i = 62651
    mol2cif(f'sdf/{i}.sdf', f'{i}.cif')