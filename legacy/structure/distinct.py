from cmath import nan
import csv
import pandas as pd
from curses.ascii import isdigit

class AfterOperation:
    def __init__(self, prp_file, save_file):
        property_list = [
            'MolecularFormula', 'MolecularWeight', 'CanonicalSMILES', 'IsomericSMILES',
            'InChI', 'InChIKey', 'IUPACName', 'Title',
            'XLogP', 'ExactMass', 'MonoisotopicMass', 'TPSA',
            'Complexity', 'Charge', 'HBondDonorCount', 'HBondAcceptorCount',
            'RotatableBondCount', 'HeavyAtomCount', 'IsotopeAtomCount', 'AtomStereoCount',
            'DefinedAtomStereoCount', 'UndefinedAtomStereoCount', 'BondStereoCount', 'DefinedBondStereoCount',
            'UndefinedBondStereoCount', 'CovalentUnitCount', 'Volume3D', 'XStericQuadrupole3D',
            'YStericQuadrupole3D', 'ZStericQuadrupole3D', 'FeatureCount3D', 'FeatureAcceptorCount3D',
            'FeatureDonorCount3D', 'FeatureAnionCount3D', 'FeatureCationCount3D', 'FeatureRingCount3D',
            'FeatureHydrophobeCount3D', 'ConformerModelRMSD3D', 'EffectiveRotorCount3D', 'ConformerCount3D',
            'Fingerprint2D'
        ]
        self.header = ['CID', 'Synonym', 'SDF'] + property_list
        self.prp_file = prp_file 
        self.save_file = save_file
        # self.cid_dt = set()
        # self.prp_data = []

    def __main__(self):
        data = pd.read_csv(self.prp_file, encoding='utf-8', header=None)
        data.columns = self.header
        df = data.drop_duplicates(subset=['CID'])
        df['Synonym'] = df['Synonym'].fillna(str(['No Synonyms']))
        
        # df = df[df.CID.apply(lambda x: x.isnumeric())]

        df.to_csv(self.save_file, encoding='utf-8', index=False)
        # with open(self.prp_file, 'r', encoding='utf-8') as f:
        #     reader = csv.DictReader(f,name=self.header)
        #     for row in reader:
        #         cid = str(row['CID'])
        #         syn = row['Synonym']
        #         if not syn:
        #             row['Synonym'] = ['No Synonyms']
        #         if (isdigit(cid[0]) and not self.cid_dt.__contains__(cid)):
        #             self.prp_data.append(row)
        #             self.cid_dt.add(cid)
            
        # with open(self.save_file, 'w+', newline='', encoding='utf-8') as f:
        #     writer = csv.DictWriter(f, self.header)
        #     writer.writeheader()
        #     writer.writerows(self.prp_data)

def call(dir):
    prp_file = dir + sep + ofn
    sve_file = dir + sep + sfn
    AfterOperation(prp_file, sve_file).__main__()


if __name__ == '__main__':
    # sep = '\\'
    sep = '/'
    ofn = 'properties.csv'
    sfn = 'dt_properties.csv'
    
    call('S_precursor')
    call('M_precursor')
    call('O_precursor')
