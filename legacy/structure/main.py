import csv
import json
import requests
import time
import random
from urllib.parse import quote_plus


class MofsCrawl:
    def __init__(self, in_path, key, out_path, sdf_dir, wrong_path, log_path):
        """Initialization function.

        Args:
            in_path (str): Input file path
            key: Key
            out_path (str): Output file path
            wrong_path: Wrong file path
        """
        with open(in_path, 'r', encoding='utf-8') as f:
            self.name_dt = json.load(f).get(key)

        self.property_list = [
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
        self.property_str = ','.join(self.property_list)
        self.down_time = 1

        self.header_list = ['CID'] + ['Synonym'] + ['SDF'] + self.property_list
        # write the head of out file
        self.out_path = out_path.split('.')[0] + '.csv'
        # with open(self.out_path, 'a+', newline='', encoding='utf-8') as f:
        #    writer = csv.DictWriter(f, self.header_list)
        #    writer.writeheader()

        # set sdf dir
        self.sdf_dir = sdf_dir

        self.wrong_header = ['Name', 'Times']
        # write the head of wrong file
        self.wrong_path = wrong_path.split('.')[0] + '.csv'
        # with open(self.wrong_path, 'a+', newline='', encoding='utf-8') as f:
        #    writer = csv.DictWriter(f, self.wrong_header)
        #    writer.writeheader()
        self.wrong_dt = []

        # set query times
        self.query_time = 1

        # set log
        self.log_path = log_path.split('.')[0] + '.txt'

        # for delete the repeated
        self.cid_set = set()

    # must use cid, if you use cas number, the ',' will be regard as part of the name
    # so that you can't get many mofs per time.
    def download(self, cid_list):
        self.write_log(f'[{time.ctime()}] : -------start downloading at times '
                       f'{self.down_time}, totally {len(cid_list)} cid numbers-------')
        self.down_time += 1

        api = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'
        cid_str = ','.join(cid_list[:])
        self.write_log(f'[{time.ctime()}] : cid_str is {cid_str}')
        
        # get properties
        url = f'{api}{cid_str}/property/{self.property_str}/json'
        try:
            res = requests.get(url).json()
            prop = res['PropertyTable']['Properties']
        except KeyError:
            try:
                res = requests.get(url).json()
                prop = res['PropertyTable']['Properties']
            except KeyError:
                res = requests.get(url).json()
                prop = res['PropertyTable']['Properties']

        # get synonyms
        url = f'{api}{cid_str}/synonyms/json'
        res = requests.get(url).json()
        syn = res['InformationList']['Information']
        # for i in range(len(syn)):
        #     if 'Synonym' not in syn[i]:
        #         syn[i]['Synonym'] = []

        self.write_log(f'[{time.ctime()}] : ---properties and synonyms downloaded---')
        # cid_list = [{'CID': sy['CID']} for sy in syn]
        # download sdf
        sdf = []
        for i in range(len(cid_list)):
            url = f'{api}/{cid_list[i]}/SDF'
            sdf_file = requests.get(url, allow_redirects=True)
            save_name = f'{cid_list[i]}.sdf'
            save_dir = f'{self.sdf_dir}{save_name}'
            open(save_dir, 'wb').write(sdf_file.content)
            sdf.append({'CID': cid_list[i], 'SDF': save_name})
            if i % 17 == 0:
                self.write_log(f'[{time.ctime()}] : ---sdf file of {i + 1} mofs downloaded---')
                time.sleep(random.uniform(1.1, 2.2))

        self.write_log(f'[{time.ctime()}] : -------all {len(cid_list)} cid numbers have downloaded-------')
        # save properties
        dt = [dict(d1, **d2) for d1, d2 in zip(
            [dict(d1, **d2) for d1, d2 in zip(syn, sdf)], prop)]

        with open(self.out_path, 'a+', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, self.header_list)
            writer.writerows(dt)

        with open(self.wrong_path, 'a+', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, self.wrong_header)
            writer.writerows(self.wrong_dt)
        self.wrong_dt.clear()

    def query_name(self, name, cid_list):
        if len(name) < 1:
            return False

        api = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/concepts/name/JSON?name='
        res = requests.get(f'{api}{quote_plus(name)}').json()
        if self.query_time % 17 == 0:
            time.sleep(random.uniform(1.1,2.2))
        self.query_time += 1
        if 'Fault' in res:
            return False
        try:
            for cid in res['ConceptsAndCIDs']['CID']:
                cid = str(cid)
                if not self.cid_set.__contains__(cid):
                    cid_list.append(cid)
                    self.cid_set.add(cid)
        except KeyError:
            with open(key_error_file, 'a+', encoding='utf-8') as f:
                f.write(name + '\n')

        return True

    def write_log(self, ss):
        with open(self.log_path, 'a+', encoding='utf-8') as f:
            f.write(ss + '\n')

    def __main__(self):
        # get information of 260 cid numbers per time
        # limit = 260

        now = 0
        cid_list = []
        self.write_log(f'[{time.ctime()}] : ----------start downloading----------')
        for name in self.name_dt.keys():
            # get valid name
            ret = self.query_name(name, cid_list)

            if not ret:              # wrong in query
                self.wrong_dt.append({'Name': name, 'Times': self.name_dt[name]})

                scl = [-1]  # special characters list
                for i in range(len(name)):
                    if not (name[i].isalpha() or name[i].isdigit()):
                        scl.append(i)
                scl.append(len(name))
                l_scl = len(scl)
                if 2 < l_scl < 8:
                    ss_list = [[]]        # sub str list
                    # ss_list[gap][i] stores all situations while start from i, gap

                    # init gap = 1
                    ss_gap_list = []
                    for i in range(l_scl - 1):
                        self.query_name(name[(scl[i]+1):scl[i+1]], cid_list)
                        ss_gap_list.append([name[(scl[i]+1):scl[i+1]]])
                    ss_list.append(ss_gap_list)

                    for gap in range(2, l_scl - 1):
                        ss_gap_list = []
                        for i in range(l_scl - 1 - gap):
                            ss_gap_i_list = []
                            for j in range(i + 1, i + gap):
                                list_j_ipg = ss_list[i + gap - j][j]
                                s_front = name[(scl[i]+1):scl[j]]
                                for s_rear in list_j_ipg:
                                    self.query_name(s_front + s_rear, cid_list)
                                    ss_gap_i_list.append(s_front + s_rear)
                            # including full substr
                            self.query_name(name[(scl[i]+1):scl[i+gap]], cid_list)
                            ss_gap_i_list.append(name[(scl[i]+1):scl[i+gap]])

                            # update ss_gap_list
                            ss_gap_list.append(ss_gap_i_list)
                        ss_list.append(ss_gap_list)
                elif l_scl >= 8:
                    with open(manually_query, 'a+', encoding='utf-8') as f:
                        f.write(name + '\n')
            # always waste a lot of name that queried
            # elif len(cid_list) >= limit:
            #     self.download(cid_list)
            #     cid_list = []

            if now % 11 == 0:
                self.write_log(f'[{time.ctime()}] : ---number {now}, name {name} just finished---')
            elif now % 50 == 0:
                print(f'[{time.ctime()}] : number {now}, name {name} just finished')
                self.download(cid_list)
                cid_list = []
            now += 1

        if len(cid_list) > 0:
            self.download(cid_list)


def call(key, fdir):
    MofsCrawl(in_file, key, f'{fdir}properties.csv', f'{fdir}sdf{seperator}',
              f'{fdir}wrong.csv', f'{fdir}log.txt').__main__()


if __name__ == '__main__':
    # S_PRECURSOR, M_PRECURSOR, O_PRECURSOR

    key_error_file = 'key_error.txt'
    # in windows
    # seperator = '\\'
    # in linux
    seperator = '/'
    # in_file = 'try.json'
    in_file = 'precursor.json'

    # f_dir = 'S_precursor' + seperator
    # manually_query = f_dir + 'man_que.txt'
    # call('S_precursor', f_dir)
    # print('-----------S_precursor done-----------')

#    f_dir = 'M_precursor' + seperator
#    manually_query = f_dir + 'man_que.txt'
#    call('M_precursor', f_dir)
#    print('-----------M_precursor done-----------')

    f_dir = 'O_precursor' + seperator
    manually_query = f_dir + 'man_que.txt'
    call('O_precursor', f_dir)
    print('-----------O_precursor done-----------')
