import os 
import json
from collections import defaultdict


name_ids = defaultdict(list)

for _dir in os.listdir('.'):
    if _dir.endswith('_precursor'):
        file = _dir + '/name_id.json'
        with open(file, 'r') as f:
            name_id = json.load(f)

            for t in name_id:
                for item in name_id[t]:
                    item = int(item)
                    if item not in name_ids[t]:
                        name_ids[t].append(item)

print(len(name_ids))

with open('name_id.json', 'w') as f:
    json.dump(name_ids, f, ensure_ascii=False)
        