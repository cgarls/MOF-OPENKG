import matplotlib.pyplot as plt
import json

with open('data/precursor_info.json', 'r') as f:
    dic = json.load(f)

fig, ax = plt.subplots(3, 1, figsize=(40, 20))
plt.xticks(rotation=270)
fig.subplots_adjust(hspace=0.4, wspace=0.2)

for ix, t in enumerate(['M_precursor', "S_precursor", "O_precursor"]):
    ax[ix].title.set_text(f'barchart {t}: {len(dic[t])}')
    ax[ix].bar(list(dic[t].keys())[:200], list(dic[t].values())[:200])
    for tick in ax[ix].get_xticklabels():
        tick.set_rotation(90)

# for ix, t in enumerate(['M_precursor', "S_precursor", "O_precursor"]):
#     ax[ix, 1].title.set_text(f'双对数 {t}: {len(dic[t])}')
#     ax[ix, 1].set_xscale('log')
#     ax[ix, 1].set_yscale('log')
#     ax[ix, 1].plot(list(dic[t].keys())[:200], list(dic[t].values())[:200])
#     for tick in ax[ix, 1].get_xticklabels():
#         tick.set_rotation(90)

plt.savefig('data/draw.png')