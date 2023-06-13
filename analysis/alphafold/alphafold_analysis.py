# To get the regions of interest in a given PDB with different PAE and PLDDT

import numpy as np
import matplotlib.gridspec as gridspec
import sys
import pickle
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
import json

path = sys.argv[1]
outf = sys.argv[2]
name = sys.argv[3]


def collate_in_order(x):
    final = []
    prev = x[0]
    curr = [x[0]]
    for i in x[1:]:
        if i == (prev + 1):
            curr.append(i)
        elif len(curr) > 1:
            final.append((curr[0], curr[-1]))
            curr = [i]
        else:
            final.append(curr[0])
            curr = [i]
        prev = i
    if len(curr) > 1:
        final.append((curr[0], curr[-1]))
    else:
        final.append(curr[0])
    return final


def get_confident_regions(f):
    bfactors = []
    locations = []
    models = PDBParser().get_structure('pdb', f)
    nums = []
    for model in models:
        for chain in model:
            nums.append(0)
            for residue in chain:
                if residue['CA']:
                    nums[-1] += 1
                    bfactors.append(residue['CA'].get_bfactor())
                    locations.append(np.array(residue['CA'].get_vector().get_array()))
    return bfactors, nums, locations


pdb = f'{path}/ranked_0.pdb'
with open(f'{path}/ranking_debug.json', 'r') as f:
    o = json.loads(f.read())
best_pkl = o['order'][0]
pkl = f'{path}/result_{best_pkl}.pkl'
with open(pkl, 'rb') as f:
    data = pickle.load(f)

b, n, locs = get_confident_regions(pdb)
locs = np.array(locs)
fig = plt.figure(figsize=(8, 8))
gs = gridspec.GridSpec(10, 11)
ax_main = plt.subplot(gs[2:, :8])
ax_xDist = plt.subplot(gs[:2, :8])
ax_yDist = plt.subplot(gs[2:, 8:10])
only_cbar = plt.subplot(gs[:, 10])
im = ax_main.imshow(data['predicted_aligned_error'][:n[0], n[0]:], aspect='auto', vmin=0, vmax=15, cmap='magma')
ax_main.set_xticks(list(range(0, n[1], 50)))
ax_main.set_yticks(list(range(0, n[0], 50)))
plt.colorbar(im, cax=only_cbar, use_gridspec=True)
new_b = [x if x > 70 else 0 for x in b]
ax_xDist.fill_between(list(range(n[1])), b[n[0]:], [0] * n[1], color='black', zorder=0)
ax_xDist.fill_between(list(range(n[1])), new_b[n[0]:], [0] * n[1], color='red', zorder=1)
ax_yDist.fill_betweenx(list(range(n[0])), b[:n[0]], [0] * n[0], color='black', zorder=0)
ax_yDist.fill_betweenx(list(range(n[0])), new_b[:n[0]], [0] * n[0], color='red', zorder=1)
ax_xDist.set_xticks([])
ax_xDist.set_xlim(ax_main.get_xlim())
ax_yDist.set_yticks([])
ax_yDist.set_ylim(ax_main.get_ylim())
plt.savefig(f'{outf}/{name}_pae_contact_map.png')
plt.close()

fig, ax = plt.subplots(1, 3, figsize=(24, 8))
ax_main = ax[0]
ax_1 = ax[1]
ax_2 = ax[2]
pae_matrix = (data['predicted_aligned_error'][:n[0], n[0]:] < 5)
plddt_matrix = np.zeros(pae_matrix.shape, dtype=bool)
for i in range(plddt_matrix.shape[0]):
    plddt_matrix[i] = b[i] > 70
for i in range(plddt_matrix.shape[1]):
    plddt_matrix[:, i] = b[n[0] + i] > 70
contact_matrix = np.zeros(pae_matrix.shape, dtype=float)
for i in range(n[0]):
    for j in range(n[1]):
        l1 = locs[i]
        l2 = locs[n[0] + j]
        contact_matrix[i, j] = np.sqrt(np.sum((l1 - l2) ** 2))
boolean = np.logical_and(plddt_matrix, pae_matrix)
im = ax_main.imshow(contact_matrix, aspect='auto', vmin=10, vmax=30, cmap='magma')
ax_1.imshow(boolean, aspect='auto', cmap='magma')
censored_matrix = contact_matrix.copy()
censored_matrix[np.logical_not(boolean)] = 100
ax_2.imshow(censored_matrix, aspect='auto', vmin=10, vmax=30, cmap='magma')
plt.colorbar(im, ax=ax_main)
ax_main.set_xticks(list(range(0, n[1], 50)))
ax_main.set_yticks(list(range(0, n[0], 50)))
ax_1.set_xticks(list(range(0, n[1], 50)))
ax_1.set_yticks(list(range(0, n[0], 50)))
ax_2.set_xticks(list(range(0, n[1], 50)))
ax_2.set_yticks(list(range(0, n[0], 50)))
plt.tight_layout()
plt.savefig(f'{outf}/{name}_pae_contact_map_binary.svg')
plt.close()

b = np.array(b)
chains = []
resnums = []
for i in range(len(n)):
    chains += [i] * n[i]
    resnums += list(range(1, n[i] + 1))
chains = np.array(chains)
resnums = np.array(resnums)
assert (np.max(np.abs(np.round(np.array(data['plddt']), 2) - b)) < 1e-3)

# f2 = open('compare.txt', 'w')
f2 = open(f'{outf}/{name}_list_high_confidence.txt', 'w')
f3 = open(f'{outf}/{name}_list_high_confidence_contacts.txt', 'w')
with open(f'{outf}/{name}_list.txt', 'w') as f:
    f.write('Residues\tChain\talign_region')
    f2.write('Residues\tChain\talign_region')
    f3.write('Residues\tChain\talign_region')
    for i in range(len(b)):
        b1 = b[i] > 70
        if b1:
            b2 = data['predicted_aligned_error'][i] <= 5
            resnum_sub = resnums[b2]
            chain_sub = chains[b2]
            resnum_sub = [resnum_sub[j] for j in range(len(resnum_sub)) if chain_sub[j] != chains[i]]
            b_sub = [b[b2][j] for j in range(len(b[b2])) if chain_sub[j] != chains[i]]
            loc_sub = [locs[b2][j] for j in range(len(locs[b2])) if chain_sub[j] != chains[i]]
            loc_sub = [np.sqrt(np.sum((j - locs[i]) ** 2)) for j in loc_sub]
            loc_sub = np.array(loc_sub)
            b_sub = np.array(b_sub)
            if len(resnum_sub) > 0:
                f.write(f'\n{resnums[i]}\t{chains[i]}\t{collate_in_order(resnum_sub)}')
            #                 f2.write(f'\n{resnums[i]}\t{chains[i]}\t{resnum_sub}')
            if np.sum(b_sub > 70) > 0:
                f2.write(f'\n{resnums[i]}\t{chains[i]}\t{collate_in_order(np.array(resnum_sub)[b_sub > 70])}')
            both_plddt_and_contact = np.logical_and(b_sub > 70, loc_sub <= 10)
            if np.sum(both_plddt_and_contact) > 0:
                f3.write(f'\n{resnums[i]}\t{chains[i]}\t{collate_in_order(np.array(resnum_sub)[both_plddt_and_contact])}')
f2.close()
f3.close()