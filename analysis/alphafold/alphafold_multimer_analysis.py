# To get the regions of interest in a given PDB with different PAE and PLDDT

import numpy as np
import sys
import pickle
from Bio.PDB import PDBParser
import json

path = sys.argv[1] # path to output of AF2 multimer
outf = sys.argv[2] # path to where output should be stored
name = sys.argv[3] # name of complex e.g. dp-pg


pdb = f'{path}/ranked_0.pdb'
with open(f'{path}/ranking_debug.json', 'r') as f:
    o = json.loads(f.read())
best_pkl = o['order'][0]
pkl = f'{path}/result_{best_pkl}.pkl'
with open(pkl, 'rb') as f:
    data = pickle.load(f)

models = PDBParser().get_structure('pdb', pdb)

f2 = open(f'{outf}/{name}_list_high_confidence.txt', 'w')

for model in models:

    chains = [c for c in model.get_chains()]

    len_chain_1 = len([r for r in chains[0]])

    for i,ra in enumerate(chains[0]):

        if ra['CA'].get_bfactor()<70:
            continue

        for j,rb in enumerate(chains[1]):

            if rb['CA'].get_bfactor()<70:
                continue

            if rb['CA'] - ra['CA'] < 10.0 and data['predicted_aligned_error'][i][len_chain_1+j-1]<5.0: # 10 A distance defining interface
                print("A",i+1,"B",j+1)
                f2.write(f'{i+1}\tA\t{j+1}\tB\n')

f2.close()
