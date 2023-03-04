######### Obtain high confidence regions from AF2 model. #########
########## ------>"May the Force serve u well..." <------#########
##################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import sys
import json

import Bio
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import MMCIFParser

cif = sys.argv[1]         # CIF file path.
pae = sys.argv[2]         # PAE json file path.

# Run the script using
# python get_high_confidence_region_from_AF2.py CHD4/CHD4.cif CHD4/CHD4_PAE.json

def get_confident_regions( cif, path ):
	# Obtain confident resgions from AF2 model.
	# A residue is confiently predicted if b-factor > 70.
	# cif --> AF2 MMCIF file.

	confident_residues = []
	models = MMCIFParser().get_structure('cif', cif)
	for model in models:
		for chain in model:
			for residue in chain:
				class ResidueSelect(Select):
				    def accept_residue(self, residue):
				        if residue["CA"] and residue["CA"].get_bfactor() > 70.0:
				            return 1
				        else:
				            return 0
				if residue["CA"] and residue["CA"].get_bfactor() > 70.0:
					confident_residues.append( residue.id[1] )
	io = PDBIO()
	io.set_structure(models)
	io.save( f"{path}_confident.pdb", ResidueSelect() )

	return confident_residues


def ranges( positions ):
	# get patches of continous confident regions.
	# positions --> list containng confident residues.

    nums = sorted( set( positions ) )
    # Get start, end positions for each continous confident patch.
    gaps = [[x, y] for x, y in zip( positions, positions[1:] ) if x+1 < y]
    edges = iter( positions[:1] + sum( gaps, [] ) + positions[-1:] )
    return list( zip( edges, edges ) )


def merge_confident_regions( confident_patches, pae ):
	# Merge confident residue patches if PAE <= 5.

	f = open( pae, "r" )
	data = json.loads( f.read() )
	x = np.array( data[0]["predicted_aligned_error"] )

	confident_regions, merged = [], []
	for patch1 in confident_patches:
		if patch1 not in merged:
			# merged = []
			r1, r2 = patch1[0], patch1[1]
			idx = confident_patches.index( patch1 )+1 if  patch1 != confident_patches[-1] else -1
			for patch2 in confident_patches[idx:]:
				c1, c2 = patch2[0], patch2[1]
				# Merge all confident patches with the interaction regions PAE <=5.
				if np.all( x[r1-1:r2, c1-1:c2] <= 5 ):
					# confident_regions.append( ( r1, c2 ) )
					merged.append( patch2 )
					patch1 = (r1, c2)
					r2 = c2
					continue
				# Skip if the adjacent patches are not merged.
				else:
					break
			confident_regions.append( patch1 )
	return confident_regions

path = cif.split( ".cif" )[0]
confident_patches = ranges( get_confident_regions( cif, path ) )
confident_regions = merge_confident_regions( confident_patches, pae )
included_in_topology = [patch for patch in confident_regions if patch[1]-patch[0] >= 20]
#TODO Step 3.
# Get the regions to print in the topo file
# min_patch = 20. 

print( f"All confident patches from pLDDT scores: \n {confident_patches} \n" )
print( f"Merged confident patches using PAE domain scores: \n {confident_regions} \n" )
print( f"Merged confident patches with atleast 20 residues: \n {included_in_topology}" )
with open(f"{path}_confident.txt", "w") as w:
	w.writelines("<<-----------Confident regions from AF2 model----------->>\n")
	w.writelines( "All confident patches from pLDDT scores: \n" )
	w.writelines( str( confident_patches ) + "\n\n" )
	w.writelines( "Merged confident patches using PAE domain scores: \n" )
	w.writelines( str( confident_regions ) + "\n\n" )
	w.writelines( "Considering only merged confident patches with atleast 20 residues: \n" )
	w.writelines( str( included_in_topology ) + "\n\n" )

print( "May the Force be with you..." )