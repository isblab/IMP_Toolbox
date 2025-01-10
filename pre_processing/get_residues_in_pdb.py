######### Obtain all the residues present in a PDB file. #########
########## ------>"May the Force serve u well..." <------#########
##################################################################

############# One above all #############
##-------------------------------------##
import sys
import requests
import wget

import Bio
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import MMCIFParser

file_name = sys.argv[1]         # PDB/CIF file path.

# Run the script using
# python get_residues_in_pdb.py *.cif True

def get_confident_regions( file ):
	# Obtain all residues present in the PDB.CIF file.
	# file --> PDB or MMCIF file.
	residues = {}
	if file.split( "." )[1] == "cif":
		models = MMCIFParser().get_structure('cif', file)
	elif file.split( "." )[1] == "pdb":
		models = PDBParser().get_structure('pdb', file)
	else:
		print( "Incorrect file format..." )
		exit()

	for model in models:
		for chain in model:
			residues[chain.id] = []
			for residue in chain:
					if residue.id[0] == " ":
						if residue["CA"]:
							residues[chain.id].append( residue.id[1] )
		break
	return residues


def ranges( positions ):
	# get patches of continous regions.
	# positions --> list residues in pdb.

    nums = sorted( set( positions ) )
    # Get start, end positions for each continous confident patch.
    gaps = [[x, y] for x, y in zip( positions, positions[1:] ) if x+1 < y]
    edges = iter( positions[:1] + sum( gaps, [] ) + positions[-1:] )
    return list( zip( edges, edges ) )

response = wget.download( f"https://files.rcsb.org/download/{file_name}.pdb", f"{file_name}.pdb" )

file = f"./{file_name}.pdb"
all_residues = get_confident_regions( file )
print( f"Residues in PDB file: {file_name} \n" )
[print(  chain, "\t", ranges( all_residues[chain] ), "\n" ) for chain in all_residues.keys()]
