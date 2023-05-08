import os, sys
import IMP
import RMF
import IMP.core
import IMP.rmf

from Bio.PDB import *


###################################################################################################
############################################# Inputs ##############################################
###################################################################################################

input_file = sys.argv[1]

pdb_files = [
    # "/home/shreyas/Projects/washburn/wdr_work/v3/modeling/data/pdbs/4h75.pdb",
    "/home/shreyas/Projects/washburn/wdr_work/v3/modeling/data/pdbs/5gt0.pdb",
    # "/home/shreyas/Projects/washburn/wdr_work/v3/modeling/data/pdbs/AF-Q9H967-F1-model_v4.pdb",
    # "/home/shreyas/Projects/washburn/wdr_work/v3/modeling/data/pdbs/AF-Q9H967-F1-model_v4.pdb",
]

# The all proteins list has the following architecture:
# [{protein:{chain_id: [copy_number, residue range]}}, {protein:{chain_id: [copy_number, residue range]}]
# The order of entries in the offset list must be the same as that in the pdb_files list

all_proteins = [
    {
        "H2A": {"C": [0, range(14, 119)], "G": [1, range(13, 118)]},
        "H2B": {"D": [0, range(28, 126)], "H": [1, range(32, 126)]},
        "H3": {"A": [0, range(37, 136)], "E": [1, range(37, 136)]},
        "H4": {"B": [0, range(24, 103)], "F": [1, range(17, 103)]},
    },
    # {"SPIN1": {"A": [0, range(45, 259)]}},
    # {"WDR76": {"A": [0, range(245, 276)]}},
    # {"WDR76": {"A": [0, range(288, 626)]}},
]


# The all proteins list has the following architecture:
# [{protein:{chain_id,residue range}}, {protein:{chain_id,residue range}]
# The order of entries in the offset list must be the same as that in the pdb_files list


###################################################################################################
##################################### Get transformations #########################################
###################################################################################################

for file_index in range(len(pdb_files)):
    print(f"Aligning: {pdb_files[file_index]}")
    pdb_file = pdb_files[file_index]
    proteins = all_proteins[file_index]
    # print(proteins)
    ccm_mdl = IMP.Model()
    ccm = RMF.open_rmf_file_read_only(input_file)
    hier = IMP.rmf.create_hierarchies(ccm, ccm_mdl)[0]
    IMP.rmf.load_frame(ccm, 0)
    ccm_mdl.update()
    pdb_ca_mdl = IMP.Model()
    pdb_ca = IMP.atom.read_pdb(pdb_file, pdb_ca_mdl, IMP.atom.CAlphaPDBSelector())
    pdb_ca_mdl.update()

    new_mdl = IMP.Model()
    reload = IMP.atom.read_pdb(pdb_file, new_mdl)

    coords_pdb_ca = {}
    coords_ccm = {}

    for prot in proteins.keys():
        for chain_id in proteins[prot]:
            protein_name = prot

            sel_ca_pdb = IMP.atom.Selection(
                pdb_ca,
                resolution=1,
                chain_id=chain_id,
                residue_indexes=[i for i in proteins[prot][chain_id][1]],
            ).get_selected_particles()
            sel_ccm = IMP.atom.Selection(
                hier,
                resolution=1,
                molecule=protein_name,
                copy_index=proteins[prot][chain_id][0],
                residue_indexes=[i for i in proteins[prot][chain_id][1]],
            ).get_selected_particles()
            print(len(sel_ca_pdb), len(sel_ccm))
            print(
                protein_name, proteins[prot][chain_id][0], proteins[prot][chain_id][1]
            )
            # print(len(sel_ccm),len(sel_ca_pdb))
            # Remove coarse grained beads
            new_ccm_sel = []
            for selection in sel_ccm:
                if not IMP.atom.Fragment.get_is_setup(selection):
                    new_ccm_sel.append(selection)
                    # print(selection)

            # print(len(sel_ca_pdb),'\n\n', len(new_ccm_sel))
            # for i in range(len(new_ccm_sel)):
            # print(new_ccm_sel[i],sel_ca_pdb[i])

            coords_pdb_ca[protein_name] = [
                IMP.core.XYZ(i).get_coordinates() for i in sel_ca_pdb
            ]
            coords_ccm[protein_name] = [
                IMP.core.XYZ(i).get_coordinates() for i in new_ccm_sel
            ]
            # print(len(coords_pdb_ca[protein_name]),len(coords_ccm[protein_name]))
    _, transformation = IMP.pmi.analysis.Alignment(
        query=coords_pdb_ca, template=coords_ccm
    ).align()
    print(transformation)

    ###################################################################################################
    #################################### Transform and write PDB ######################################
    ###################################################################################################

    IMP.atom.transform(reload, transformation)
    IMP.atom.write_pdb(reload, f"./aligned_{pdb_file.split('/')[-1]}.pdb")
