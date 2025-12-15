############## Create contact maps for protein pairs ###############
########## ------>"May the Force serve u well..." <------###########
####################################################################


############# One above all #############
##-------------------------------------##
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import IMP
import RMF
import IMP.rmf
import os, sys
from multiprocessing import Process, Manager

__doc__ = "Get contact maps for a given protein pair"


def parse_args():
    parser = argparse.ArgumentParser(description="Get contact maps for a given protein pair")
    parser.add_argument('--inputA', '-ia', dest="ia", help='cluster list of sample A RMFs. (eg. cluster.0.sample_A.txt)', required=True)
    parser.add_argument('--inputB', '-ib', dest="ib", help='cluster list of sample B RMFs. (eg. cluster.0.sample_B.txt)', required=True)
    parser.add_argument('--rmfA', '-ra', dest="ra", help='rmfA file. (eg. A_gsm_clust0.rmf3)', required=True)
    parser.add_argument('--rmfB', '-rb', dest="rb", help='rmfB file. (eg. B_gsm_clust0.rmf3)', required=True)
    parser.add_argument('--clmdl', '-c', dest="cmdl", help='Cluster center model rmf3 file. (eg. cluster_center_model.rmf3)', required=True)
    parser.add_argument('--textA', '-ta', dest="ta", help='Text file associated with rmfA. (eg. A_gsm_clust0.txt)', required=True)
    parser.add_argument('--prots', '-p', dest="prots", help='Proteins (eg. MTA1.0,MTA1.1)', required=True)
    # parser.add_argument('--size', '-s', dest="size", help='Size of each protein', required=True)
    # parser.add_argument('--rb_dimer', '-rbd', dest="rbd", help='Residue range in the rigid body dimer.', required=True)
    # parser.add_argument('--prot2', '-p2', dest="prot2", help='Protein2 (eg. MTA1.1)', required=True)
    return parser.parse_args()


def get_bead_name(particle):

    '''
    Input: particle
    Output: bead name in the format molecule_name:copy_number:start_residue:end_residue
    '''

    mol_name = IMP.atom.get_molecule_name(IMP.atom.Hierarchy(particle))
    copy_number=IMP.atom.get_copy_index(IMP.atom.Hierarchy(particle))

    if IMP.atom.Fragment.get_is_setup(particle): # CG bead
        residues_in_bead = IMP.atom.Fragment(particle).get_residue_indexes()
        bead_name = mol_name+":"+str(copy_number)+":"+str(min(residues_in_bead))+":"+str(max(residues_in_bead))
    else:
        residue_in_bead = str(IMP.atom.Residue(particle).get_index())
        bead_name = mol_name+":"+str(copy_number)+":"+residue_in_bead+":"+residue_in_bead

    return bead_name


def get_nmodels_in_A(ta_file):
    with open(ta_file,'r') as taf:
        ln_count = 0
        for ln in taf.readlines():
            ln_count += 1
    return ln_count


def get_distances(sample,rmf_fh,mdl,hier,distances_dict):
    i = 0
    for sa in sample:
        IMP.rmf.load_frame(rmf_fh, int(sa))
        mdl.update()
        sel1 = IMP.atom.Selection(hier, resolution=RESOLUTION, molecule=p1name, copy_index=int(copy1))   
        sel2 = IMP.atom.Selection(hier, resolution=RESOLUTION, molecule=p2name, copy_index=int(copy2)) #, residue_indexes = range(1,180)

        print(f'Reading frame {i} out of {len(sample)}')
        i+=1

        for bead1 in sel1.get_selected_particles():
            for bead2 in sel2.get_selected_particles():
                # print(bead1,bead2,'##############')
                dist = IMP.core.get_distance(IMP.core.XYZR(bead1),IMP.core.XYZR(bead2))
                if dist<0:
                    dist=0
                key = get_bead_name(bead1) + "--" + get_bead_name(bead2)
                if key in distances_dict.keys():
                    distances_dict[key] += dist
                else:
                    distances_dict[key] = dist


def get_interacting_beads( dist_dict, threshold = 20 ):
    # Given a distance matrix, return all the bead pairs with distance < threshold.

    interaction_dict = {bead.split("--")[0]:[] for bead in dist_dict.keys()}

    for k in dist_dict.keys():
        if dist_dict[k] <= threshold:
            bead1, bead2 = k.split("--")
            p1, c1, r11, r12 = bead1.split( ":" )
            p2, c2, r21, r22 = bead2.split( ":" )
            r11, r12, r21, r22 = int( r11 ), int( r12 ), int( r21 ), int( r22 )
            if f"{p1}-{p2}" in proteins.keys():
                if not get_intersection( proteins[f"{p1}-{p2}"][p1], ( r11, r12 ) ) or not get_intersection( proteins[f"{p1}-{p2}"][p2], ( r21, r22 ) ):
                    interaction_dict[bead1].append( bead2 )
            
            elif f"{p2}-{p1}" in proteins.keys():
                if not get_intersection( proteins[f"{p2}-{p1}"][p1], ( r11, r12 ) ) or not get_intersection( proteins[f"{p2}-{p1}"][p2], ( r21, r22 ) ):
                    interaction_dict[bead1].append( bead2 )
            
            else:
                interaction_dict[bead1].append( bead2 )
    return interaction_dict


def convert_beads_to_patches( interacting_beads ):
    # Given a list of interacting beads, create a list of interacting patches.

    prot1 = ""
    patches = []
    for bead in interacting_beads:
        prot1, cp1, res1, res2 = bead.split( ":" )
        patches.append( (res1, res2) )

    return prot1, patches


def merge_adjacent_patches( patches ):
    # Merge adjacent patches.

    adjacent_patches, merged = [], []
    for patch1 in patches:
        if patch1 not in merged:
            r1, r2 = patch1[0], patch1[1]
            idx = patches.index( patch1 )+1 if  patch1 != patches[-1] else -1
            for patch2 in patches[idx:]:
                c1, c2 = patch2[0], patch2[1]
                # Merge patches if they are adjacent.
                if int( c1 ) == int( r2 )+1:
                    merged.append( patch2 )
                    patch1 = (r1, c2)
                    r2 = c2
                    continue
                # Skip if the adjacent patches are not merged.
                else:
                    break
            adjacent_patches.append( patch1 )
    return adjacent_patches

def get_intersection( x, y ):
    # Obtain the intersecting positions in the input lists.
    # x --> region modeled as a rigid body dimer.
    # y --> interacting bead.
    # x = [1,2,3,4]; y = [3,4,5,6] --> intersection --> [3,4]

    print( x, "\t\t", y)
    intersect = []
    flag = True
    for x_ in x:
        x_ = np.arange( x_[0], x_[1]+1 ) if x != [] else []
        y_ = np.arange( y[0], y[1]+1 )
        if x != [] and len( x_ ) > len( y_ ):
            intersect.append( set( x_ ).intersection( y_ ) )
        
        elif x != [] and len( x_ ) < len( y_ ):
            intersect.append( set( y_ ).intersection( x_ ) )
        
        else:
            intersect.append( set() )
            continue
    # False if iteracting beads are not part of a rigid body dimer else True.
    flag = False if all( [i == set() for i in intersect] ) else True

    return flag

################################################################################
############################# Main #############################################
################################################################################
if __name__ == '__main__':
    args = parse_args()
    # print('############################\n',args.prots,'\n############################')
    prot1,prot2 = (args.prots).split(',')
    nA = get_nmodels_in_A(args.ta)
    RESOLUTION = 10

    # Step 1. Go through each RMF
    # Step 2. In each RMF, get the two selections corresponding to the 2 proteins.
    #          Create a matrix of appropriate size, get the res ranges per bead in a list, if reading first RMF.
    # Step 3. For each bead pair, get the distance between the bead centers. Add it to the distances in the appropriate matrix entry.
    # Step 4. After 1-3, get the average for all matrix entries by dividing by num of models

    # Step 5. Fix the matrix for plotting.
    # Get the size and rigid body dimers dict.
    # size of protein = protein_length + 1.
    sizes_dict = {"HDAC1":483, "SAP30": 221, "SUDS3":329, "SIN3A":1274}
    proteins = {
    "SAP30-SIN3A": {
    "SAP30": [(130,219)],
    "SIN3A": [(456,527)],
    }, 
    "SUDS3-SIN3A": {
    "SUDS3": [(205,228)],
    "SIN3A": [(604,728)],
    }
    }
    
    # Create list of model indices for sampleA, sample_B
    sample_A_models = []
    sample_B_models = []

    with open(args.ia,'r') as iaf:
         for ln in iaf.readlines():
             sample_A_models.append(int(ln.strip()))
    with open(args.ib,'r') as ibf:
         for ln in ibf.readlines():
             sample_B_models.append(int(ln.strip()))

    sample_A_models.sort()
    sample_B_models.sort()

    sample_B_models_id_corrected = []
    for bmodel in sample_B_models:
        sample_B_models_id_corrected.append(bmodel-nA)

    # print(sample_B_models_id_corrected)
    # print(sample_B_models)
    nModels = len(sample_A_models)+len(sample_B_models)
    print(f'Total number of models:\t {nModels}')


    # Step 1. Go through each RMF
    sample_A_rmf = args.ra
    sample_B_rmf = args.rb


    mdl_a = IMP.Model()
    rmf_fh_a = RMF.open_rmf_file_read_only(sample_A_rmf)
    h_a = IMP.rmf.create_hierarchies(rmf_fh_a, mdl_a)[0]

    mdl_b = IMP.Model()
    rmf_fh_b = RMF.open_rmf_file_read_only(sample_B_rmf)
    h_b = IMP.rmf.create_hierarchies(rmf_fh_b, mdl_b)[0]


    p1name = prot1.split('.')[0]
    p2name = prot2.split('.')[0]

    copy1 = prot1.split('.')[1]
    copy2 = prot2.split('.')[1]

    size1 = sizes_dict[prot1.split('.')[0]]
    size2 = sizes_dict[prot2.split('.')[0]]
    # print(prot1, copy1, prot2, copy2)

    print(prot1,'\t',prot2)
    
    manager = Manager()
    distances = manager.dict()
    

    i = 0

    #Creates an object of class Process. Specify the target (function to be executed) and args (arguments to be passed to the function).
    sample_a_process = Process(target=get_distances, args=(sample_A_models, rmf_fh_a, mdl_a, h_a, distances))
    sample_b_process = Process(target=get_distances, args=(sample_B_models_id_corrected, rmf_fh_b, mdl_b, h_b, distances))             


    sample_a_process.start()
    print('Reading rmfA file')
    sample_b_process.start()
    print('Reading rmfB file')

    # get_distances(sample_B_models, rmf_fh_b, mdl_b, h_b, distances) 
    sample_a_process.join()
    sample_b_process.join()


    for k in distances.keys():
        distances[k] = distances[k]/nModels
    # nBeadsA = len(sel1.get_selected_particles())+1
    # nBeadsB = len(sel2.get_selected_particles())+1

    print('Distance Matrix Size: ',size1,size2,'\n')
    mat = np.zeros((size1,size2))

    for k in distances.keys():
        # print(k)
        bead1 = k.split('--')[0]
        bead2 = k.split('--')[1]
        # print(bead2)
        res_start_1 = int(bead1.split(':')[2])
        res_end_1 =   int(bead1.split(':')[3])
        res_start_2 = int(bead2.split(':')[2])
        res_end_2 = int(bead2.split(':')[3])
        
        for i in range(res_start_1,res_end_1+1):
            for j in range(res_start_2,res_end_2+1):
                mat[i,j] = distances[k]

    np.savetxt(f'{p1name}.{copy1}-{p2name}.{copy2}_Distance-matrix.csv', mat, delimiter = ',')

    
    for threshold in [10, 20]:
        interacting_beads = get_interacting_beads( distances, threshold )
        interactions = {}
        for bead1 in interacting_beads.keys():
            p2, patches = convert_beads_to_patches( interacting_beads[bead1] )
            merged_patches = merge_adjacent_patches( patches )
            interactions[bead1] = merged_patches

        prot1_dict = {f"{p1name}.{copy1}": []}
        prot2_dict = {f"{p2name}.{copy2}": []}

        for key in interactions.keys():
            if interactions[key] != []:
                for inter in interactions[key]:
                    p1, cp1, r1, r2 = key.split( ":" )
                    prot1_dict[f"{p1name}.{copy1}"].append( f"{r1}-{r2}" )
                    prot2_dict[f"{p2name}.{copy2}"].append( f"{inter[0]}-{inter[1]}" )
        
        # Save interacting beads to a csv file.
        df = pd.DataFrame()
        if f"{p1name}.{copy1}" == f"{p2name}.{copy2}":
            df[f"{p1name}.{copy1}"] = prot1_dict[f"{p1name}.{copy1}"]
            df[f"{p2name}.{copy2}_"] = prot2_dict[f"{p2name}.{copy2}"]
        else:
            df[f"{p1name}.{copy1}"] = prot1_dict[f"{p1name}.{copy1}"]
            df[f"{p2name}.{copy2}"] = prot2_dict[f"{p2name}.{copy2}"]     
        df.to_csv( f"{p1name}.{copy1}-{p2name}.{copy2}_Contacts_{threshold}.csv", index = False )


    plt.figure( dpi = 600 )
    plt.imshow(mat, cmap='hot')
    plt.xlabel(prot2)
    plt.ylabel(prot1)
    plt.colorbar()
    plt.savefig(f'{p1name}.{copy1}-{p2name}.{copy2}_contact-map.png',dpi=600)

    # plt.show()
print("May the Force serve u well!!!")
exit(0)
