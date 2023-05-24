import argparse
from multiprocessing import *
import numpy as np
from matplotlib import pyplot as plt
import IMP
import RMF
import IMP.rmf
import os, sys

################################################################################
#################################### Inputs ####################################
################################################################################

num_copies = {"H2A": 2, "H2B": 2, "H3": 2, "H4": 2, "SPIN1": 1, "WDR76": 1}


################################################################################
################################## Functions ###################################
################################################################################

__doc__ = "Given a cluster of models and validation crosslinks file, find the minimum distance for each crosslink in the cluster and classify the crosslinks as satisfied or violated."


def parse_args():
    parser = argparse.ArgumentParser(
        description="Given a cluster of models and validation crosslinks file, find the minimum distance for each crosslink in the cluster and classify the crosslinks as satisfied or violated."
    )
    parser.add_argument(
        "--inputA",
        "-ia",
        dest="ia",
        help="cluster list of sample A RMFs. (eg. cluster.0.sample_A.txt)",
        required=True,
    )
    parser.add_argument(
        "--inputB",
        "-ib",
        dest="ib",
        help="cluster list of sample B RMFs. (eg. cluster.0.sample_B.txt)",
        required=True,
    )
    parser.add_argument(
        "--rmfA",
        "-ra",
        dest="ra",
        help="rmfA file. (eg. A_gsm_clust0.rmf3)",
        required=True,
    )
    parser.add_argument(
        "--rmfB",
        "-rb",
        dest="rb",
        help="rmfB file. (eg. B_gsm_clust0.rmf3)",
        required=True,
    )
    parser.add_argument(
        "--clmdl",
        "-c",
        dest="cmdl",
        help="Cluster center model rmf3 file. (eg. cluster_center_model.rmf3)",
        required=True,
    )
    parser.add_argument(
        "--textA",
        "-ta",
        dest="ta",
        help="Text file associated with rmfA. (eg. A_gsm_clust0.txt)",
        required=True,
    )
    parser.add_argument(
        "--xl",
        "-x",
        dest="xlfile",
        help="Validation crosslinks file. (eg. input/mhr/xlms/mhr_bs3dss_validation.dat)",
        required=True,
    )
    parser.add_argument(
        "--threshold",
        "-t",
        dest="threshold",
        help="Distance threshold for the given crosslink type",
        required=True,
    )

    return parser.parse_args()


def get_xls_from_file(input_file):
    xls = []
    with open(input_file, "r") as inf:
        for ln in inf.readlines():
            xls.append(ln.strip())
    return xls


def get_nmodels_in_A(ta_file):
    with open(ta_file, "r") as taf:
        ln_count = 0
        for ln in taf.readlines():
            ln_count += 1
    return ln_count


def get_xl_min_distances(hier, mdl, xl_list, copy_nums, xl_min_dist_dict):
    for xl in xl_list:
        if not xl.startswith("Protein1") and not xl.startswith("Linker"):
            prot1, res1, prot2, res2 = (
                xl.split(",")[0],
                int(xl.split(",")[1]),
                xl.split(",")[2],
                int(xl.split(",")[3]),
            )

            if (prot1 not in copy_nums) or (prot2 not in copy_nums):
                continue
            num_cp1 = copy_nums[prot1]
            num_cp2 = copy_nums[prot2]

            for cp1 in range(num_cp1):
                for cp2 in range(num_cp2):
                    # xlink = f"{prot1}.{cp1},{res1},{prot2}.{cp2},{res2}"
                    # print(prot1, "\t", prot2, "\t", cp1, "\t", cp2)
                    # try:
                    particle1 = IMP.atom.Selection(
                        hier, molecule=prot1, copy_index=cp1, residue_index=res1
                    ).get_selected_particle_indexes()[0]
                    particle2 = IMP.atom.Selection(
                        hier, molecule=prot2, copy_index=cp2, residue_index=res2
                    ).get_selected_particle_indexes()[0]
                    # except Exception as e:
                    #     xlink = f"{prot1}.{cp1},{res1},{prot2}.{cp2},{res2}"
                    #     print(xlink)
                    #     print(IMP.atom.Selection(hier, molecule=prot1, copy_index=cp1, residue_index=res1))
                    #     print(IMP.atom.Selection(hier, molecule=prot2, copy_index=cp2, residue_index=res1))
                    #     exit()
                    coords1 = IMP.core.XYZ(mdl, particle1)
                    coords2 = IMP.core.XYZ(mdl, particle2)
                    dist = IMP.core.get_distance(coords1, coords2)

                    if xl not in xl_min_dist_dict.keys():
                        xl_min_dist_dict[xl] = dist
                    else:
                        if dist < xl_min_dist_dict[xl]:
                            xl_min_dist_dict[xl] = dist


################################################################################
##################################### Main #####################################
################################################################################

args = parse_args()

val_xls = get_xls_from_file(args.xlfile)
nA = get_nmodels_in_A(args.ta)

# Create list of model indices for sampleA, sample_B
sample_A_models = []
sample_B_models = []

with open(args.ia, "r") as iaf:
    for ln in iaf.readlines():
        sample_A_models.append(int(ln.strip()))
with open(args.ib, "r") as ibf:
    for ln in ibf.readlines():
        sample_B_models.append(int(ln.strip()))

sample_A_models.sort()
sample_B_models.sort()

nModels = len(sample_A_models) + len(sample_B_models)
print(f"Total number of models:\t {nModels}")

mdl_a = IMP.Model()
rmf_fh_a = RMF.open_rmf_file_read_only(args.ra)
h_a = IMP.rmf.create_hierarchies(rmf_fh_a, mdl_a)[0]

mdl_b = IMP.Model()
rmf_fh_b = RMF.open_rmf_file_read_only(args.rb)
h_b = IMP.rmf.create_hierarchies(rmf_fh_b, mdl_b)[0]


# ---------------------------------------------------------------------------
# Running get_xl_min_distances() function using 2 parallel processes
# ---------------------------------------------------------------------------

max_models = max(len(sample_A_models), len(sample_B_models))
mp_manager = Manager()
xl_min_dist_a = mp_manager.dict()
xl_min_dist_b = mp_manager.dict()

"""
The two parallel processes read models from sampleA and sampleB independently and add the minimum distance crosslink entries
to two separate dictionaries. Thus, the processes never interact with each other. mp_manager.dict() manages the global variable
modifications carried out by these function (in our case these global variable modifications are the entries to min_dist dictionaries).
"""

for i in range(max_models):
    jobs = []
    statement_p1 = "Process-1: Not active"
    statement_p2 = "Process-2: Not active"

    if i < len(sample_A_models):
        sa = sample_A_models[i]
        IMP.rmf.load_frame(rmf_fh_a, int(sa))
        mdl_a.update()
        jobs.append(
            Process(
                target=get_xl_min_distances,
                args=(h_a, mdl_a, val_xls, num_copies, xl_min_dist_a),
            )
        )
        statement_p1 = (
            f"Process-1: Working on {i} out of {len(sample_A_models)} on SampleA"
        )

    if i < len(sample_B_models):
        sb = sample_B_models[i]
        fname = sb - nA
        IMP.rmf.load_frame(rmf_fh_b, fname)
        mdl_b.update()
        jobs.append(
            Process(
                target=get_xl_min_distances,
                args=(h_b, mdl_b, val_xls, num_copies, xl_min_dist_b),
            )
        )
        statement_p2 = (
            f"Process-2: Working on {i} out of {len(sample_B_models)} on SampleB"
        )

    # We will start the processes now:
    for job in jobs:
        job.start()
    print(statement_p1, "\t", statement_p2)

    # job.join() will halt the code execution till the processes finishes
    for job in jobs:
        job.join()

overall_xl_min_dist = {}
violated_crosslinks = []
with open(
    f"xlviol_validation_set_{(args.xlfile).split('/')[-1].split('.')[0]}.txt", "w"
) as outf:
    for key in xl_min_dist_a.keys():
        overall_xl_min_dist[key] = min(
            float(xl_min_dist_a[key]), float(xl_min_dist_b[key])
        )
        overall_min_dist = min(float(xl_min_dist_a[key]), float(xl_min_dist_b[key]))
        out_string = f"{key} : {overall_min_dist}\n"
        outf.write(out_string)
        if overall_min_dist > float(args.threshold):
            violated_crosslinks.append(key)

n_violated = len(violated_crosslinks)
percent_violated = (n_violated / len(val_xls)) * 100

with open(
    f"xlviol_validation_set_{(args.xlfile).split('/')[-1].split('.')[0]}.log", "w"
) as logf:
    log_str = f"Total number of crosslinks in the validation set are: {len(val_xls)}\nNumber of violated crosslinks from the validation set are: {n_violated}\nPercent \
violation from the validation set is: {percent_violated}"
    logf.write(log_str)
"""
with open('all_xl_details.txt', 'w') as outf:
    for entry in overall_xl_min_dist:
        outstr = f"{entry} : {overall_xl_min_dist[entry]}\n"
        outf.write(outstr)
"""
print("May the Force b with u...")
