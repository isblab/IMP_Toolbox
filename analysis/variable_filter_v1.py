"""
Algorithm
Examine a set of multipliers (mean, mean-0.25 std, .. ) for set of data restraints 
Choose the largest number of models (lowest multiplier) 
    [such that nA and nB are each less than 15k: earlier ] about the same as : nA +nB < 20k
    It passes the KS test on total score and (A,B) looks similar on score distribution
    If you dont find a multiplier even after the narrower search in point 1, [take a random subset of the nearest multiplier that passes the KS test] OR take a single lenient cutoff on EV (less than mean) along with score multipliers on the other restraints. 

"""

import os, sys
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tabulate import tabulate

sys.path.append("/home/shreyas/imp-clean/imp/modules/sampcon/pyext/src/")
import scores_convergence


# This function will take the individual dataframes and compare them with the mean-multiplier*std of common dataframe and output the dataframe of the models that satisfy the filter
def variable_filter(std_multiplier, df):
    if args.evr:
        if args.use_gemr:
            selected_df = df[
                (
                    df["XLs_sum"]
                    <= common_df_mean["XLs_sum"]
                    + std_multiplier * common_df_std["XLs_sum"]
                )
                & (
                    df["GaussianEMRestraint_None"]
                    <= common_df_mean["GaussianEMRestraint_None"]
                    + std_multiplier * common_df_std["GaussianEMRestraint_None"]
                )
                & (df["EV_sum"] <= common_df_mean["EV_sum"])
            ]
        else:
            selected_df = df[
                (
                    df["XLs_sum"]
                    <= common_df_mean["XLs_sum"]
                    + std_multiplier * common_df_std["XLs_sum"]
                )
                & (df["EV_sum"] <= common_df_mean["EV_sum"])
            ]

    else:
        if args.use_gemr:
            selected_df = df[
                (
                    df["XLs_sum"]
                    <= common_df_mean["XLs_sum"]
                    + std_multiplier * common_df_std["XLs_sum"]
                )
                & (
                    df["GaussianEMRestraint_None"]
                    <= common_df_mean["GaussianEMRestraint_None"]
                    + std_multiplier * common_df_std["GaussianEMRestraint_None"]
                )
            ]
        else:
            selected_df = df[
                (
                    df["XLs_sum"]
                    <= common_df_mean["XLs_sum"]
                    + std_multiplier * common_df_std["XLs_sum"]
                )
            ]
    return selected_df


"""
Enter the names of data restraints as they are in the column headers of
selected_models_*_cluster*_detailed.csv files
"""
parser = argparse.ArgumentParser()

parser.add_argument(
    "-e",
    "--evr",
    action="store_true",
    default=False,
    help="Shall I use the Excluded Volume Restraint in the filter?",
)
parser.add_argument(
    "-gemr",
    "--use_gemr",
    action="store_true",
    default=False,
    help="Shall I use the Gaussian EM Restraint in the filter?",
)
parser.add_argument(
    "-c", "--cluster_num", default=0, help="On which cluster shall I run the filter?"
)
parser.add_argument(
    "-lc",
    "--lowest_cutoff",
    default=-2,
    help="What is the standard deviation multiplier for the most stringent cutoff?",
)
parser.add_argument(
    "-hc",
    "--highest_cutoff",
    default=10,
    help="What is the standard deviation multiplier for the most lenient cutoff?",
)
parser.add_argument(
    "-ss", "--step_size", default=0.1, help="What step size shall I use?"
)
parser.add_argument(
    "-n",
    "--num_models",
    default=30000,
    type=int,
    help="What is the maximum number of models you want?",
)
parser.add_argument(
    "-g",
    "--gsmsel",
    default=os.getcwd() + "/model_analysis/",
    type=str,
    help="Where is the model_analysis directory",
)


print("If you are using the EVR flag, consider reducing the step size")

args = parser.parse_args()

cluster_num = str(args.cluster_num)
lowest_cutoff = float(args.lowest_cutoff)
highest_cutoff = float(args.highest_cutoff)
step_size = float(args.step_size)
num_models = int(args.num_models)
gsm_sel_dir = str(args.gsmsel)
print(args)
sys_name = "ignore"

data_restraint_names = []

cluster_csv_fileA = (
    gsm_sel_dir + "selected_models_A_cluster" + cluster_num + "_detailed.csv"
)
cluster_csv_fileB = (
    gsm_sel_dir + "selected_models_B_cluster" + cluster_num + "_detailed.csv"
)


std_mult_dtrst = []
for i in np.arange(highest_cutoff, lowest_cutoff - step_size, -step_size):
    std_mult_dtrst.append(round(i, 3))
# print(std_mult_dtrst)


# Reading the CSV files and combining them to apply a common cutoff
dfA = pd.read_csv(cluster_csv_fileA)
dfB = pd.read_csv(cluster_csv_fileB)
print("Loaded the csv files")
df_list = [dfA, dfB]
common_df = pd.concat(df_list, ignore_index=True)
print("Calculating mean")
common_df_mean = common_df.mean()
print("Calculating standard deviation")
common_df_std = common_df.std()


print("Running variable filter")
# out = np.array([['DRest_Multiplier', 'nModelsA', 'nModelsB', 'KS_D-value', 'KS_p-value']])
out = []
out_str = f"Cluster number: {cluster_num} \nLowest cutoff: {lowest_cutoff} \nHighest cutoff: {highest_cutoff} \nStep size: {step_size} \nMaximum number of models to be selected: {num_models} \n\n"
mult_found = False

for multiplier in std_mult_dtrst:
    sel_dfA = variable_filter(multiplier, dfA)
    sel_dfB = variable_filter(multiplier, dfB)
    # Combining the score files for checking for run representation.
    sel_df_list = [sel_dfA, sel_dfB]
    sel_common_df = pd.concat(sel_df_list, ignore_index=True)
    nModelsT = len(sel_common_df.index)
    nModelsA = len(sel_dfA.index)
    nModelsB = len(sel_dfB.index)
    nRunsA = sel_dfA.traj.nunique()
    nRunsB = sel_dfB.traj.nunique()

    # Obtaining the scoresA and scoresB for sampling convergence
    scoresA = list(sel_dfA["Total_Score"])
    scoresB = list(sel_dfB["Total_Score"])
    scores = list(sel_common_df["Total_Score"])

    # Check if the two score distributions are similar

    if len(scoresA) == 0 or len(scoresB) == 0:
        continue
    ksd_pval = scores_convergence.get_scores_distributions_KS_Stats(
        scoresA, scoresB, 100, sys_name
    )
    ksd = ksd_pval[0]
    ksp = ksd_pval[1]
    # out = np.append(out,[multiplier, nModelsA, nModelsB, ksd, ksp])
    results = [multiplier, nModelsA, nModelsB, nRunsA, nRunsB, ksd, ksp]
    out.append(results)
    if nModelsA + nModelsB <= args.num_models:
        if ksp > 0.05 or (ksp <= 0.05 and ksd < 0.3):
            mult_found = True
            break

out = tabulate(
    out,
    headers=[
        "DRest_Multiplier",
        "nModelsA",
        "nModelsB",
        "nRunsA",
        "nRunsB",
        "KS_D-value",
        "KS_p-value",
    ],
)
print(out)

if mult_found == True:
    print(f"\nOptimal filter found.\nExtracted at {multiplier}")
    with open("var_filt_out.log", "w") as outf:
        outf.write(out_str)
        outf.write(out)
        outf.write(f"\n\nOptimal filter found.\nExtracted at {multiplier}")

    nBins = int(max(scores) - min(scores))
    plt.figure()
    plt.hist(scoresA, bins=nBins, histtype="step", label="ScoresA")
    plt.hist(scoresB, bins=nBins, histtype="step", label="ScoresB")
    plt.title("Scores of sampleA and sampleB")
    plt.xlabel("Total Score")
    plt.ylabel("nModels")
    plt.legend()
    plt.savefig("var_filt_out.png")
    # plt.show()
    sel_dfA.to_csv(cluster_csv_fileA.replace("selected_models", "good_scoring_models"))
    sel_dfB.to_csv(cluster_csv_fileB.replace("selected_models", "good_scoring_models"))

else:
    print("\nOptimal multiplier not found")
