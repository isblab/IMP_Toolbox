# Originally written by Shreyas Arvindekar. Modified by stochastic13 (Satwik)
'''
Algorithm

Examine a set of multipliers (mean, mean-0.25 std, .. )
for a set of data restraints

Choose the largest number of models (lowest multiplier)
[such that nA and nB are each less than 15k: earlier ] about the same as :
    nA +nB < 20k

It passes the KS test on total score and (A,B) looks similar on score
distribution

If you dont find a multiplier even after the narrower search in point 1,
[take a random subset of the nearest multiplier that passes the KS test]
    OR
take a single lenient cutoff on EV (less than mean) along with score
multipliers on the other restraints.
'''

import os
import time
import argparse
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from matplotlib import pyplot as plt

start_time = time.time()

# This function will take the individual dataframes and compare them with the
# mean-multiplier*std of common dataframe and output the dataframe of the
# models that satisfy the filter
def variable_filter(
    std_multiplier: float,
    df: pd.DataFrame,
    score_list: list,
    df_mean: pd.Series,
    df_std: pd.Series,
):
    """ Apply a variable filter based on standard deviation multipliers.

    Args:
        std_multiplier (float): Multiplier for standard deviation.
        df (pd.DataFrame): DataFrame to filter.
        score_list (list): List of score columns to consider.
        df_mean (pd.Series): Mean values for each score column in combined df.
        df_std (pd.Series): Standard deviation values for each score column in combined df.

    Returns:
        pd.DataFrame: Filtered DataFrame.
    """

    temp = None

    for i in score_list:

        if temp is None:
            temp = df[i] <= (df_mean[i] + std_multiplier * df_std[i])
        else:
            temp = temp & (df[i] <= (df_mean[i] + std_multiplier * df_std[i]))

    return df[temp]

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c',
        '--cluster_num',
        required=True,
        type=int,
        help='On which cluster shall I run the filter?'
    )
    parser.add_argument(
        '-lc',
        '--lowest_cutoff',
        default=-2,
        help='What is the standard deviation multiplier for the most stringent cutoff?'
    )
    parser.add_argument(
        '-hc',
        '--highest_cutoff',
        default=3,
        help='What is the standard deviation multiplier for the most lenient cutoff?'
    )
    parser.add_argument(
        '-ss',
        '--step_size',
        default=0.01,
        help='What step size shall I use?'
    )
    parser.add_argument(
        '-n',
        '--num_models',
        default=30000,
        type=int,
        help='What is the maximum number of models you want?'
    )
    parser.add_argument(
        '-g',
        '--gsmsel',
        default=os.path.join("/data/omkar/imp_toolbox_test/analysis"),
        type=str,
        help='Where is the gsm_sel directory'
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        default=os.path.join("/data/omkar/imp_toolbox_test/analysis", 'variable_filter_output'),
        type=str,
        help="Path to save theoutput files (default: analysis_output_path/set1/variable_filter_output)",
    )
    parser.add_argument(
        "--restraint_handles",
        nargs="+",
        default=[
            'EV_sum',
            'EM_sum',
        ],
        help="List of restraint handles to analyze (default: ['EV_sum', 'EM_sum'])"
    )

    args = parser.parse_args()

    cluster_num = str(args.cluster_num)
    lowest_cutoff = float(args.lowest_cutoff)
    highest_cutoff = float(args.highest_cutoff)
    step_size = float(args.step_size)
    num_models = int(args.num_models)
    gsm_sel_dir = str(args.gsmsel)

    print(args)

    data_restraint_names = []

    cluster_csv_fileA = os.path.join(
        gsm_sel_dir, f'selected_models_A_cluster{cluster_num}_detailed.csv'
    )
    cluster_csv_fileB = os.path.join(
        gsm_sel_dir, f'selected_models_B_cluster{cluster_num}_detailed.csv'
    )

    std_mult_dtrst = []
    for i in np.arange(highest_cutoff, lowest_cutoff - step_size, -step_size):
        std_mult_dtrst.append(round(i, 3))

    # Reading the CSV files and combining them to apply a common cutoff
    dfA = pd.read_csv(cluster_csv_fileA)
    dfB = pd.read_csv(cluster_csv_fileB)
    print('Loaded the csv files')

    df_list = [dfA, dfB]
    common_df = pd.concat(df_list, ignore_index=True)
    print(common_df.head())

    print('Calculating mean')
    common_df_mean = common_df.mean(numeric_only=True)

    print('Calculating standard deviation')
    common_df_std = common_df.std(numeric_only=True)

    print('Running variable filter')
    out = []
    out_str = (
        f'Cluster number: {cluster_num} \
        \nLowest cutoff: {lowest_cutoff} \
        \nHighest cutoff: {highest_cutoff} \
        \nStep size: {step_size} \
        \nMaximum number of models to be selected: {num_models} \n\n'
    )
    mult_found = False

    for multiplier in std_mult_dtrst:

        default_restraints = args.restraint_handles + ['Total_Score']

        sel_dfA = variable_filter(
            multiplier, dfA, default_restraints, common_df_mean, common_df_std
        )
        sel_dfB = variable_filter(
            multiplier, dfB, default_restraints, common_df_mean, common_df_std
        )

        # Combining the score files for checking for run representation.
        sel_df_list = [sel_dfA, sel_dfB]
        sel_common_df = pd.concat(sel_df_list, ignore_index=True)
        nModelsT = len(sel_common_df.index)
        nModelsA = len(sel_dfA.index)
        nModelsB = len(sel_dfB.index)
        nRunsA = sel_dfA.traj.nunique()
        nRunsB = sel_dfB.traj.nunique()

        # Obtaining the scoresA and scoresB for sampling convergence
        scoresA = list(sel_dfA['Total_Score'])
        scoresB = list(sel_dfB['Total_Score'])
        scores = list(sel_common_df['Total_Score'])

        # Check if the two score distributions are similar
        ksd_pval = ks_2samp(scoresA, scoresB)
        ksd = ksd_pval[0]
        ksp = ksd_pval[1]

        # out = np.append(out,[multiplier, nModelsA, nModelsB, ksd, ksp])
        results = [multiplier, nModelsA, nModelsB, nRunsA, nRunsB, ksd, ksp]
        out.append(results)
        if nModelsA + nModelsB <= args.num_models:
            if (ksp > 0.05) or ((ksp <= 0.05) and (ksd < 0.3)):
                mult_found = True
                break

    out_df = pd.DataFrame(
        out,
        columns=[
            'DRest_Multiplier',
            'nModelsA',
            'nModelsB',
            'nRunsA',
            'nRunsB',
            'KS_D-value',
            'KS_p-value'
        ]
    )
    print(out_df.head())

    if mult_found == True:
        print(f'\nOptimal filter found.\nExtracted at {multiplier}')
        print(out_str)
        print(out_df)
        nBins = int(max(scores) - min(scores))
        print(nBins)
        plt.figure()
        plt.hist(scoresA, bins=nBins, histtype='step', label='ScoresA')
        plt.hist(scoresB, bins=nBins, histtype='step', label='ScoresB')
        plt.title('Scores of sampleA and sampleB')
        plt.xlabel('Total Score')
        plt.ylabel('nModels')
        plt.legend()
        os.makedirs(args.output_dir, exist_ok=True)
        plt.savefig(os.path.join(args.output_dir, 'var_filt_out.png'))
        # plt.show()
        fileA_basename = os.path.basename(cluster_csv_fileA)
        fileB_basename = os.path.basename(cluster_csv_fileB)
        sel_dfA.to_csv(
            # cluster_csv_fileA.replace('selected_models', 'good_scoring_models')
            os.path.join(
                args.output_dir,
                fileA_basename.replace('selected_models', 'good_scoring_models')
            )
        )
        sel_dfB.to_csv(
            # cluster_csv_fileB.replace('selected_models', 'good_scoring_models')
            os.path.join(
                args.output_dir,
                fileB_basename.replace('selected_models', 'good_scoring_models')
            )
        )
    else:
        print('\nOptimal multiplier not found')

    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")