import os
import argparse
import mrcfile
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from scipy.interpolate import RegularGridInterpolator as rgi
from pprint import pprint
from cardiac_desmosome.constants.host_system_constants import CHIMERAX_RUN_CMD


def calculate_metrics(v1, v2, filter_zero=False):
    """ Calculate various metrics between two vectors.

    Args:
        v1 (np.ndarray): Vector 1.
        v2 (np.ndarray): Vector 2.
        filter_zero (bool, optional): If True, zero values in v2 will be filtered out from both vectors. Defaults to False.

    Returns:
        tuple: A tuple containing:
            - overlap (float): The sum of element-wise products of v1 and v2.
            - corr_over_zero (float): Correlation coefficient over non-zero elements.
            - corr_over_mean (float): Correlation coefficient after mean subtraction.
            - (corr_pearson, corr_spearman) (tuple): Pearson and Spearman correlation coefficients.
    """

    v1, v2 = v1.flatten(), v2.flatten()

    if filter_zero:
        v1 = v1[v2 != 0]
        v2 = v2[v2 != 0]

    overlap = np.sum(v1 * v2)
    overlap_mean_sub = np.sum((v1 - v1.mean()) * (v2 - v2.mean()))

    mag1 = np.sqrt(np.sum(v1 ** 2))
    mag2 = np.sqrt(np.sum(v2 ** 2))

    mag_meaned1 = np.sqrt(np.sum((v1 - v1.mean()) ** 2))
    mag_meaned2 = np.sqrt(np.sum((v2 - v2.mean()) ** 2))

    corr_over_zero = overlap / mag1 / mag2
    corr_over_mean = overlap_mean_sub / mag_meaned1 / mag_meaned2

    corr_pearson = pearsonr(v1, v2).statistic
    corr_spearman, pval = spearmanr(v1, v2)

    return overlap, corr_over_zero, corr_over_mean, (corr_pearson, corr_spearman)


def calculate_with_external_grid_with_addition_two_lists(
    list_of_g,
    list_of_v1,
    list_of_v2,
    name,
    voxel_spacing=5
):
    """ Calculate metrics between two lists of grid points and their corresponding values.

    Args:
        list_of_g (list of tuples): List of grid points, where each point is a tuple of (x, y, z).
        list_of_v1 (list of np.ndarray): List of values corresponding to the first set of grid points.
        list_of_v2 (list of np.ndarray): List of values corresponding to the second set of grid points.
        name (str): Name for the output, used in print statements.
        voxel_spacing (int, optional): Spacing for the voxel grid. Defaults to 5.
    """

    # the two "v" lists are added and compared
    interpolators = []
    for g, v in zip(list_of_g, list_of_v1 + list_of_v2):
        interpolators.append(rgi(g, v, bounds_error=False, fill_value=0))
    xs, ys, zs = [], [], []
    for g in list_of_g:
        xs.append(g[0])
        ys.append(g[1])
        zs.append(g[2])

    xfinal = np.arange(np.min(np.hstack(xs)), np.max(np.hstack(xs)), voxel_spacing)
    yfinal = np.arange(np.min(np.hstack(ys)), np.max(np.hstack(ys)), voxel_spacing)
    zfinal = np.arange(np.min(np.hstack(zs)), np.max(np.hstack(zs)), voxel_spacing)

    xfinal, yfinal, zfinal = np.meshgrid(xfinal, yfinal, zfinal, indexing='ij')
    desired_points = np.array([
        (i, j, k) for i, j, k in zip(xfinal.flatten(), yfinal.flatten(), zfinal.flatten())
    ])

    values = []
    for i in interpolators:
        values.append(i(desired_points).reshape(xfinal.shape))

    sum1 = np.sum(values[:len(list_of_v1)], axis=0)
    sum2 = np.sum(values[len(list_of_v1):], axis=0)

    o, c0, cm, cp = calculate_metrics(sum1, sum2)
    print(f'{name} (all points)')
    print(f'\tOverlap: {o:.4f}, C0: {c0:.4f}, Cam: {cm:.4f} ({cp[0]:.4f}, {cp[1]:.4f})')

    o, c0, cm, cp = calculate_metrics(sum1, sum2, False)
    print(f'{name} (all points, filter=False)')
    print(f'\tOverlap: {o:.4f}, C0: {c0:.4f}, Cam: {cm:.4f} ({cp[0]:.4f}, {cp[1]:.4f})')


threshold = {
    "DP-N": 0.02925, # 25% of 0.117
    "DP-S": 0.1084, # 20% of 0.542
    "DSC-C1": 0.02943, # 27% of 0.109
    "DSC-C2": 0.0244, # 25% of 0.0813
    "DSC-C3": 0.026, # 25% of 0.104
    "DSG-C1": 0.023328, # 24% of 0.0972
    "DSG-C2": 0.017444, # 28% of 0.0623
    "DSG-C3": 0.0202, # 20% of 0.101
    "DSG-C4": 0.02544, # 16% of 0.159
    "PG-C": 0.017992, # 26% of 0.0692
    "PG-N": 0.015525, # 25% of 0.0621
    "PG-S": 0.17955, # 35% of 0.513
    "PKP-C": 0.1908, # 30% of 0.636
    "PKP-N1": 0.07128, # 25% of 0.264
    "PKP-N2": 0.0465, # 30% of 0.155
    "PKP-N3": 0.028192, # 32% of 0.0881
    "PKP-N4": 0.024875, # 25% of 0.0995
}

def parse_chimerax_log(chimerax_log: str):

    log_lines = []
    with open(chimerax_log, 'r') as f:
        log_lines = f.readlines()

    correlation_list = []

    for idx, line in enumerate(log_lines):
        if (
            idx + 1 >= len(log_lines) or
            line.startswith("Correlation of") is False
        ):
            continue
        line_parts = line.split(" ")
        nxt_line_parts = log_lines[idx + 1].split(" ")
        model_mrc = line_parts[2]
        ref_mrc = line_parts[-2]
        correlation = nxt_line_parts[2].split(",")[0].strip()
        cam = nxt_line_parts[-1].strip()
        correlation_list.append({
            "model_mrc": model_mrc,
            "ref_mrc": ref_mrc,
            "correlation": correlation,
            "cam": cam,
        })

    return correlation_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--sampcon_cluster_path',
        type=str,
        required=True,
        help='Path to the directory containing Sample_A and Sample_B'
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["pasani", "chimerax"],
        default="pasani",
        help="Method to use for correlation calculation. Default is 'chimerax'."
    )
    args = parser.parse_args()


    pa = os.path.join(args.sampcon_cluster_path, 'Sample_A')
    pb = os.path.join(args.sampcon_cluster_path, 'Sample_B')

    fa = [os.path.join(pa,f"{f}") for f in os.listdir(pa) if 'LPD' in f]
    fb = [os.path.join(pb,f"{f}") for f in os.listdir(pb) if 'LPD' in f]

    # pprint(fa)
    # pprint(fb)

    if args.mode == "pasani":

        mrc_a = [mrcfile.open(x, 'r', permissive=True) for x in fa]
        mrc_b = [mrcfile.open(x, 'r', permissive=True) for x in fb]

        grid_points = []
        value_points = []

        for mrc in mrc_a + mrc_b:
            xvals = mrc.voxel_size.x * np.arange(mrc.data.shape[2]) + mrc.header.origin.x
            yvals = mrc.voxel_size.y * np.arange(mrc.data.shape[1]) + mrc.header.origin.y
            zvals = mrc.voxel_size.z * np.arange(mrc.data.shape[0]) + mrc.header.origin.z
            grid_points.append((xvals, yvals, zvals))
            value_points.append(mrc.data.transpose(2, 1, 0).copy())

        calculate_with_external_grid_with_addition_two_lists(
            grid_points,
            value_points[:len(mrc_a)],
            value_points[len(mrc_a):],
            'Sample A <-> Sample B'
        )

    elif args.mode == "chimerax":

        log_path = os.path.join(args.sampcon_cluster_path, 'correlation_logs.txt')

        for idx, (f1, f2) in enumerate(zip(fa, fb)):
            id1 = os.path.basename(f1).split('LPD_')[1].split('.mrc')[0]
            id2 = os.path.basename(f2).split('LPD_')[1].split('.mrc')[0]
            thr1 = threshold[id1]
            thr2 = threshold[id2]
            commands_run = [
                f"open {f1}",
                f"open {f2}",
                f"volume #1 sdLevel {thr1}",
                f"volume #2 sdLevel {thr2}",
                f"measure correlation #1 in_map #2 envelope false",
                f"measure correlation #2 in_map #1 envelope false",
                f"close all",
            ]
            log_mode = ">>" if idx > 0 else ">"
            os.system(
                f"{CHIMERAX_RUN_CMD} --exit --nogui --cmd " +
                f"'{"; ".join(commands_run)}' {log_mode} {log_path}"
            )
        correlation_list = parse_chimerax_log(log_path)
        # pprint(correlation_list)
        df = pd.DataFrame(correlation_list)
        df.to_csv(os.path.join(args.sampcon_cluster_path, 'correlation_results.csv'), index=False)