import re
import os
import argparse
import mrcfile
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from scipy.interpolate import RegularGridInterpolator as rgi
from cardiac_desmosome.constants.host_system_constants import CHIMERAX_RUN_CMD


def calculate_metrics(
    vector1: np.ndarray,
    vector2: np.ndarray, 
    filter_zero: bool=False,
) -> tuple:
    """ Calculate metrics between two vectors

    Author: Satwik Pasani

    ## Arguments:

    - **vector1 (_type_)**:<br />
        Vector 1.

    - **vector2 (_type_)**:<br />
        Vector 2.

    - **filter_zero (bool, optional):**:<br />
        If True, zero values in vector2 will be filtered out from both vectors.

    ## Returns:

    - **tuple**:<br />
        A tuple containing:
        - **overlap (float)**:<br />
            The sum of element-wise products of vector1 and vector2.
        - **corr_over_zero (float)**:<br />
            Correlation coefficient over non-zero elements.
        - **corr_over_mean (float)**:<br />
            Correlation coefficient after mean subtraction.
        - **(corr_pearson, corr_spearman) (tuple)**:<br
            Pearson and Spearman correlation coefficients.
    """

    vector1, vector2 = vector1.flatten(), vector2.flatten()

    if filter_zero:
        vector1 = vector1[vector2 != 0]
        vector2 = vector2[vector2 != 0]

    overlap = np.sum(vector1 * vector2)
    overlap_mean_sub = np.sum(
        (vector1 - vector1.mean()) * (vector2 - vector2.mean())
    )

    mag1 = np.sqrt(np.sum(vector1 ** 2))
    mag2 = np.sqrt(np.sum(vector2 ** 2))

    mag_meaned1 = np.sqrt(np.sum((vector1 - vector1.mean()) ** 2))
    mag_meaned2 = np.sqrt(np.sum((vector2 - vector2.mean()) ** 2))

    corr_over_zero = overlap / mag1 / mag2
    corr_over_mean = overlap_mean_sub / mag_meaned1 / mag_meaned2

    corr_pearson = pearsonr(vector1, vector2).statistic
    corr_spearman, pval = spearmanr(vector1, vector2)

    return overlap, corr_over_zero, corr_over_mean, (corr_pearson, corr_spearman)

def calculate_with_external_grid_with_addition_two_lists(
    grid_points: list,
    grid1_values: list,
    grid2_values: list,
    name: str,
    voxel_spacing: int,
):
    """ Calculate metrics between two lists of grid points and their 
    corresponding values.

    Author: Satwik Pasani
    Mofified by: OMG

    ## Arguments:

    - **grid_points (list)**:<br />
        List of grid points, where each point is a tuple of (x, y, z).

    - **grid1_values (list)**:<br />
        List of values corresponding to the first set of grid points.

    - **grid2_values (list)**:<br />
        List of values corresponding to the second set of grid points.

    - **name (str)**:<br />
        Name for the output, used in print statements.

    - **voxel_spacing (int, optional):**:<br />
        Spacing for the voxel grid. Defaults to 5.
    """

    interpolators = []

    for g, v in zip(grid_points, grid1_values + grid2_values):
        interpolators.append(rgi(g, v, bounds_error=False, fill_value=0))

    xs, ys, zs = [], [], []
    for g in grid_points:
        xs.append(g[0])
        ys.append(g[1])
        zs.append(g[2])

    xfinal = np.arange(
        np.min(np.hstack(xs)), np.max(np.hstack(xs)), voxel_spacing
    )
    yfinal = np.arange(
        np.min(np.hstack(ys)), np.max(np.hstack(ys)), voxel_spacing
    )
    zfinal = np.arange(
        np.min(np.hstack(zs)), np.max(np.hstack(zs)), voxel_spacing
    )

    xfinal, yfinal, zfinal = np.meshgrid(xfinal, yfinal, zfinal, indexing='ij')
    desired_points = np.array([
        (i, j, k) for i, j, k in zip(
            xfinal.flatten(), yfinal.flatten(), zfinal.flatten()
        )
    ])

    values = []
    for i in interpolators:
        values.append(i(desired_points).reshape(xfinal.shape))

    sum1 = np.sum(values[:len(grid1_values)], axis=0)
    sum2 = np.sum(values[len(grid1_values):], axis=0)

    overlap1, corr1, cam1, (pearson_c1, spearman_c1) = calculate_metrics(
        vector1=sum1, 
        vector2=sum2, 
        filter_zero=False
    )

    overlap2, corr2, cam2, (pearson_c2, spearman_c2) = calculate_metrics(
        vector1=sum2, 
        vector2=sum1, 
        filter_zero=False
    )

    correlation_list = [{
        "model_mrc": name,
        "ref_mrc": name,
        "correlation": corr1,
        "cam": cam1,
    }, {
        "model_mrc": name,
        "ref_mrc": name,
        "correlation": corr2,
        "cam": cam2,
    }]

    return correlation_list

def parse_chimerax_correlation_log(
    chimerax_log: str, 
    command: str="fitmap",
) -> list:
    """ Parse the log file from ChimeraX output.

    ## Arguments:

    - **chimerax_log (str)**:<br />
        Path to the ChimeraX log file.

    - **command (str, optional):**:<br />
        The command for which to parse the log. 
        Options are "fitmap" and "measure correlation". Default is "fitmap".

    ## Returns:

    - **list**:<br />
        A list of dictionaries containing the parsed correlation results.
    """

    attrs = {
        "fitmap": {
            "look_for": "Fit map",
            "regex1": r"Fit map (.+) in map (.+) using (\d+) points",
            "regex2": r"  correlation\s*=\s*([0-9]+\.[0-9]+),\s*correlation\s+about\s+mean\s*=\s*([0-9]+\.[0-9]+),\s*overlap\s*=\s*([0-9]+\.[0-9]+)"
        },
        "measure correlation": {
            "look_for": "Correlation of",
            "regex1": r"Correlation\s+of\s+(\S+)\s+#(\d+)\s+above\s+level\s+([0-9.]+e[+-]?[0-9]+)\s+in\s+(\S+)\s+#(\d+)",
            "regex2": r"correlation\s*=\s*([0-9]+\.[0-9]+),\s*correlation\s+about\s+mean\s*=\s*([0-9]+\.[0-9]+)"
        }
    }

    assert command in attrs, f"Command {command} not recognized. Valid options are: {list(attrs.keys())}"

    log_lines = []
    with open(chimerax_log, 'r') as f:
        log_lines = f.readlines()

    correlation_list = []

    for idx, line in enumerate(log_lines):
        if (
            idx + 1 >= len(log_lines) or
            line.startswith(attrs[command]["look_for"]) is False
        ):
            continue

        regex1 = attrs[command]["regex1"]
        regex2 = attrs[command]["regex2"]

        match1 = re.match(regex1, line)
        match2 = re.match(regex2, log_lines[idx + 1].strip())

        if command == "fitmap":
            model_map, ref_map, num_points = match1.groups()
            correlation, cam, overlap = match2.groups()

        elif command == "measure correlation":
            model_map, ref_map = match1.groups()
            correlation, cam = match2.groups()    

        correlation_list.append({
            "model_mrc": model_map.replace("mrc", "").strip(),
            "ref_mrc": ref_map.replace("mrc", "").strip(),
            "correlation": correlation,
            "cam": cam,
        })

    return correlation_list

if __name__ == "__main__":

    parser = argparse.ArgumentParser("""
    Correlate LPD maps from Sample_A and Sample_B directories.
    1. If --use_combined_map is set, combines all LPD maps in Sample_A and Sample_B
       respectively and computes correlation between the combined maps.
    2. If --use_combined_map is not set, computes correlation for each pair of LPD maps
       from Sample_A and Sample_B.
    3. Supports two modes for correlation calculation: 'pasani' and 'chimerax'.
       - 'pasani': Uses custom interpolation and correlation calculation.
       - 'chimerax': Uses ChimeraX commands to compute correlations.
    """)

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

    parser.add_argument(
        "--use_combined_map",
        action='store_true',
        help="Whether to use the combined map for correlation calculation. If not set, individual LPD maps will be compared."
    )

    args = parser.parse_args()

    pa = os.path.join(args.sampcon_cluster_path, 'Sample_A')
    pb = os.path.join(args.sampcon_cluster_path, 'Sample_B')

    fa = [os.path.join(pa,f"{f}") for f in os.listdir(pa) if 'LPD' in f]
    fb = [os.path.join(pb,f"{f}") for f in os.listdir(pb) if 'LPD' in f]

    if args.mode == "pasani":

        if args.use_combined_map:

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

            correlation_list = calculate_with_external_grid_with_addition_two_lists(
                grid_points=grid_points,
                grid1_values=value_points[:len(mrc_a)],
                grid2_values=value_points[len(mrc_a):],
                name='combined_map',
                voxel_spacing=6,
            )

        else:
            correlation_list = []
            for f1, f2 in zip(fa, fb):
                mrc_a = mrcfile.open(f1, 'r', permissive=True)
                mrc_b = mrcfile.open(f2, 'r', permissive=True)

                grid_points = []
                value_points = []

                for mrc in [mrc_a, mrc_b]:
                    xvals = mrc.voxel_size.x * np.arange(mrc.data.shape[2]) + mrc.header.origin.x
                    yvals = mrc.voxel_size.y * np.arange(mrc.data.shape[1]) + mrc.header.origin.y
                    zvals = mrc.voxel_size.z * np.arange(mrc.data.shape[0]) + mrc.header.origin.z
                    grid_points.append((xvals, yvals, zvals))
                    value_points.append(mrc.data.transpose(2, 1, 0).copy())

                correlation_list += calculate_with_external_grid_with_addition_two_lists(
                    grid_points=grid_points,
                    grid1_values=value_points[:1],
                    grid2_values=value_points[1:],
                    name=os.path.basename(f1).replace('.mrc', ''),
                    voxel_spacing=6,
                )

    elif args.mode == "chimerax":

        log_path = os.path.join(args.sampcon_cluster_path, 'correlation_logs.txt')

        if args.use_combined_map:

            commands_run = []

            for f1 in fa:
                commands_run += [f"open {f1}"]

            for f2 in fb:
                commands_run += [f"open {f2}"]

            commands_run += [
                f"volume add #1-{len(fa)}",
                f"volume add #{len(fa)+1}-{len(fa)+len(fb)}",
                f"fitmap #{len(fa)+len(fb)+1} in_map #{len(fa)+len(fb)+2} shift false rotate false envelope false zeros true maxSteps 0",
                f"fitmap #{len(fa)+len(fb)+2} in_map #{len(fa)+len(fb)+1} shift false rotate false envelope false zeros true maxSteps 0",
                f"close all",
            ]

            os.system(
                f"{CHIMERAX_RUN_CMD} --exit --nogui --cmd " +
                f"'{"; ".join(commands_run)}' > {log_path}"
            )

        else:

            for idx, (f1, f2) in enumerate(zip(fa, fb)):

                commands_run = [
                    f"open {f1}",
                    f"open {f2}",
                    f"fitmap #1 in_map #2 shift false rotate false envelope false zeros true",
                    f"fitmap #2 in_map #1 shift false rotate false envelope false zeros true",
                    f"close all",
                ]

                log_mode = ">>" if idx > 0 else ">"

                os.system(
                    f"{CHIMERAX_RUN_CMD} --exit --nogui --cmd " +
                    f"'{"; ".join(commands_run)}' {log_mode} {log_path}"
                )

        correlation_list = parse_chimerax_correlation_log(log_path)

    df = pd.DataFrame(correlation_list)
    # print(df)

    corr_ = df["correlation"].astype(float)
    cam_ = df["cam"].astype(float)
    print(f"Average Correlation: {corr_.mean():.4f} ± {corr_.std():.4f}")
    print(f"Average CAM: {cam_.mean():.4f} ± {cam_.std():.4f}")

    out_csv = os.path.join(args.sampcon_cluster_path, 'correlation_results.csv')
    df.to_csv(out_csv, index=False)
    print(f"Correlation results saved to {out_csv}")