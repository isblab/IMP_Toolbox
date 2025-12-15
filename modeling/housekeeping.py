import os
import re
import time
import traceback
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from tqdm import tqdm
from argparse import ArgumentParser
from equilibration import detectEquilibration

def get_stat_file_names(path: str, to_sort: bool=True) -> tuple:
    """ Get the names of the stat and stat_replica files in the given path.

    Args:
        path (str):
            Path to the directory containing stat and stat_replica files.
        to_sort (bool, optional):
            Whether to sort the files based on the number in their names.
            Defaults to True.

    Returns:
        tuple:
            A tuple containing two lists:
            - stat_files (list): List of stat file names.
            - stat_replica_files (list): List of stat_replica file names.
    """

    files = os.listdir(path)

    stat_files = [
        x for x in files if re.search(r"stat[.][0-9]+[.]out", x)
    ]
    stat_replica_files = [
        x for x in files if re.search(r"stat_replica[.][0-9]+[.]out", x)
    ]

    if to_sort:
        stat_files = sorted(
            stat_files, key=lambda x: int(x.split(".")[1])
        )
        stat_replica_files = sorted(
            stat_replica_files, key=lambda x: int(x.split(".")[1])
        )

    return stat_files, stat_replica_files

def get_stat_head(path: str, stat_files: list) -> tuple:
    """ Get the header information from the first line of the first stat file.
    The header contains the index -> heading mapping and the environ details
    ```
    head_dict = {
        0: 'MonteCarlo_Nframe',
        1: 'Restraint_Score',
    }
    inv_head_dict = {
        'MonteCarlo_Nframe': 0,
        'Restraint_Score': 1,
    }
    ```

    Args:
        path (str):
            Path to the directory containing stat and stat_replica files.
        stat_files (list):
            List of stat or stat_replica file names.

    Returns:
        tuple:
            A tuple containing two dictionaries:
            - head_dict (dict): Dictionary mapping index to heading.
            - inv_head_dict (dict): Dictionary mapping heading to index.
    """

    with open(os.path.join(path, stat_files[0])) as f:

        first_line = f.readline().strip()

        head_info = r"^.+\'STAT2HEADER_IMP_VERSIONS\'[:][ ]*([\"\'])[{].*[}]\1"
        match = re.match(
            head_info,
            first_line,
            flags=re.DOTALL
        )

        assert match is not None, "Header structure different"

        head_dict = eval(first_line)

    inv_head_dict = dict()

    for idx in head_dict:
        if isinstance(idx, int):
            inv_head_dict[head_dict[idx]] = idx

    return head_dict, inv_head_dict

def check_order(flat_frame_order: list) -> None:
    """ Check the order of the collated frame order.
    Check that not more than 1 frame reports a 1.0 temperature at a time

    Args:
        flat_frame_order (list):
            Frame indices where the temperature is 1.0 across all replicas.
    """

    assert len(flat_frame_order) == len(
        set(flat_frame_order)
    ), "Some frames reported as 1.0 temperature in multiple replicas"

    assert max(flat_frame_order) == (
        len(flat_frame_order) - 1
    ), "Some frames missing in the frame order"

def missing_stat_file_check(stat_files: list, stat_replica_files: list) -> None:
    """ Check for missing stat or stat_replica files.
    1. Check if there is no missing number in the stat files
    2. Check if a corresponding replica file is present for each stat file
    3. Check if the number of stat and stat_replica files are the same

    Args:
        stat_files (list):
            List of stat files
        stat_replica_files (list):
            List of stat_replica files
    """

    set_a = set([int(i.split(".")[1]) for i in stat_files])
    set_b = set(list(range(len(stat_files))))
    assert set_a == set_b, "Stat files are not as expected"

    for i in stat_files:
        x = ["stat_replica"] + i.split(".")[1:]
        assert (".".join(x) in stat_replica_files), (
            f"Corresponding stat_replica file missing for stat.{i}.out"
        )

    assert len(stat_files) == len(stat_replica_files), (
        "Number of stat and stat_replica files are different"
    )

def get_stat_stack(path: str, stat_files: list) -> list:
    """ Get the stack from the given stat files.
    replica_stack_list = replicas x frames x header_fields
    stack_list = non_empty_replicas x temp_1_frames x header_fields

    Args:
        path (str):
            Path to the directory containing stat and stat_replica files.
        stat_files (list):
            List of stat or stat_replica file names.

    Returns:
        list:
            A list of numpy arrays, each array corresponds to a stat file.
    """

    stack_list = []

    for stat_file in stat_files:
        with open(os.path.join(path, stat_file)) as f:
            _rd = f.readline() # discard line 1
            frame_info_list = []

            for line in f:
                try:
                    line_dict = eval(line)

                except NameError:
                    print(line.strip())
                    line_dict = {}
                    for key_val_pair in line.split(","):
                        key, val = key_val_pair.split(":")
                        key = key.replace("{", "").replace("}", "").strip()
                        val = val.replace("{", "").replace("}", "").strip()
                        print(key, val)

                        if val == "inf":
                            line_dict[int(key)] = np.inf
                        else:
                            line_dict[int(key)] = eval(val)

                frame_info_list.append(line_dict)

            if len(frame_info_list) == 0:
                print(f"No lines parsed in {stat_file}, skipping this file.")
                continue

            stack_list.append(np.vstack([
                list(line_dict.values()) for line_dict in frame_info_list
            ]))

    return stack_list

def get_frame_order(
    inv_replica_head_dict: dict,
    replica_stack_list: list
) -> list:
    """ Get the frame order from the replica stack list.
    frame_order = replicas x temp_1_frames

    Args:
        inv_replica_head_dict (dict):
            Dictionary mapping heading to index for replica files.
        replica_stack_list (list):
            List of numpy arrays, each array corresponds to a stat_replica file.

    Returns:
        list:
            A list of numpy arrays, each array contains indices of frames at
            temperature 1.0 for each replica.
    """

    frame_order = ([])

    # we only take the frames where the temperature is 1.0
    for replica_stack in replica_stack_list:
        # remember replica_stack is a matrix of frames x header fields
        # temperature at each frame for this replica
        current_temp_idx = inv_replica_head_dict["ReplicaExchange_CurrentTemp"]
        current_temp_arr = replica_stack[:, current_temp_idx].flatten()

        # boolean array marking the frames where temperature is 1.0
        temp_bool_arr = (current_temp_arr == "1.0")
        temp_1_idxs = np.where(temp_bool_arr)[0]
        # indices of the frames where temperature is 1.0
        frame_order.append(temp_1_idxs)

    return frame_order

def sanity_check_replica_stack_list(
    replica_stack_list: list,
    stat_replica_files: list,
) -> None:
    """ Sanity check the replica stack list.
    1. Check if the number of replica files is equal to the number of stacks
    2. Check if all arrays have the same number of frames
    3. Check if all arrays have the same number of header fields

    Args:
        replica_stack_list (list):
            List of numpy arrays, each array corresponds to a stat_replica file.
        stat_replica_files (list):
            List of stat_replica file names.
    """

    # number of stacks should be equal to number of replica files
    assert len(replica_stack_list) == len(stat_replica_files), (
        "Some replica files could not be parsed"
    )

    # all arrays should have the same number of frames
    assert len(
        set([replica_arr.shape[0] for replica_arr in replica_stack_list])
    ) == 1, "Replica stat files have different number of frames"

    # all arrays should have the same number of header fields
    assert len(
        set([replica_arr.shape[1] for replica_arr in replica_stack_list])
    ) == 1, "Replica stat files have different number of header fields"

def parser(run_dir: str) -> tuple:
    """ Parse the stat and stat_replica files in the given path.

    Args:
        path (str):
            Path to the directory containing stat and stat_replica files.

    Returns:
        tuple: A tuple containing:
            - per_frame_stats (list): List of dictionaries with stats for each frame.
            - per_frame_replica_stats (list): List of dictionaries with replica stats for each frame.
            - stack_list (list): List of numpy arrays, each array corresponds to a stat file.
            - replica_stack_list (list): List of numpy arrays, each array corresponds to a stat_replica file.
            - inv_head_dict (dict): Dictionary mapping header names to their indices for stat files.
            - inv_replica_head_dict (dict): Dictionary mapping header names to their indices for stat_replica files.
            - frame_order (list): List of numpy arrays, each array contains indices of frames at temperature 1.0 for each replica.
    """

    stat_files, stat_replica_files = get_stat_file_names(run_dir)

    missing_stat_file_check(stat_files, stat_replica_files)

    ###########################################################################
    # getting header information
    # assumes header is same across all stat files
    ###########################################################################
    _head_dict, inv_head_dict = get_stat_head(
        run_dir, stat_files
    )
    _replica_head_dict, inv_replica_head_dict = get_stat_head(
        run_dir, stat_replica_files
    )

    ###########################################################################
    # parsing the stat_replica files
    ###########################################################################

    replica_stack_list = get_stat_stack(run_dir, stat_replica_files)

    sanity_check_replica_stack_list(replica_stack_list, stat_replica_files)

    frame_order = get_frame_order(inv_replica_head_dict, replica_stack_list)
    # all replicas should be represented in the frame order (even those with no temp 1.0 frames)
    assert len(frame_order) == len(replica_stack_list), (
        "Frame order length does not match number of replicas"
    )

    # Flattened list version of frame_order
    flat_frame_order = [
        frame_idx
        for frame_idxs in frame_order
        for frame_idx in frame_idxs.tolist()
    ]

    check_order(flat_frame_order)

    ###########################################################################
    # parsing the stat files
    ###########################################################################

    stack_list = get_stat_stack(run_dir, stat_files)
    # contains the parsed dictionaries according to the stat_files
    # each entry in stack_list corresponds to a non-empty stat file
    # each entry is a matrix of temp_1_frames x header fields

    assert len(flat_frame_order) == sum(
        [replica_arr.shape[0] for replica_arr in stack_list]
    ), "All frames not included in frame_order/array"

    assert len(flat_frame_order) == len(
        replica_stack_list[0]
    ); "All frames not included in frame_order/array_replica"


    flat_stack_list = [
        frame_stat_info
        for idx in range(len(stack_list))
        for frame_stat_info in stack_list[idx].tolist()
    ]

    # per_frame_stats is sorted flat_stack_list by frame idx
    # flat_stack_list = frames x header fields
    per_frame_stats = sorted(zip(flat_frame_order, flat_stack_list))
    per_frame_stats = [frame_stats[1] for frame_stats in per_frame_stats]

    monte_carlo_nframes = np.array([
        int(frame_stats[inv_head_dict["MonteCarlo_Nframe"]])
        for frame_stats in per_frame_stats
    ])
    frame_diff = np.diff(monte_carlo_nframes)

    # confirm that MonteCarlo_Nframe matches the order based on temperature 1
    assert len(np.unique(frame_diff)) == 1 and np.unique(frame_diff)[0] == 1, (
        "Temperature based frame ordering does not match MonteCarlo_Nframe"
    )

    sorted_replicas, _x = sort_the_replica_exchanges_lowest_temp(frame_order)

    flat_replica_stack_list = []

    for idx in range(len(flat_stack_list)):

        temp_1_frame_idx = flat_frame_order[idx]
        replica_idx = sorted_replicas[idx]
        per_replica_info = replica_stack_list[replica_idx]
        # this contains info for all frames
        # we only want the frames where temperature is 1.0
        req_replica_info = per_replica_info[temp_1_frame_idx]
        flat_replica_stack_list.append(req_replica_info.tolist())

    # per_frame_replica_stats is sorted flat_replica_stack_list by frame idx
    per_frame_replica_stats = sorted(
        zip(flat_frame_order, flat_replica_stack_list)
    )

    per_frame_replica_stats = [
        frame_rep_stats[1] for frame_rep_stats in per_frame_replica_stats
    ]

    return (
        per_frame_stats,
        per_frame_replica_stats,
        stack_list,
        replica_stack_list,
        inv_head_dict,
        inv_replica_head_dict,
        frame_order,
    )

def is_float(s):
    """
    Checks if a string can be successfully converted to a float.

    Args:
        s: The string to check.

    Returns:
        True if the string can be converted to a float, False otherwise.
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def avg_stats(
    per_frame_stats: list,
    frame_order: list,
    inv_head_dict: dict,
    burn_in: float=0.8,
):
    """ Calculate average stat from per-frame stats after a burn-in period.

    Args:
        per_frame_stats (list):
            scores and other stats per frame (frames x header fields)
        frame_order (list):
            frame indices where the temperature is 1.0 across all replicas
        inv_head_dict (dict):
            Dictionary mapping stat header names to their indices
        burn_in (float, optional):
            Fraction of frames to discard as burn-in. Defaults to 0.8.

    Returns:
        pd.DataFrame:
            DataFrame containing average stat for each key after burn-in.
    """

    sel_inv_head_dict = {}
    per_frame_stats = np.array(per_frame_stats)

    # only consider stat headers which are indexed
    sel_inv_head_dict = {
        key: h_idx
        for key, h_idx in inv_head_dict.items()
        if is_float(per_frame_stats[0, h_idx])
    }

    nframes = len(per_frame_stats)
    n_remove = int(burn_in * nframes)

    _s, exchange_indices = sort_the_replica_exchanges_lowest_temp(frame_order)

    avg_stats = {}
    rb_acceptance_scores = []
    srb_acceptance_scores = []
    ball_mover_acceptance_scores = []

    for stat_head in list(sel_inv_head_dict.keys()):

        if "Acceptance" not in stat_head:
            score = parse_key(
                stat_head,
                per_frame_stats,
                sel_inv_head_dict,
                exchange_indices,
                n_remove,
                adjust=False,
                exact_match=True
            )
            avg_stats[stat_head] = np.mean(score)

        else:
            score = parse_key(
                stat_head,
                per_frame_stats,
                sel_inv_head_dict,
                exchange_indices,
                n_remove,
                adjust=True,
                exact_match=True
            )
            avg_stats[stat_head] = np.nanmean(score)

            if "RigidBody" in stat_head or "_structure_" in stat_head:
                rb_acceptance_scores.extend(score)

            if "Super" in stat_head or "_srb_" in stat_head:
                srb_acceptance_scores.extend(score)

            if "_BallMover-" in stat_head:
                ball_mover_acceptance_scores.extend(score)

    avg_stats = {
        "Average_BallMover_Acceptance_Rate": np.nanmean(ball_mover_acceptance_scores),
        "Average_Rigid_Body_Acceptance_Rate": np.nanmean(rb_acceptance_scores),
        "Average_Super_Rigid_Body_Acceptance_Rate": np.nanmean(srb_acceptance_scores),
        **avg_stats
    }

    avg_stats_df = pd.DataFrame(
        avg_stats.items(),
        columns=["Key", f"Avg_last_{100-int(burn_in*100)}_percent"]
    )

    return avg_stats_df

def sort_the_replica_exchanges_lowest_temp(frame_order: list) -> tuple:
    """ Sort the replica exchanges based on the lowest temperature (1.0).

    Args:
        frame_order (list):
            List of numpy arrays, each array contains indices of frames at
            temperature 1.0 for each replica.

    Returns:
        tuple:
            A tuple containing:
            - sorted_replicas (list):
                List of replica indices sorted by frame order.
            - sorted_replicas_bool (list):
                List of boolean values indicating replica exchanges at each
                frame.
    """

    flat_frame_order = []
    replica_idxs = []

    for replica_idx in range(len(frame_order)):
        flat_frame_order += frame_order[replica_idx].tolist()
        replica_idxs += [
            replica_idx for _ in range(len(frame_order[replica_idx]))
        ]

    # Sort the replica indices according the frame_order
    sorted_replicas = sorted(zip(flat_frame_order, replica_idxs))
    sorted_replicas = [zipped_idx[1] for zipped_idx in sorted_replicas]

    # Boolean array marking exchanges
    sorted_replicas_bool = [False] + (np.diff(sorted_replicas) != 0).tolist()

    assert len(sorted_replicas_bool) == len(
        sorted_replicas
    ), "Length of exchange array does not match length of sorted replicas"

    return sorted_replicas, sorted_replicas_bool

def parse_key(
    search_string: str | list,
    per_frame_stats: list,
    sel_inv_head_dict: dict,
    exchange_indices: list,
    n_remove: int,
    adjust: bool=False,
    exact_match: bool=False,
    return_type: str="mean",
) -> np.ndarray:
    """ Parse a specific key from the per-frame stats and compute the mean.

    Args:
        stat_head (str):
            The key to search for in the stat headers.
        per_frame_stats (list):
            scores and other stats per frame (frames x header fields)
        sel_inv_head_dict (dict):
            Dictionary mapping selected stat header names to their indices
        exchange_indices (list):
            List of boolean values indicating replica exchanges at each frame.
        n_remove (int):
            Number of initial frames to discard as burn-in.
        adjust (bool, optional):
            Whether to adjust the values based on replica exchanges.
            Defaults to False.

    Returns:
        np.ndarray:
            Array of mean values for the specified key after burn-in.
    """

    if isinstance(search_string, str):
        if exact_match:
            assert search_string in sel_inv_head_dict, (
                f"{search_string} not found in stat headers"
            )
            stat_head_idxs = [sel_inv_head_dict[search_string]]
        else:
            stat_head_idxs = [
                sel_inv_head_dict[stat_head] for stat_head in sel_inv_head_dict.keys()
                if search_string in stat_head
            ]

    elif isinstance(search_string, list):
        exact_match = True
        assert all(
            [s in sel_inv_head_dict for s in search_string]
        ), "Some search strings not found in stat headers"
        stat_head_idxs = [
            sel_inv_head_dict[search_str] for search_str in search_string
        ]

    else:
        raise ValueError("search_string must be a string or a list of strings")

    stat_vals_list = [
        [float(frame_stats[stat_head_idx]) for frame_stats in per_frame_stats]
        for stat_head_idx in stat_head_idxs
    ]

    if adjust:
        adjusted_stat_vals_list = correct_mc_cumulative(
            stat_vals_list, exchange_indices
        )
        filtered_stat_vals = np.array(
            [x[n_remove:] for x in adjusted_stat_vals_list]
        )
        # filtered_stat_vals = adjusted_stat_vals_list[n_remove:]

    else:
        filtered_stat_vals = np.array(
            [x[n_remove:] for x in stat_vals_list]
        )
        # filtered_stat_vals = np.array(stat_vals_list[n_remove:])

    if return_type == "sum":
        return np.nansum(filtered_stat_vals, axis=0).flatten()
    elif return_type == "mean":
        return np.nanmean(filtered_stat_vals, axis=0).flatten()
    else:
        raise ValueError("return_type must be 'sum' or 'mean'")

def correct_mc_cumulative(
    mc_array: list,
    exchange_indices: list,
) -> list:
    """Change the MC acceptance ratios from cumulative to instantaneous (per frame)

    Args:
        mc_array (list):
            List of numpy arrays containing cumulative MC acceptance ratios.
        exchange_indices (list):
            List of boolean values indicating replica exchanges at each frame.

    Returns:
        np.ndarray:
            Array of instantaneous MC acceptance ratios with NaNs at exchange frames.
    """
    adjusted_mc_array = []

    for sub_mc in mc_array:
        num_frames = np.arange(1, len(sub_mc) + 1)
        corrected_array = [sub_mc[0]]

        # sub_mc is cumulative acceptance ratio until that frame
        total_skip_first = sub_mc[1:] * num_frames[1:]
        total_skip_last = sub_mc[:-1] * num_frames[:-1]
        acceptances = total_skip_first - total_skip_last

        corrected_array = corrected_array + acceptances.tolist()
        corrected_array = np.array(corrected_array)

        # At all exchanges, the replica changes, and hence, the numbers are not
        # valid for the frames where an exchange happened
        # TODO: Does the exchange happen before or after the MC steps?
        # This changes which all entries to NaN out
        corrected_array[np.array(exchange_indices)] = np.nan

        # if the exchange happened too many times, take the original array
        if np.sum(np.isnan(corrected_array)) > 0.8 * len(corrected_array):
            adjusted_mc_array.append(np.array(sub_mc))
            # return np.array(mc_array)
        else:
            adjusted_mc_array.append(corrected_array)
            # return corrected_array

    return adjusted_mc_array

def prod_run_worker(run_dir: str, burn_in: float=0.8) -> pd.DataFrame:
    """ Worker function to parse a run directory and compute average stats.

    Args:
        run_dir (str):
            Path to the directory containing stat and stat_replica files.
        burn_in (float, optional):
            Fraction of frames to discard as burn-in. Defaults to 0.8.

    Returns:
        pd.DataFrame:
            DataFrame containing average stats for the run directory.
    """

    (
        per_frame_stats,
        per_frame_replica_stats,
        stack_list,
        replica_stack_list,
        inv_head_dict,
        inv_replica_head_dict,
        frame_order,
    ) = parser(run_dir)

    avg_stats_df = avg_stats(
        per_frame_stats,
        frame_order,
        inv_head_dict,
        burn_in=burn_in,
    )

    return avg_stats_df

def sort_the_replica_exchanges_all_temp(
    replica_stack_list: list,
    inverted_dict_replica: dict
):
    """ Sort the replica exchanges based on all temperatures.
    Args:
        replica_stack_list (list): _description_
        inverted_dict_replica (dict): _description_

    Returns:
        _type_: _description_
    """

    x = np.array([0 for i in range(len(replica_stack_list[0]))], dtype=np.int32)
    # Total exchanges at each step

    for replica_stack in replica_stack_list:

        # Get the temperatures for each frame for this particular replica
        curr_replica_temp = replica_stack[
            :, inverted_dict_replica["ReplicaExchange_CurrentTemp"]
        ].flatten()
        curr_replica_temp = [float(x) for x in curr_replica_temp]

        # Check whenever the temperature changes indicating an exchange
        curr_replica_exchanges = [0] + (
            np.abs(np.diff(curr_replica_temp)) > 1e-5
        ).tolist()

        curr_replica_exchanges = np.array(curr_replica_exchanges, dtype=np.int32)

        x += curr_replica_exchanges

    # Since each exchange is added twice in the two exchanging replicas
    return x / 2

def get_moving_sd(score, frac=0.1):
    """ Calculate moving window standard deviation.

    Args:
        score (list): Score values.
        frac (float, optional): Fraction of total samples for window size. Defaults to 0.1.

    Returns:
        list: List of standard deviation values for each position.
    """

    # moving window sd calculation. each window is frac * all_samples
    sds = []
    for i in range(len(score)):
        inds = [
            max(0, int(i - (frac / 2) * len(score))),
            min(len(score) - 1, int(i + (frac / 2) * len(score))),
        ]
        sds.append(np.std(score[inds[0] : inds[1]]))
    return sds

def plot_trajectory_score(
    save_path,
    per_frame_stats,
    n_remove,
    exchange_indices,
    tot_score,
):
    traj_name = os.path.basename(save_path)
    n = len(per_frame_stats)
    n_keep = n - n_remove

    fitline = sm.nonparametric.lowess(
        tot_score[n_remove:],
        np.arange(n_remove, n),
        frac=0.2,
        return_sorted=False
    )

    plt.figure(figsize=(10, 10))
    plt.scatter(
        np.arange(n_remove, n),
        tot_score[n_remove:],
        color="black",
        zorder=0,
        s=10,
    )
    plt.scatter(
        np.arange(n_remove, n)[exchange_indices[n_remove:]],
        np.array(tot_score[n_remove:])[exchange_indices[n_remove:]],
        color="red",
        s=5,
        zorder=10,
    )
    plt.plot(
        np.arange(n_remove, n),
        fitline,
        color="green",
        lw=3,
        zorder=100,
    )
    plt.ylim(
        *np.sort(tot_score[n_remove:])[
            np.array([int(0.01 * n_keep), int(0.99 * n_keep)])
        ]
    )
    plt.xlabel("Frame Index")
    plt.ylabel("Total Score")
    plt.title(f"Trajectory {traj_name} Total Score")
    plt.legend(["Total Score", "Replica Exchange Points", "LOWESS Fit"], loc="upper right")
    plt.savefig(f"{save_path}/equilibriation_total_score.png")
    plt.close("all")

def geyer_equilibriation_test(score):
    """ Perform Geyer equilibration test on the given score trajectory.

    Args:
        score (list): Score values.

    Returns:
        int: Equilibration time index.
    """
    try:
        t_eq = detectEquilibration(score, 1, "geyer")[0]
    except ValueError:
        print(
            "Geyer equilibriation test failed. See the traceback below:"
        )
        traceback.print_exc()
        t_eq = None

    return t_eq

def multiscale_equilibriation_test(score):
    """ Perform Multiscale equilibration test on the given score trajectory.

    Args:
        score (list): Score values.

    Returns:
        int: Equilibration time index.
    """
    try:
        t_eq = detectEquilibration(score, 1, "multiscale")[0]
    except ValueError:
        print(
            "Multiscale equilibriation test failed. See the traceback below:"
        )
        traceback.print_exc()
        t_eq = None

    return t_eq

def get_restraint_wise_score_trajectory(
    per_frame_stats,
    inv_head_dict,
    n_remove,
    exchange_indices,
    n,
    restraints,
    eq,
):
    all_res_info = dict()
    t_eq = None
    t_eq_2 = None

    for restraint in restraints:
        reskey, resname = restraint.split(":")
        score = parse_key(
                reskey,
                per_frame_stats,
                inv_head_dict,
                exchange_indices,
                n_remove,
                adjust=False,
                exact_match=False,
                return_type="sum",
            )

        if eq:
            t_eq = geyer_equilibriation_test(score)
            t_eq_2 = multiscale_equilibriation_test(score)

        fitline = sm.nonparametric.lowess(
            score, np.arange(n_remove, n), frac=0.2, return_sorted=False
        )

        all_res_info[resname] = (
            fitline - fitline.max(),
            get_moving_sd(score),
            t_eq,
            t_eq_2,
        )

    return all_res_info

def plot_restraint_wise_score_trajectory(
    savepath,
    n_remove,
    n,
    all_res_info,
):
    for resname in all_res_info:
        plt.figure()
        plt.plot(np.arange(n_remove, n), all_res_info[resname][0], color="red", lw=3)
        plt.fill_between(
                np.arange(n_remove, n),
                all_res_info[resname][0] + all_res_info[resname][1],
                all_res_info[resname][0] - all_res_info[resname][1],
                alpha=0.5,
                zorder=-10,
                color="black",
            )

        legend_labels = ["Restraint Score (shifted)", "Moving SD"]

        if not (all_res_info[resname][2] is None):
            plt.axvline(
                n_remove + all_res_info[resname][2],
                0,
                1,
                color="green",
                lw=2,
            )
            legend_labels.append("Geyer Equilibration Point")

        if not (all_res_info[resname][3] is None):
            plt.axvline(
                    n_remove + all_res_info[resname][3],
                    0,
                    1,
                    color="green",
                    lw=2,
                    linestyle=":",
                )
            legend_labels.append("Multiscale Equilibration Point")

        plt.xlabel("Frame Index")
        plt.ylabel("Restraint Score (shifted)")
        plt.title(f"Restraint: {resname}")
        plt.legend(legend_labels, loc="upper right")
        plt.savefig(f"{savepath}/equilibriation_{resname}.png")
        plt.close("all")

if __name__ == "__main__":
    start_time = time.time()
    args = ArgumentParser()
    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to the directory containing stat and stat_replica files",
    )
    args.add_argument(
        "-b",
        "--burn_in",
        type=float,
        default=0.9,
        help="Fraction of frames to discard as burn-in. Default is 0.1",
    )
    args.add_argument(
        "-r",
        "--run_type",
        type=str,
        default="replica",
        help="Type of run: replica or single. Default is replica",
        choices=["replica", "single"]
    )
    args.add_argument(
        "-c",
        "--num_proc",
        type=int,
        default=4,
        help="Number of processes to use. Do not use too many, you might run out of memory.",
    )
    args.add_argument(
        "-s",
        "--sysname",
        type=str,
        default="cardiac_desmosome",
        help="System name for specific analysis. Default is cardiac_desmosome",
    )
    args.add_argument(
        "--restraints",
        nargs="+",
        type=str,
        default=["ExcludedVolumeSphere:EVR", "GaussianEMRestraint:EMR"],
        help="List of restraint names for specific analysis.",
    )
    args.add_argument(
        "--detect_equilibration",
        action="store_true",
        help="Whether to perform equilibration detection.",
    )
    args = args.parse_args()

    if args.run_type == "single":

        result = parser(args.input)

        (
            per_frame_stats,
            per_frame_replica_stats,
            stack_list,
            replica_stack_list,
            inv_head_dict,
            inv_replica_head_dict,
            frame_order,
        ) = result

        avg_stats_df = avg_stats(
            per_frame_stats,
            frame_order,
            inv_head_dict,
            burn_in=args.burn_in,
        )

        avg_stats_df.to_csv(
            os.path.join(args.input, "average_stats.csv"),
            header=True,
            index=False
        )

        n_remove = int(args.burn_in * len(per_frame_stats))
        exchange_indices = sort_the_replica_exchanges_lowest_temp(frame_order)[1]
        tot_score = [float(j[inv_head_dict["Total_Score"]]) for j in per_frame_stats]

        n = len(per_frame_stats)
        n_keep = n - n_remove

        plot_trajectory_score(
            args.input,
            per_frame_stats,
            n_remove,
            exchange_indices,
            tot_score
        )
        exchange_counts_all = sort_the_replica_exchanges_all_temp(
            replica_stack_list,
            inv_replica_head_dict
        )

        restraints = args.restraints
        eq = True if args.detect_equilibration else False

        all_res_info = get_restraint_wise_score_trajectory(
            per_frame_stats,
            inv_head_dict,
            n_remove,
            exchange_indices,
            n,
            restraints,
            eq,
        )

        plot_restraint_wise_score_trajectory(
            args.input,
            n_remove,
            n,
            all_res_info
        )

    elif args.run_type == "replica":

        from concurrent.futures import ProcessPoolExecutor, as_completed
        from tqdm import tqdm

        traj_dirs = os.listdir(args.input)
        traj_dirs = [
            os.path.join(args.input, d) for d in traj_dirs
            if os.path.isdir(os.path.join(args.input, d))
        ]
        traj_dirs = sorted(traj_dirs)

        num_proc = np.min([args.num_proc, len(traj_dirs)])

        results = []
        with ProcessPoolExecutor(max_workers=num_proc) as executor:
            futures = {
                executor.submit(
                    prod_run_worker,
                    run_dir,
                    args.burn_in
                ): run_dir for run_dir in traj_dirs
            }

            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Processing trajectories"
            ):
                run_dir = futures[future]
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    print(f"Error processing {run_dir}: {e}")

        # average across all trajectories
        all_avg_stats_df = pd.concat(results)
        all_avg_stats_df = all_avg_stats_df.groupby("Key").mean().reset_index()
        all_avg_stats_df.to_csv(
            os.path.join(args.input, "average_stats_all_traj.csv"),
            header=True,
            index=False
        )

    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.2f} seconds")
