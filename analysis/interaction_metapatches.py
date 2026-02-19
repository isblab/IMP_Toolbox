import os
import re
import getpass
import pandas as pd
import matplotlib.pyplot as plt
from pprint import pprint
from argparse import ArgumentParser
from plotly import graph_objects as go
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import (
    get_res_range_from_key,
    get_key_from_res_range
)

_user = getpass.getuser()

def parse_patch_csv_filename(csv: str):
    """ Parse the patch CSV filename to extract the protein names and residue ranges for both proteins.
    The expected filename format is: patches_{p1_name}:{p1_start}-{p1_end}_{p2_name}:{p2_start}-{p2_end}.csv
    or patches_{p1_name}:{p1_start}_{p2_name}:{p2_start}.csv

    e.g.: patches_PKP2a:1-837_DSC2a:720-901.csv, patches_PKP2a:1-837_DSC2a:901.csv

    ## Arguments:

    - **csv (str)**:<br />
        Path to the patch CSV file.
        The filename should follow the expected format to extract the necessary information.

    ## Returns:

    - **tuple**:<br />
        Tuple of four strings: (p1_name, p1_range_full, p2_name, p2_range_full)
    """

    groups = re.match(
        r"patches_(\w+):([\d\-]+)_(\w+):([\d\-]+)\.csv",
        os.path.basename(csv)
    ).groups()

    if not groups or len(groups) != 4:
        raise ValueError(f"Filename {csv} does not match the expected format.")

    p1_name = groups[0]
    p1_range_full = groups[1]

    p2_name = groups[2]
    p2_range_full = groups[3]

    return p1_name, p1_range_full, p2_name, p2_range_full

def extract_interacting_patches_info(csv: str):
    """ Extract the following information from the patches csv file.
    - protein names
    - residues ranges (start, end)

    ## Arguments:

    - **csv (str)**:<br />
        Path to the patch CSV file.
        The filename should follow the expected format to extract the necessary information.

    ## Returns:

    - **dict**:<br />
        Dictionary containing the extracted information:
        - p1_name: name of the first protein
        - p1_start: start residue of the first protein
        - p1_end: end residue of the first protein
        - p2_name: name of the second protein
        - p2_start: start residue of the second protein
        - p2_end: end residue of the second protein
    """

    _p1_name, p1_range_full, _p2_name, p2_range_full = parse_patch_csv_filename(csv)

    p1_range_ = get_res_range_from_key(p1_range_full)
    p2_range_ = get_res_range_from_key(p2_range_full)

    p1_start, p1_end = min(p1_range_), max(p1_range_)
    p2_start, p2_end = min(p2_range_), max(p2_range_)

    # Analyze the patch CSV to identify interacting metapatches
    df = pd.read_csv(csv)
    p1_name = df.columns[0].split("_")[0]
    p2_name = df.columns[1].split("_")[0]

    # swap if reversed in the filename vs the CSV columns
    if _p1_name == p2_name and _p2_name == p1_name:
        p1_range_full, p2_range_full = p2_range_full, p1_range_full
        p1_start, p1_end, p2_start, p2_end = p2_start, p2_end, p1_start, p1_end

    info_dict = {
        "p1_name": p1_name,
        "p1_start": p1_start,
        "p1_end": p1_end,
        "p2_name": p2_name,
        "p2_start": p2_start,
        "p2_end": p2_end,
    }

    return df, info_dict

def extract_metapatches(p1_res: set, p2_res: set, xy_data: pd.DataFrame):
    """ Extract interacting metapatches.

    ## Arguments:

    - **p1_res (set)**:<br />
        Interface residues for protein 1 + any filled gaps.

    - **p2_res (set)**:<br />
        Interface residues for protein 2 + any filled gaps.

    - **xy_data (pd.DataFrame)**:<br />
        DataFrame containing the interacting residue pairs from the contact map.
        It should have two columns:
            "x" for residues of protein 1
            "y" for residues of protein 2

    ## Returns:

    - **list**:<br />
        List of tuples containing the interacting metapatches (p1_patch, p2_patch)
    """

    # update the patches based on the filled gaps
    p1_patches = get_key_from_res_range(list(p1_res), as_list=True)
    p2_patches = get_key_from_res_range(list(p2_res), as_list=True)

    possible_metapatches = []
    for p1_patch in p1_patches:
        for p2_patch in p2_patches:
            possible_metapatches.append((p1_patch, p2_patch))

    interacting_metapatches = []

    for p1_patch, p2_patch in possible_metapatches:

        patch1_start, patch1_end = min(get_res_range_from_key(p1_patch)), max(get_res_range_from_key(p1_patch))
        patch2_start, patch2_end = min(get_res_range_from_key(p2_patch)), max(get_res_range_from_key(p2_patch))

        # check if any of the points in the contact map fall within this patch region
        in_patch = xy_data[
            (xy_data["x"] >= patch1_start) & (xy_data["x"] <= patch1_end) &
            (xy_data["y"] >= patch2_start) & (xy_data["y"] <= patch2_end)
        ]
        if in_patch.empty:
            continue

        # reduce the patch region to tightly fit the points in the contact map that fall within it
        patch1_start = in_patch["x"].min()
        patch1_end = in_patch["x"].max()
        patch2_start = in_patch["y"].min()
        patch2_end = in_patch["y"].max()

        interacting_metapatches.append((f"{patch1_start}-{patch1_end}", f"{patch2_start}-{patch2_end}"))

    return interacting_metapatches

def fill_gaps(
    threshold: int,
    p_gaps: list,
    p_res: set,
):
    """ Fill gaps within the list of residues based on threshold.

    ## Arguments:

    - **threshold (int)**:<br />
        Maximum gap size (in residues) to fill when identifying metapatches.

    - **p_gaps (list)**:<br />
        List of gap ranges (as strings like "10-15") to be filled.

    - **p_res (set)**:<br />
        Set of residues that are currently identified as part of patches, to which the filled gaps will be added.

    ## Returns:

    - **set**:<br />
        Updated set of residues with the filled gaps added.
    """

    for gap in p_gaps:
        gap_start, gap_end = min(get_res_range_from_key(gap)), max(get_res_range_from_key(gap))
        if len(list(range(gap_start, gap_end+1))) <= threshold:
            p_res.update(list(range(gap_start, gap_end+1)))

    return p_res

def identify_gaps(
    p_start: int,
    p_end: int,
    p1_res: set
):
    """ Identify the gaps within a set of residues.

    ## Arguments:

    - **p_start (int)**:<br />
        Start residue number of the protein.

    - **p_end (int)**:<br />
        End residue number of the protein.

    - **p1_res (set)**:<br />
        Set of residues that are currently identified as part of patches.

    ## Returns:

    - **list**:<br />
        List of gap ranges (as strings like "10-15") that are within the protein's range but not part of any patch.
    """

    p_gaps = set(range(p_start, p_end+1)) - p1_res
    p_gaps = get_key_from_res_range(list(p_gaps), as_list=True)

    # remove start and end gaps
    p_gaps = [gap for gap in p_gaps if not (min(get_res_range_from_key(gap)) == p_start or max(get_res_range_from_key(gap)) == p_end)]

    return p_gaps

def get_p1_p2_patch_residues(patches_df: pd.DataFrame):
    """ Extract interacting residues from both proteins from patches dataframe.

    ## Arguments:

    - **patches_df (pd.DataFrame)**:<br />
        DataFrame containing the interacting patches information.
        It should have two columns corresponding to the two proteins, with patch ranges as values.

    ## Returns:

    - **tuple**:<br />
        A tuple of two sets,
        - interacting residues of the first protein
        - interacting residues of the second protein
    """

    p1_patches = set(patches_df[patches_df.columns[0]].dropna().tolist())
    p2_patches = set(patches_df[patches_df.columns[1]].dropna().tolist())

    p1_res = set()
    p2_res = set()

    for res_range in p1_patches:
        p1_res.update(get_res_range_from_key(res_range))

    for res_range in p2_patches:
        p2_res.update(get_res_range_from_key(res_range))

    return p1_res, p2_res

def get_interacting_residue_pairs(patches_df: pd.DataFrame):
    """ Extract interacting residue pairs from patches dataframe.

    ## Arguments:

    - **patches_df (pd.DataFrame)**:<br />
        DataFrame containing the interacting patches information.
        It should have two columns corresponding to the two proteins, with patch ranges as values.

    ## Returns:

    - **pd.DataFrame**:<br />
        DataFrame with columns "x" and "y" representing interacting residue pairs.
    """

    xy_data = []
    for _, row in patches_df.iterrows():
        if pd.isna(row[patches_df.columns[0]]) or pd.isna(row[patches_df.columns[1]]):
            continue

        x_vals = get_res_range_from_key(row[patches_df.columns[0]])
        y_vals = get_res_range_from_key(row[patches_df.columns[1]])

        for x in x_vals:
            for y in y_vals:
                xy_data.append((x, y))

    xy_data = pd.DataFrame(xy_data, columns=["x", "y"])

    return xy_data

def plot_metapatches(
    xy_data: pd.DataFrame,
    interacting_metapatches: list,
    info_dict: dict,
    out_dir: str,
    plotting_library: str = "matplotlib",
):
    """ Plot the metapatches.

    ## Arguments:

    - **xy_data (pd.DataFrame)**:<br />
        DataFrame with columns "x" and "y" representing interacting residue pairs.

    - **interacting_metapatches (list)**:<br />
        List of tuples, each tuple containing two strings representing patch ranges for the interacting metapatches.

    - **info_dict (dict)**:<br />
        Dictionary containing the extracted information about the proteins and their residue ranges:
        - p1_name: name of the first protein
        - p1_start: start residue of the first protein
        - p1_end: end residue of the first protein
        - p2_name: name of the second protein
        - p2_start: start residue of the second protein
        - p2_end: end residue of the second protein

    - **out_dir (str)**:<br />
        Directory where the plot will be saved.

    - **plotting_library (str, optional):**:<br />
        Plotting library to use for visualizations.
        Default is "matplotlib". Other option is "plotly".
    """

    p1_name = info_dict["p1_name"]
    p2_name = info_dict["p2_name"]
    p1_start = info_dict["p1_start"]
    p1_end = info_dict["p1_end"]
    p2_start = info_dict["p2_start"]
    p2_end = info_dict["p2_end"]

    file_name = f"metapatches_{p1_name}:{p1_start}-{p1_end}_{p2_name}:{p2_start}-{p2_end}"

    if plotting_library == "matplotlib":

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(xy_data["x"], xy_data["y"], color='blue', alpha=0.5, s=10)

        for p1_patch, p2_patch in interacting_metapatches:
            p1_start_, p1_end_ = min(get_res_range_from_key(p1_patch)), max(get_res_range_from_key(p1_patch))
            p2_start_, p2_end_ = min(get_res_range_from_key(p2_patch)), max(get_res_range_from_key(p2_patch))
            width = p1_end_ - p1_start_
            height = p2_end_ - p2_start_
            if width > 0 and height > 0:
                ax.add_patch(plt.Rectangle((p1_start_, p2_start_), width, height, color='red', alpha=0.3))

        ax.set_xlim(p1_start, p1_end)
        ax.set_ylim(p2_start, p2_end)

        ax.set_title(f"Interacting Metapatches for {p1_name} and {p2_name}")
        ax.set_xlabel(f"{p1_name} Residue Number")
        ax.set_ylabel(f"{p2_name} Residue Number")

        plt.savefig(os.path.join(out_dir, f"{file_name}.png"), dpi=300)
        plt.close()

    elif plotting_library == "plotly":

        fig = go.Figure()

        for p1_patch, p2_patch in interacting_metapatches:
            p1_start_, p1_end_ = min(get_res_range_from_key(p1_patch)), max(get_res_range_from_key(p1_patch))
            p2_start_, p2_end_ = min(get_res_range_from_key(p2_patch)), max(get_res_range_from_key(p2_patch))
            width = p1_end_ - p1_start_
            height = p2_end_ - p2_start_
            if width > 0 and height > 0:
                fig.add_shape(
                    type="rect",
                    x0=p1_start_,
                    y0=p2_start_,
                    x1=p1_end_,
                    y1=p2_end_,
                    line=dict(color='red', width=2),
                    fillcolor='rgba(255, 0, 0, 0.3)',
                )

                # hove text for the metapatch
                fig.add_trace(
                    go.Scatter(
                        x=[p1_start_, p1_start_, p1_end_, p1_end_, p1_start_],
                        y=[p2_start_, p2_end_, p2_end_, p2_start_, p2_start_],
                        mode='lines',
                        fill='toself',
                        name='',
                        opacity=0,
                        text=f"Metapatch: {p1_name}:{p1_patch} - {p2_name}:{p2_patch}",
                        hoverinfo='text',
                    )
                )

        fig.add_trace(
            go.Scatter(
                x=xy_data["x"],
                y=xy_data["y"],
                mode='markers',
                marker=dict(color='blue', opacity=0.5, size=5),
                text=[f"Residue {x} ({p1_name}) - Residue {y} ({p2_name})" for x, y in zip(xy_data["x"], xy_data["y"])],
                name='Interacting residue pairs',
                hovertemplate='%{text}<extra></extra>',
                showlegend=True
            )
        )
        fig.update_layout(
            title=f"Interacting Metapatches for {p1_name} and {p2_name}",
            xaxis_title=f"{p1_name} Residue Number",
            yaxis_title=f"{p2_name} Residue Number",
            xaxis=dict(range=[p1_start, p1_end]),
            yaxis=dict(range=[p2_start, p2_end]),
        )

        fig.write_html(os.path.join(out_dir, f"{file_name}.html"))

if __name__ == "__main__":

    parser = ArgumentParser(description="Extract interacting metapatches.")

    parser.add_argument(
        "--interaction_map_dir",
        type=str,
        default=f"/data/{_user}/Projects/cardiac_desmosome/analysis/prod_runs/analysis_cis_dsc_5716/set4_analysis3/contact_map",
        help="Path to the contact map output directory"
    )

    parser.add_argument(
        "--threshold",
        type=int,
        default=40,
        help="Maximum gap size (in residues) to fill when identifying metapatches"
    )

    parser.add_argument(
        "--plotting",
        type=str,
        default="plotly",
        choices=["matplotlib", "plotly", "no_plot"],
        help="Plotting library to use for visualizations"
    )

    args = parser.parse_args()

    interaction_map_dir = args.interaction_map_dir
    threshold = args.threshold

    patches_csvs = [
        os.path.join(interaction_map_dir, f)
        for f in os.listdir(interaction_map_dir)
        if f.startswith("patches_") and f.endswith(".csv")
    ]

    for _, csv in enumerate(patches_csvs):

        patches_df, info_dict = extract_interacting_patches_info(csv)
        p1_start = info_dict["p1_start"]
        p1_end = info_dict["p1_end"]
        p2_start = info_dict["p2_start"]
        p2_end = info_dict["p2_end"]
        p1_name = info_dict["p1_name"]
        p2_name = info_dict["p2_name"]

        # extract interacting residue pairs from the CSV
        xy_data = get_interacting_residue_pairs(patches_df)

        # identify gaps in the patches and fill them if they are below the threshold size
        p1_res, p2_res = get_p1_p2_patch_residues(patches_df)

        p1_gaps = identify_gaps(p1_start, p1_end, p1_res)
        p2_gaps = identify_gaps(p2_start, p2_end, p2_res)

        p1_res = fill_gaps(threshold, p1_gaps, p1_res)
        p2_res = fill_gaps(threshold, p2_gaps, p2_res)

        interacting_metapatches = extract_metapatches(p1_res, p2_res, xy_data)

        # pprint(interacting_metapatches)
        # print("-"*50)

        if args.plotting == "no_plot":
            continue

        plot_metapatches(
            xy_data=xy_data,
            interacting_metapatches=interacting_metapatches,
            info_dict=info_dict,
            out_dir=interaction_map_dir,
            plotting_library=args.plotting
        )