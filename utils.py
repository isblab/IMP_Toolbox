from collections import defaultdict
from Bio.PDB import Select
import Bio
import json
import numpy as np
import pandas as pd
import requests
import platform
import os

def get_dropbox_path():
    system_spec_paths = {
    "Windows": f"C:\\Users\\{os.getlogin()}\\Dropbox",
    "Linux": f"/home/{os.getlogin()}/Dropbox"
}

    if platform.system() in system_spec_paths:
        path_to_dropbox = system_spec_paths[platform.system()]
    else:
        path_to_dropbox = input("Enter the path to your Dropbox folder: ")
    return path_to_dropbox

def request_session(max_retries=3):
    """Create a request session with set max retries

    Args:
        max_retries (int, optional): Defaults to 3.

    Returns:
        req_sess (requests.sessions.Session): request session
    """
    req_sess = requests.Session()
    req_sess.mount(
        "https://",
        requests.adapters.HTTPAdapter(max_retries=max_retries)
    )
    return req_sess

def request_result(get_request, uniprot_id, ignore_error=False):
    """Get the result of a get request

    Args:
        get_request (requests.models.Response): get request
        uniprot_id (str): valid entrypoint id
        ignore_error (bool, optional): Defaults to False.

    """
    if get_request.status_code == 200:
        try:
            return get_request.json()
        except json.decoder.JSONDecodeError:
            return get_request.content
    else:
        print(f"Error while requesting {uniprot_id}") if not ignore_error else None
        return None

def write_json(file_path, data):
    """Write data to a json file

    Args:
        file_path (str): path to json file
        data (dict): data to write
    """
    with open(file_path, "w") as f:
        json.dump(data, f)

def read_json(file_path):
    """Load a json file

    Args:
        file_path (str): path to json file

    Returns:
        data (dict): data from json file
    """
    with open(file_path, "r") as f:
        data = json.load(f)
    return data


def read_fasta(fasta_file):
    """
    Read a fasta file and return a dictionary of sequences
    """

    all_sequences = {}

    with open(fasta_file, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith(">"):
            seq_id = line[1:].strip()
        else:
            seq = line.strip()
            all_sequences[seq_id] = seq if seq_id not in all_sequences else all_sequences[seq_id] + seq

    return all_sequences

##%%

def get_key_from_res_range(res_range, as_list=False):
    """Returns a residue range string from a list of residue numbers.

    Args:
        res_range (list): List of residue numbers, e.g., [1, 2, 3, 5, 6, 7]

    Returns:
        str: Residue range string, e.g., "1-3,5-7"
    """

    if not res_range:
        return ""
    res_range = sorted(res_range)
    ranges = []
    start = prev = res_range[0]
    for num in res_range[1:]:
        if num == prev + 1:
            prev = num
        else:
            ranges.append(f"{start}-{prev}") if start != prev else ranges.append(str(start))
            start = prev = num
    if start == prev:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{prev}")
    if as_list:
        return ranges
    else:
        return ",".join(ranges)

def get_ones(matrix, axis=0):
    """Get the indices of 1s in a matrix rowwise or columnwise
    example: matrix = np.array(
        [1, 0, 1],
        [0, 1, 0],
        [1, 1, 0]
    )
    get_ones(matrix, axis=0) -> {0: {0, 2}, 1: {1}, 2: {0, 1}}
    get_ones(matrix, axis=1) -> {0: {0, 2}, 1: {1, 2}, 2: {0}}

    Args:
        matrix (_type_): _description_
        axis (int, optional): _description_. Defaults to 0.

    Returns:
        _type_: _description_
    """    
    one_sets = {}
    if axis == 0:
        for i in range(matrix.shape[0]):
            one_sets[i] = set(np.where(matrix[i] == 1)[0])
    elif axis == 1:
        for j in range(matrix.shape[1]):
            one_sets[j] = set(np.where(matrix[:, j] == 1)[0])
    return one_sets

def split_sets(one_sets):
    """Split the sets of 1s into subsets
    example:
    one_sets = {0: {0, 2, 3, 5, 6}, 1: {1}, 2: {0, 1}}
    split_sets(one_sets) -> {0: [[0, 2], [3], [5, 6]], 1: [[1]], 2: [[0, 1]]}

    Args:
        one_sets (dict): dictionary of sets of 1s

    Returns:
        _type_: _description_
    """
    new_one_sets = {}
    for i, one_set in one_sets.items():
        sub_sets = []
        one_set = sorted(list(one_set))
        for idx, val in enumerate(one_set):
            if idx == 0:
                sub_sets.append([val])
            elif val - one_set[idx-1] == 1:
                sub_sets[-1].append(val)
            else:
                sub_sets.append([val])
        new_one_sets[i] = sub_sets
    return new_one_sets

# split_sets(
#     {
#         0: {0, 2, 3, 5, 6},
#         1: {1},
#         2: {0, 1}
#     }
# )

def remove_duplicate_sets(list_of_sets):
    """Remove duplicate sets from a list of sets
    example:
    list_of_sets = [{0, 1, 2}, {0, 1, 2}, {1, 2, 3}]
    remove_duplicate_sets(list_of_sets) -> [{0, 1, 2}, {1, 2, 3}]

    Args:
        list_of_sets (_type_): _description_

    Returns:
        _type_: _description_
    """    
    
    new_list_of_sets = []
    for set1 in list_of_sets:
        new_list_of_sets.append(set1) if set1 not in new_list_of_sets else None
    return new_list_of_sets

# remove_duplicate_sets([{0, 1, 2}, {0, 1, 2}, {1, 2, 3}])



def add_subsets_of_other_sets(one_sets, list_of_sets):
    """Add the subsets of the sets in list_of_sets to the one_sets
    example:
    one_sets = {0: {0, 1, 2, 3, 4, 5, 6}, 1: {1}, 2: {0, 1}}
    list_of_sets = [{0, 1}, {1}, {0, 1, 2, 3, 4, 5, 6}]
    add_subsets_of_other_sets(one_sets, list_of_sets) -> {0: [{0, 1}, {0, 1, 2, 3, 4, 5, 6}, {1}], 1: [{1}], 2: [{0, 1}, {1}]}

    Args:
        one_sets (_type_): _description_
        list_of_sets (_type_): _description_
    """
    from collections import defaultdict
    new_one_sets = defaultdict(list)
    for set1 in list_of_sets:
        for idx, one_set in one_sets.items():
            if set1.issubset(one_set):
                new_one_sets[idx].append(set1) if set1 not in new_one_sets[idx] else None
    return new_one_sets

# add_subsets_of_other_sets(
#     {
#         0: {0, 1, 2, 3, 4, 5, 6},
#         1: {1},
#         2: {0, 1}
#     },
#     [{0, 1}, {1}, {0, 1, 2, 3, 4, 5, 6}]
# )



def special_dict_to_df(dict_, columns):
    """Convert a dictionary to a pandas DataFrame
    example:
    dict_ = {
        1: [{1, 2}, {5}],
        2: [{4, 5}, {6}]
        }
    columns = ["A", "B"]
    special_dict_to_df(dict_, columns) ->
    A   B
    1   {1, 2}
    1   {5}
    2   {4, 5}
    2   {6}

    Args:
        dict_ (dict): dictionary to convert
        columns (list): column names

    Returns:
        pd.DataFrame: DataFrame from the dictionary
    """
    import pandas as pd
    if all([isinstance(val, list) for val in dict_.values()]):
        df = pd.DataFrame(dict_, columns=columns)
        df_rows = []
        for k, v in dict_.items():
            for val in v:
                df_rows.append([str(k), val])
        df = pd.DataFrame(df_rows, columns=columns)
    else:
        raise ValueError("All values in the dictionary must be lists.")
    return df

# special_dict_to_df(
#     {
#         1: [{1, 2}, {5}],
#         2: [{4, 5}, {6}]
#     },
#     ["A", "B"]
# )



def df_set_to_res_range(df, cols_to_convert):
    """Convert sets in a DataFrame to residue ranges
    example:
    df = pd.DataFrame({
        "A": [1, 1, 2, 2],
        "B": [{1, 2}, {5}, {4, 5}, {6}]
    })
    df_set_to_res_range(df, "B") ->
    A   B
    1   1-2
    1   5
    2   4-5
    2   6

    Args:
        df (_type_): _description_
        groupby_col (_type_): _description_

    Returns:
        _type_: _description_
    """    
    for col in cols_to_convert:
        for idx, row in df.iterrows():
            one_set = row[col]
            if len(one_set) == 1:
                one_set = str(list(one_set)[0])
            else:
                one_set = sorted(list(one_set))
                one_set = f"{one_set[0]}-{one_set[-1]}"
            df.at[idx, col] = one_set
    return df

# df = special_dict_to_df(
#     {
#         1: [{1, 2}, {5}],
#         2: [{4, 5}, {6}]
#     },
#     ["A", "B"]
# )
# df_set_to_res_range(df, ["B"])



def groupby_df(df, groupby_col, agg_col):
    """Group a DataFrame by a column and aggregate another column
    example:
    df = pd.DataFrame({
        "A": [1, 1, 2, 2, 3, 4],
        "B": ['1-2', '5', '4-5', '6', '1-2', '1-2']
    })
    groupby_df(df, "B", "A") ->
    A       B
    [1,3,4] 1-2
    [2]     4-5
    [1]     5
    [2]     6

    Args:
        df (_type_): _description_
        groupby_col (_type_): _description_
        agg_col (_type_): _description_

    Returns:
        _type_: _description_
    """    
    df_group = df.groupby([groupby_col])[agg_col].apply(','.join).reset_index()
    for idx, row in df_group.iterrows():
        one_set = row[agg_col].split(",")
        one_set = [int(x) for x in one_set]
        one_set = sorted(one_set)
        df_group.at[idx, agg_col] = one_set
    return df_group

# df = special_dict_to_df(
#     {
#         1: [{1, 2}, {5}],
#         2: [{4, 5}, {6}],
#         3: [{1, 2}],
#         4: [{1, 2}]
#     },
#     ["A", "B"]
# )
# df = df_set_to_res_range(df, "B")
# groupby_df(df, "B", "A")


def combine_dfs(df1, df2, chain1, chain2):
    combined_df = pd.concat([df2, df1], axis=0)
    combined_df.reset_index(drop=True, inplace=True)

    new_df = pd.DataFrame(columns=[chain1, chain2])
    df_rows = []

    for idx, row in combined_df.iterrows():
        if isinstance(row[chain1], list):
            for res_range in get_key_from_res_range(row[chain1], as_list=True):
                df_rows.append([res_range, row[chain2]])
        elif isinstance(row[chain2], list):
            for res_range in get_key_from_res_range(row[chain2], as_list=True):
                df_rows.append([row[chain1], res_range])

    new_df = pd.DataFrame(df_rows, columns=[chain1, chain2])
    new_df.drop_duplicates(inplace=True, keep="first")
    new_df.reset_index(drop=True, inplace=True)
    return new_df


def df_res_range_to_set(df, cols_to_convert):
    for col in cols_to_convert:
        for id, row in df.iterrows():
            a_range = row[col].split("-")
            a_range = set(list(range(int(a_range[0]), int(a_range[1])+1))) if len(a_range) > 1 else set([int(a_range[0])])
            df.at[id, col] = a_range

    return df

def remove_subset_rows(df, chain1, chain2):
    del_idx = []
    from copy import deepcopy
    for idx1, row1 in df.iterrows():
        a_range1 = df.at[idx1, chain1]
        b_range1 = df.at[idx1, chain2]
        rem_df = deepcopy(df)
        rem_df.drop(index=idx1, inplace=True)
        for idx2, row2 in rem_df.iterrows():
            a_range2 = df.at[idx2, chain1]
            b_range2 = df.at[idx2, chain2]
            if a_range1.issubset(a_range2) and b_range1.issubset(b_range2):
                del_idx.append(idx1) if idx1 not in del_idx else None

    for idx in del_idx:
        df.drop(idx, inplace=True)
    return df

def get_patches_from_matrix(matrix, chain1, chain2):
    """Get all interacting patches from a binary matrix

    Args:
        matrix (_type_): _description_
        chain1 (_type_): _description_
        chain2 (_type_): _description_

    Returns:
        _type_: _description_
    """    
    row_sets = get_ones(matrix, axis=0)
    col_sets = get_ones(matrix, axis=1)
    # print("row_sets", row_sets)
    # print("col_sets", col_sets)

    split_row_sets = split_sets(row_sets)
    split_col_sets = split_sets(col_sets)
    # print("split_row_sets", split_row_sets)
    # print("split_col_sets", split_col_sets)

    all_row_sets = [set(x) for xs in split_row_sets.values() for x in xs]
    all_col_sets = [set(x) for xs in split_col_sets.values() for x in xs]
    # print("all_row_sets", all_row_sets)
    # print("all_col_sets", all_col_sets)

    all_row_sets = remove_duplicate_sets(all_row_sets)
    all_col_sets = remove_duplicate_sets(all_col_sets)
    # print("all_row_sets", all_row_sets)
    # print("all_col_sets", all_col_sets)

    split_row_sets = add_subsets_of_other_sets(row_sets, all_row_sets)
    split_col_sets = add_subsets_of_other_sets(col_sets, all_col_sets)
    # print("split_row_sets", split_row_sets)
    # print("split_col_sets", split_col_sets)

    df_row = special_dict_to_df(split_row_sets, [chain1, chain2])
    df_col = special_dict_to_df(split_col_sets, [chain2, chain1])
    # print("df_row", df_row)
    # print("df_col", df_col)

    # now we have to use groupby to group the patches
    # so convert sets to strings
    df_row = df_set_to_res_range(df_row, [chain2])
    df_col = df_set_to_res_range(df_col, [chain1])
    # print("df_row", df_row)
    # print("df_col", df_col)

    df_row = groupby_df(df_row, chain2, chain1)
    df_col = groupby_df(df_col, chain1, chain2)

    combined_df = combine_dfs(df_row, df_col, chain1, chain2)
    combined_df = df_res_range_to_set(combined_df, [chain1, chain2])
    combined_df = remove_subset_rows(combined_df, chain1, chain2)
    combined_df = df_set_to_res_range(combined_df, [chain1, chain2])

    return combined_df
# import numpy as np
# binary_matrix = np.array(
#     [
#         [0, 0, 1, 1, 1],
#         [0, 1, 1, 0, 1],
#         [0, 1, 1, 0, 0],
#         [0, 0, 0, 1, 1],
#         [0, 0, 0, 1, 1]
#     ]
# )
# print(binary_matrix)
# get_patches_from_matrix(binary_matrix, "A", "B")
##%%

def save_df(df, file_path, index=False, header=True, sep=","):
    """Save a pandas DataFrame to a file

    Args:
        df (pandas.DataFrame): DataFrame to save
        file_path (str): path to save the DataFrame
    """
    df.to_csv(file_path, index=index, header=header, sep=sep)

def generate_cmap(n, scheme="soft-warm"):
    colors = set()  # Use a set to prevent duplicates
    import random
    import time
    start = time.time()
    while len(colors) < n:
        if scheme == "standard":
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = random.randint(0, 255)

        elif scheme == "non-bright":
            r = random.randint(50, 180)
            g = random.randint(50, 180)
            b = random.randint(50, 180)

        elif scheme == "earth-tone":
            r = random.randint(100, 180)
            g = random.randint(60, 140)
            b = random.randint(40, 120)

        elif scheme == "cool-tone":
            r = random.randint(50, 120)
            g = random.randint(100, 180)
            b = random.randint(120, 255)

        elif scheme == "soft-warm":
            r = random.randint(180, 255)
            g = random.randint(130, 200)
            b = random.randint(90, 160)

        elif scheme == "contrasting-non-bright":
            if len(colors) % 2 == 0:  # Alternate between darker and lighter muted tones
                r = random.randint(40, 120)
                g = random.randint(40, 120)
                b = random.randint(40, 120)
            else:
                r = random.randint(140, 200)
                g = random.randint(140, 200)
                b = random.randint(140, 200)

        else:
            raise ValueError("Invalid scheme. Choose from 'non-bright', 'earth-tone', 'cool-tone', or 'soft-warm'.")

        color = "#{:02x}{:02x}{:02x}".format(r, g, b)
        colors.add(color)  # Add to set (automatically avoids duplicates)

        if time.time() - start > 10:
            break

    return list(colors)  # Convert set back to list before returning

# def generate_cmap(n):
#     """
#     Generate a custom colormap with n colors.
#     """

#     import random
#     import time
#     colors = []

#     r = random.randint(50, 180)  # Avoid very bright reds
#     g = random.randint(50, 180)  # Avoid very bright greens
#     b = random.randint(50, 180)  # Avoid very bright blues

#     start = time.time()

#     while len(colors) < n:

#         # color = "%06x" % random.randint(0, 0xFFFFFF)
#         color = "#{:02x}{:02x}{:02x}".format(r, g, b)

#         if f"#{color}" not in colors:
#             colors.append(f"#{color}")

#         if time.time() - start > 10:
#             break

#     return colors

def str_join( prots: list, sep: str ):
    """
    Given a list of strings, join them by the provided separator (sep).
    """
    return f"{sep}".join( prots )


def str_split( prot_pair: str, sep: str ):
    """
    Split the given string by the  provided separator (sep).
    """
    return prot_pair.split( sep )

def plot_map(contact_map: np.array, chain1: str, chain2: str, p1_region: tuple, p2_region: tuple, plot_type: str):
    """Plot the contact map

    Args:
        contact_map (np.array): binary contact map
    """

    # p1_region, p2_region = interacting_region.values()
    # chain1, chain2 = interacting_region.keys()

    xtick_vals = np.arange(
        0, p2_region[1] - p2_region[0] + 1#, p2_region[1] - p2_region[0]
    )
    xtick_labels = [str(x+p2_region[0]) for x in xtick_vals]

    ytick_vals = np.arange(
        0, p1_region[1] - p1_region[0] + 1#, p1_region[1] - p1_region[0]
    )
    ytick_labels = [str(x+p1_region[0]) for x in ytick_vals]

    num_unique_patches = len(np.unique(contact_map))
    colorscale = generate_cmap(num_unique_patches, scheme="standard")

    if plot_type == "interactive":

        import plotly.graph_objects as go


        fig = go.Figure(
            data=go.Heatmap(
                z=contact_map,
                colorscale=colorscale,
                xgap=0.2,
                ygap=0.2,
            )
        )

        fig.update_layout(
            title="Contact Map",
            yaxis_title=f"Residue number of {chain1}",
            xaxis_title=f"Residue number of {chain2}",
            xaxis=dict(
                tickmode="array",
                tickformat=".0f",
                tickvals=xtick_vals,
                ticktext=xtick_labels,
            ),
            yaxis=dict(
                tickmode="array",
                tickformat=".0f",
                tickvals=ytick_vals,
                ticktext=ytick_labels,
            ),
        )

    elif plot_type == "static":

            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors

            cmap = mcolors.ListedColormap(colorscale)
            fig, ax = plt.subplots()
            # cax = ax.matshow(contact_map, cmap="viridis")

            fig = ax.imshow(
                contact_map,
                cmap=cmap,
                interpolation="nearest",
            )

            xtick_labels = [xtick_labels[0], xtick_labels[-1]]
            ytick_labels = [ytick_labels[0], ytick_labels[-1]]
            xtick_vals = [xtick_vals[0], xtick_vals[-1]]
            ytick_vals = [ytick_vals[0], ytick_vals[-1]]

            ax.set_xticks(xtick_vals)
            ax.set_xticklabels(xtick_labels)
            ax.set_yticks(ytick_vals)
            ax.set_yticklabels(ytick_labels)

            ax.set_xlabel(f"Residue number of {chain2}")
            ax.set_ylabel(f"Residue number of {chain1}")
            ax.set_title("Contact Map")

    return fig


def save_map(
    contact_map: np.array,
    patches: dict,
    # interacting_region: dict,
    chain1: str,
    chain2: str,
    p1_region: tuple,
    p2_region: tuple,
    out_file: str,
    save_plot=False,
    plot_type="static",
):
    """Save the interacting patches and the segmented map to a file.

    Args:
        contact_map (np.array): binary contact map or segmented map
        patches (dict): interacting patches from the map
        interacting_region (dict): interacting region specified by the user
        out_file (str): path to save the output file
        save_plot (bool, optional): save the plot. Defaults to False.
    """

    out_dir = os.path.dirname(out_file)
    file_name = os.path.basename(out_file).split(".")[0]

    # txt_outfile = os.path.join(out_dir, f"{file_name}.txt")

    # print(f"Writing interacting patches to {txt_outfile}")

    # # with open(txt_outfile, "w") as f:
    # #     for patch_id, patch in patches.items():
    # #             f.write(f"Patch {patch_id}\n")
    # #             for chain, res_range in patch.items():
    # #                 f.write(f"{chain}: {res_range}\n")

    csv_outfile = os.path.join(out_dir, f"{file_name}.csv")

    print(f"Writing interacting patches to {csv_outfile}")

    import pandas as pd
    df_rows = []
    for _, patch in patches.items():
        ch1_res_range = patch[chain1]
        ch2_res_range = patch[chain2]
        df_rows.append([ch1_res_range, ch2_res_range])

    df = pd.DataFrame(df_rows, columns=[f"{chain1}", f"{chain2}"])

    save_df(df, csv_outfile)

    if save_plot==True:
        if plot_type == "interactive":
            from utils import plot_map
            fig = plot_map(
                contact_map=contact_map,
                # interacting_region=interacting_region,
                chain1=chain1,
                chain2=chain2,
                p1_region=p1_region,
                p2_region=p2_region,
                plot_type="interactive",
            )
            out_file = os.path.join(out_dir, f"{file_name}.html")
            fig.write_html(
                out_file,
                full_html=False,
            )

        elif plot_type == "static":

            from utils import plot_map
            fig = plot_map(
                contact_map=contact_map,
                # interacting_region=interacting_region,
                chain1=chain1,
                chain2=chain2,
                p1_region=p1_region,
                p2_region=p2_region,
                plot_type="static",
            )

            out_file = os.path.join(out_dir, f"{file_name}.png")
            fig.figure.savefig(out_file)

_select = Select()

def save_pdb(structure: Bio.PDB.Structure, out_file: str, res_select_obj: Bio.PDB.Select = _select):
    """
    Given the ResidueSelect object, save the structure as a PDB file.
    """
    from Bio.PDB import PDBIO
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_file, res_select_obj)

# utils
def get_distance_map(coords1: np.array, coords2: np.array):
    """
    Create an all-v-all distance map.

    Returns a matrix of distances between all pairs of atoms/residues in the two sets of coordinates.
    """

    from scipy.spatial import distance_matrix

    distance_map = distance_matrix(coords1, coords2)

    return distance_map

def get_contact_map(distance_map: np.array, contact_threshold: float):
    """
    Given the distance map, create a binary contact map by thresholding distances.

    Returns a binary matrix, where 1 indicates a contact and 0 indicates no contact.
    """

    contact_map = np.where(
        distance_map <= contact_threshold, 1, 0
    )

    return contact_map

def get_interaction_map(
    coords1: np.array,
    coords2: np.array,
    contact_threshold: float,
    map_type: str
    ):
    """
    Create an interaction map, given the input coordinates.

    Returns a distance map or a contact map, based on the map_type specified.
    """

    distance_map = get_distance_map(coords1, coords2)

    if map_type == "distance":
        return distance_map

    elif map_type == "contact":
        contact_map = get_contact_map(distance_map, contact_threshold)
        return contact_map

    else:
        raise Exception("Invalid map_type specified...")
