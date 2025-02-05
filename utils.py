from Bio.PDB import Select
import Bio
import json
import numpy as np
import pandas as pd
import requests
import os

## API related functions
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


## read and write functions
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

    Args:
        fasta_file (str): Path to fasta file

    Returns:
        all_sequences (dict): dictionary in the format {sequence_header: sequence}
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

# uncomment the line below to test in an interactive shell
##%%

def get_key_from_res_range(res_range: list, as_list=False):
    """Returns a residue range string from a list of residue numbers.

    Args:
        res_range (list): List of residue numbers, e.g., [1, 2, 3, 5, 6, 7]

    Returns:
        str: Residue range string, e.g., "1-3,5-7"

    Example:
    get_key_from_res_range([1, 2, 3, 5, 6, 7]) -> "1-3,5-7" \n
    get_key_from_res_range([1, 2, 3, 5, 6, 7], as_list=True) -> ["1-3", "5-7"]
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
            ranges.append(
                f"{start}-{prev}") if start != prev else ranges.append(str(start)
            )
            start = prev = num

    if start == prev:
        ranges.append(str(start))

    else:
        ranges.append(f"{start}-{prev}")

    if as_list:
        return ranges

    else:
        return ",".join(ranges)

#%%

def get_patches_from_matrix(matrix, chain1, chain2):
    """Get all interacting patches from a binary matrix

    Args:
        matrix (np.array): binary matrix
        chain1 (str): chain 1 identifier
        chain2 (str): chain 2 identifier

    Returns:
        patches (dict): dictionary of interacting patches

    Example:
    matrix = np.array([
        [0, 0, 0, 1],
        [0, 1, 1, 1],
        [0, 0, 1, 1],
        [0, 1, 0, 0],
        [0, 1, 0, 1]
    ])
    get_patches_from_matrix(matrix, "A", "B") ->
    A           B
    {0, 1, 2}   {3}
    {1}         {1, 2, 3}
    {1, 2}      {2, 3}
    {3, 4}      {1}
    {4}         {3}

    """

    assert np.unique(matrix).tolist() == [0, 1]; "Matrix must be binary and non-empty"

    # Required inner functions
    ##%%
    def get_ones(matrix, axis=0):
        """Get the indices of 1s in a binary matrix rowwise or columnwise

        Args:
            matrix (np.array): binary matrix
            axis (int, optional): 0 for rowwise, 1 for columnwise. Defaults to 0.

        Returns:
            one_sets (dict): {k:v} where v is a set of indices of 1s for key k

        Example:
            matrix = np.array(
            [1, 0, 1],
            [0, 1, 0],
            [1, 1, 0]
        ) \n
        get_ones(matrix, axis=0) -> {0: {0, 2}, 1: {1}, 2: {0, 1}} \n
        get_ones(matrix, axis=1) -> {0: {0, 2}, 1: {1, 2}, 2: {0}}
        """
        import numpy as np
        assert np.unique(matrix).tolist() == [0, 1]; "Matrix must be binary"

        one_sets = {}

        if axis == 0: # row_sets
            for i in range(matrix.shape[0]):
                one_sets[i] = set(np.where(matrix[i] == 1)[0])

        elif axis == 1: # col_sets
            for j in range(matrix.shape[1]):
                one_sets[j] = set(np.where(matrix[:, j] == 1)[0])

        return one_sets
    # test
    # import numpy as np
    # matrix = np.array([[1, 0, 1], [0, 1, 0], [1, 1, 0]])
    # print(get_ones(matrix, axis=0))
    # print(get_ones(matrix, axis=1))

   # #%%

    def split_sets(one_sets):
        """Split the sets in one_sets into sub-sets such that
        each subset only contains consecutive indices

        Args:
            one_sets (dict): {k:v} where v is a set of indices of 1s for key k

        Returns:
            new_one_sets (dict): dictionary of lists of lists where each list contains the indices of 1s

        Example:
        one_sets = {0: {0, 1, 2, 3, 5, 6}, 1: {1}, 2: {0, 1}} \n
        split_sets(one_sets) -> {0: [[0, 1, 2, 3], [5, 6]], 1: [[1]], 2: [[0, 1]]}
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
    # test
    # one_sets = {0: {0, 1, 2, 3, 5, 6}, 1: {1}, 2: {0, 1}}
    # print(split_sets(one_sets))
    ##%%

    def extend_sets_by_subsets(one_sets):
        """Add the subsets of the sets in list_of_sets to the one_sets

        Args:
            one_sets (dict): {k:v} where v is a set of indices of 1s for key k

        Returns:
            new_one_sets (dict): {k:v} where v is a list of sets of indices of 1s for key k
                                each set is a subset of the original set and is present in
                                the values of one_sets

        Example:
        one_sets = {
            0: {0, 1, 2, 3, 5, 6},
            1: {1},
            2: {0, 1}
        }
        extend_sets_by_subsets(one_sets) -> {
                0: [{0, 1, 2, 3}, {5, 6}, {1}, {0, 1}],
                1: [{1}],
                2: [{0, 1}, {1}]
            }
        """

        from collections import defaultdict
        split_one_sets = split_sets(one_sets)

        new_one_sets = defaultdict(list)
        list_of_sets = []
        list_of_sets = [
            set(x)
            for xs in split_one_sets.values()
            for x in xs
            if set(x) not in list_of_sets
        ]

        for set1 in list_of_sets:
            for idx, one_set in one_sets.items():
                if set1.issubset(one_set):
                    new_one_sets[idx].append(set1) if set1 not in new_one_sets[idx] else None

        return new_one_sets
    # test
    # one_sets = {
    #     0: {0, 1, 2, 3, 5, 6},
    #     1: {1},
    #     2: {0, 1}
    # }
    # print(extend_sets_by_subsets(one_sets))
    ##%%

    def special_dict_to_df(one_sets, columns):
        """Convert a dictionary to a pandas DataFrame

        Args:
            one_sets (dict): dictionary to convert
            columns (list): column names

        Returns:
            df (pd.DataFrame): DataFrame with the dictionary keys as first column
                               and values as second column in columns

        Example:
        one_sets = {
            1: [{1, 2}, {5}],
            2: [{4, 5}, {6}]
            }
        columns = ["A", "B"]
        special_dict_to_df(one_sets, columns) ->
        A     B
        "1"   {1, 2}
        "1"   {5}
        "2"   {4, 5}
        "2"   {6}

        """

        import pandas as pd

        if all([isinstance(val, list) for val in one_sets.values()]):
            df = pd.DataFrame(one_sets, columns=columns)
            df_rows = []
            for k, v in one_sets.items():
                for val in v:
                    df_rows.append([str(k), val])
            df = pd.DataFrame(df_rows, columns=columns)

        else:
            raise ValueError("All values in the dictionary must be lists.")

        return df
    # test
    # one_sets = {
    #     1: [{1, 2}, {5}],
    #     2: [{4, 5}, {6}]
    # }
    # special_dict_to_df(one_sets, ["A", "B"])
    ##%%

    def groupby_df(df, groupby_col, agg_col):
        """Group a DataFrame by a column and aggregate another column

        Args:
            df (pd.DataFrame): DataFrame with groupby_col and agg_col
            groupby_col (str): column to group by (each value is a set)
            agg_col (str): column to aggregate (each value is a string)

        Returns:
            df_group (pd.DataFrame): grouped DataFrame with both columns as a set

        Example:
        df = pd.DataFrame({
            "A": ["1", "1", "1", "2", "3", "4"],
            "B": [{1}, {1,2}, {5}, {4,5}, {1,2}, {1,2}]
        })
        groupby_df(df, "B", "A") ->
        B       A
        {1}     {1}
        {1, 2}  {1, 3, 4}
        {4, 5}  {2}
        {5}     {1}
        """

        df_group = df.groupby(df[groupby_col].map(tuple))[agg_col].apply(','.join).reset_index()
        for idx, row in df_group.iterrows():
            one_set = row[agg_col].split(",")
            one_set = [int(x) for x in one_set]
            one_set = sorted(one_set)
            df_group.at[idx, agg_col] = set(one_set)

        df_group[groupby_col] = df_group[groupby_col].apply(lambda x: set(x))

        return df_group
    # test
    # import pandas as pd
    # df = pd.DataFrame({
    #     "A": ["1", "1", "1", "2", "3", "4"],
    #     "B": [{1}, {1,2}, {5}, {4,5}, {1,2}, {1,2}]
    # })
    # groupby_df(df, "B", "A")
    ##%%

    def combine_dfs(df1, df2, chain1, chain2):
        """Combine two DataFrames with columns chain1 and chain2
        into a new DataFrame with interacting residues ranges without duplicates

        Args:
            df1 (pd.DataFrame): DataFrame 1
            df2 (pd.DataFrame): DataFrame 2
            chain1 (str): column name 1
            chain2 (str): column name 2

        Returns:
            new_df (pd.DataFrame): combined DataFrame of interacting residues ranges without duplicates

        Example:
        df1 = pd.DataFrame({
            "A":[{1, 3, 4}, {1}, {1, 2}, {0, 1, 2, 4}],
            "B":[{1}, {1, 2, 3}, {2, 3}, {3}]
        })
        df2 = pd.DataFrame({
            "A": [{0, 1, 2}, {1}, {1, 2}, {3, 4}, {4}],
            "B": [{3}, {1, 2, 3}, {2, 3}, {1}, {1, 3}]
        })

        combine_dfs(df1, df2, "A", "B") -> pd.DataFrame
        A   B
        0-2 3
        1   1-3
        1-2 2-3
        3-4 1
        4   1
        4   3
        1   1
        """

        combined_df = pd.concat([df2, df1], axis=0)
        combined_df.reset_index(drop=True, inplace=True)

        new_df = pd.DataFrame(columns=[chain1, chain2])
        df_rows = []

        for _, row in combined_df.iterrows():
            if isinstance(row[chain1], set) and isinstance(row[chain2], set):
                for res_range1 in get_key_from_res_range(row[chain1], as_list=True):
                    for res_range2 in get_key_from_res_range(row[chain2], as_list=True):
                        df_rows.append([res_range1, res_range2])

        new_df = pd.DataFrame(df_rows, columns=[chain1, chain2])
        new_df.drop_duplicates(inplace=True, keep="first")
        new_df.reset_index(drop=True, inplace=True)

        return new_df
    # # test
    # import pandas as pd
    # from utils import get_key_from_res_range
    # df1 = pd.DataFrame({
    #     "A":[{1, 3, 4}, {1}, {1, 2}, {0, 1, 2, 4}],
    #     "B":[{1}, {1, 2, 3}, {2, 3}, {3}]
    # })
    # df2 = pd.DataFrame({
    #     "A": [{0, 1, 2}, {1}, {1, 2}, {3, 4}, {4}],
    #     "B": [{3}, {1, 2, 3}, {2, 3}, {1}, {1, 3}]
    # })
    # chain1 = "A"
    # chain2 = "B"
    # combine_dfs(df1, df2, chain1, chain2)
    ##%%

    def remove_subset_rows(df, chain1, chain2):
        """Remove rows that are subsets of other rows
        (from chatgpt)

        Args:
            df (pd.DataFrame): DataFrame with columns chain1 and chain2
            chain1 (str): column name 1
            chain2 (str): column name 2

        Returns:
            filtered_df (pd.DataFrame): DataFrame with subset rows removed

        Example:
        import pandas as pd
        df = pd.DataFrame({
            "A": [{0, 1, 2}, {1}, {1, 2}, {3, 4}, {4}, {4}, {1}],
            "B": [{3}, {1, 2, 3}, {2, 3}, {1}, {1}, {3}, {1}]
        })
        remove_subset_rows(df, "A", "B")
        A           B
        {0, 1, 2}   {3}
        {1}         {1, 2, 3}
        {1, 2}      {2, 3}
        {3, 4}      {1}
        {4}         {3}
        """

        def is_subset(row, other_row):
            return row[chain1].issubset(other_row[chain1]) and row[chain2].issubset(other_row[chain2])

        rows_to_keep = []
        for i, row in df.iterrows():
            if not any(is_subset(row, df.iloc[j]) for j in range(len(df)) if i != j):
                rows_to_keep.append(i)

        filtered_df = df.loc[rows_to_keep].reset_index(drop=True)

        return filtered_df
    # test
    # import pandas as pd
    # df = pd.DataFrame({
    #     "A": [{0, 1, 2}, {1}, {1, 2}, {3, 4}, {4}, {4}, {1}],
    #     "B": [{3}, {1, 2, 3}, {2, 3}, {1}, {1}, {3}, {1}]
    # })
    # remove_subset_rows(df, "A", "B")
    ##%%

    # actual function
    row_sets = get_ones(matrix, axis=0)
    col_sets = get_ones(matrix, axis=1)

    split_row_sets = extend_sets_by_subsets(row_sets)
    split_col_sets = extend_sets_by_subsets(col_sets)

    df_row = special_dict_to_df(split_row_sets, [chain1, chain2])
    df_col = special_dict_to_df(split_col_sets, [chain2, chain1])

    df_row = groupby_df(df_row, chain2, chain1)
    df_col = groupby_df(df_col, chain1, chain2)

    combined_df = combine_dfs(df_row, df_col, chain1, chain2)

    def res_range_to_set(range_str):
        """ Convert a residue range string to a set of residue numbers

        Args:
            range_str (str): residue range string

        Returns:
            range_set (set): set of residue numbers

        Example:
        res_range_to_set("1-3") -> {1, 2, 3}
        res_range_to_set("1") -> {1}
        """

        range_set = range_str.split("-")

        if len(range_set) == 2:
            range_set_start, range_set_end = map(int, range_set)

        else:
            range_set_start = range_set_end = list(map(int, range_set))[0]

        range_set = set(list(range(range_set_start, range_set_end+1)))

        return range_set

    for col in [chain1, chain2]:
        combined_df[col] = combined_df[col].apply(res_range_to_set)

    combined_df = remove_subset_rows(combined_df, chain1, chain2)

    return combined_df

# test
# binary_matrix = np.array([
#     [0, 0, 1, 1, 0, 0, 0, 1],
#     [0, 1, 1, 1, 0, 1, 1, 1],
#     [0, 1, 1, 1, 1, 1, 1, 1],
#     [0, 0, 0, 0, 0, 1, 1, 0],
#     [0, 0, 1, 1, 1, 0, 0, 0],
#     [0, 0, 1, 1, 0, 0, 1, 0]
# ])

# chain1 = "A"
# chain2 = "B"

# get_patches_from_matrix(binary_matrix, chain1, chain2)

#%%

def generate_cmap(n, scheme="soft-warm"):
    """ Generate a list of n colors
    (modfied from chatgpt)

    Args:
        n (int): number of colors
        scheme (str, optional): Defaults to "soft-warm".

    Returns:
        colors (list): list of n colors
    """

    import random
    import time

    colors = set()
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
            raise ValueError(
                "Invalid scheme. Choose from 'non-bright', 'earth-tone', 'cool-tone', or 'soft-warm'."
            )

        color = "#{:02x}{:02x}{:02x}".format(r, g, b)
        colors.add(color)

        if time.time() - start > 10:
            break

    return list(colors)


def plot_map(contact_map: np.array, chain1: str, chain2: str, p1_region: tuple, p2_region: tuple, plot_type: str):
    """Plot the contact map

    Args:
        contact_map (np.array): binary contact map or segmented map with labels
    """

    xtick_vals = np.arange(
        0, p2_region[1] - p2_region[0] + 1
    )
    xtick_labels = [str(x+p2_region[0]) for x in xtick_vals]

    ytick_vals = np.arange(
        0, p1_region[1] - p1_region[0] + 1
    )
    ytick_labels = [str(x+p1_region[0]) for x in ytick_vals]

    num_unique_patches = len(np.unique(contact_map))

    colorscale = generate_cmap(
        n=num_unique_patches,
        scheme="standard"
    )

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
    chain1: str,
    chain2: str,
    p1_region: tuple,
    p2_region: tuple,
    out_file: str,
    save_plot=False,
    plot_type="static",
):
    """Save the interacting patches and the contact map to a file.

    Args:
        contact_map (np.array): binary contact map or contact map
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
    df.to_csv(csv_outfile, index=False, header=True, sep=",")

    if save_plot==True:
        if plot_type == "interactive" or plot_type == "both":
            from utils import plot_map
            fig = plot_map(
                contact_map=contact_map,
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

            if plot_type == "both":
                plot_type = "static"

        if plot_type == "static":

            from utils import plot_map
            fig = plot_map(
                contact_map=contact_map,
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
