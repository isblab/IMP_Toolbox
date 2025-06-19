```python
def get_patches_from_matrix(matrix, chain1, chain2):
    """Get all interacting patches from a binary matrix

    Args:
        matrix (np.ndarray): binary matrix
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
            matrix (np.ndarray): binary matrix
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
                    (
                        new_one_sets[idx].append(set1)
                        if set1 not in new_one_sets[idx]
                        else None
                    )

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
```
