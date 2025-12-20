from collections import defaultdict
import os
from Bio.PDB import Select
import Bio
import Bio.PDB
import Bio.PDB.Structure
import warnings
import pandas as pd
from ruamel.yaml import YAML
import numpy as np
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import (
    get_key_from_res_range,
    get_res_range_from_key,
)
from IMP_Toolbox.utils_imp_toolbox.file_helpers import read_fasta
import importlib.util
import sys
import getpass

_user = getpass.getuser()

# Define the full path to your module file
psa_path = f"/home/{_user}/anaconda3/lib/python3.12/site-packages/psa.py"
module_name = 'pairwise_seq_aln'  # This name will be used to refer to the module

# Create a module specification
spec = importlib.util.spec_from_file_location(
    module_name, psa_path, submodule_search_locations=[psa_path]
)

# Create a new module object based on the spec
pairwise_seq_aln = importlib.util.module_from_spec(spec)

# Add the module to sys.modules to make it available for future imports
sys.modules[module_name] = pairwise_seq_aln

# Execute the module's code in its own namespace
spec.loader.exec_module(pairwise_seq_aln)

def update_config(
    input_file: str,
    updates: dict = None,
    mode: str = "replace",
):
    """Update config file with a new field or update an existing field

    Args:
        input_file (str): Path to input config file
        updates (dict, optional): Fields to update in the config file.
            Defaults to None.
        mode (str, optional): Mode to update the config file.
            Defaults to "replace". ("append" or "replace")
    """

    yaml = YAML()

    update_fields = list(updates.keys()) if updates else []

    if len(update_fields) == 0:

        print("No fields to update in config")
        return None

    yaml.preserve_quotes = True

    with open(input_file, "r") as f:
        config_yaml = yaml.load(f)

    existing_fields = list(config_yaml.keys())

    for field in update_fields:

        add_field = False

        if field in existing_fields:
            if mode == "replace":
                config_yaml[field] = updates[field]
            elif mode == "append":
                # need to change this, not working as expected
                config_yaml[field].update(updates[field])
            else:
                raise ValueError("Invalid mode. Use 'replace' or 'append")

        else:
            print(f"{field} not found in config")
            print("Adding field to config")
            add_field = True

        if add_field:
            config_yaml[field] = updates[field]
            add_field = False

    with open(input_file, "w") as f:
        yaml.dump(config_yaml, f)

    print(f"Config file updated with {update_fields}")


def save_structure_obj(
    structure: Bio.PDB.Structure,
    out_file: str,
    res_select_obj: Bio.PDB.Select = Select(),
    save_type: str = "cif",
    preserve_header_footer=False,
    **kwargs,
):
    """
    Given the ResidueSelect object, save the structure as a PDB file.
    """

    if save_type == "pdb":
        from Bio.PDB import PDBIO

        io = PDBIO()
        io.set_structure(structure)
        io.save(out_file, res_select_obj)

        if preserve_header_footer:
            warnings.warn(
                """
                PDB files do not support headers and footers.
                Saving without headers and footers.
                """
            )

    elif save_type == "cif":

        from Bio.PDB.mmcifio import MMCIFIO

        io = MMCIFIO()
        io.set_structure(structure)
        io.save(out_file, res_select_obj)

        if preserve_header_footer:
            # save structure with headers and footers
            extra_header = structure.header_footer
            if "header" in extra_header and "footer" in extra_header:
                with open(out_file, "r") as f:
                    lines = f.readlines()

                lines.insert(0, "\n")
                lines.append("\n")

                for header_line in extra_header["header"][::-1]:
                    lines.insert(0, f"{header_line}")

                for footer_line in extra_header["footer"][::-1]:
                    lines.append(f"{footer_line}")

                with open(out_file, "w") as f:
                    for line in lines:
                        f.write(line)

            else:
                warnings.warn(
                    """
                    No header or footer information found in the structure.
                    Saving without headers and footers.
                    """
                )

        if "af_offset" in kwargs:
            print(
                f"Adding af_offset to the end of the file: \
                {kwargs['af_offset']}"
            )
            # add af_offset to the end of the file
            af_offset = kwargs["af_offset"]
            af_offset_lines = "".join(
                [
                    f"{key} {" ".join(map(str, val))}\n"
                    for key, val in af_offset.items()
                ]
            )
            with open(out_file, "a") as f:
                f.write(
                    f"\n# \nloop_ \n_af_offset.chain_id \n_af_offset.start \
                    \n_af_offset.end \n{af_offset_lines}\n#"
                )

        if "uniprot_ids" in kwargs:
            # add uniprot_ids to the end of the file
            uniprot_ids = kwargs["uniprot_ids"]
            uniprot_ids_lines = "\n".join(uniprot_ids)
            with open(out_file, "a") as f:
                f.write(f"\n# \n_uniprot_ids {uniprot_ids_lines}\n#")


class MatrixPatches:
    """Class to get interacting patches from a binary matrix"""

    def __init__(
        self,
        matrix: np.ndarray,
        row_obj: str = "row_obj",
        col_obj: str = "col_obj",
    ):
        self.matrix = matrix
        self.row_obj = row_obj
        self.col_obj = col_obj
        self.santiy_check_row_col_obj()

    def santiy_check_row_col_obj(self):

        assert self.row_obj != self.col_obj, (
            """`pd.grouby()` will fail if row_obj and col_obj are the same
            Add suffixes to row_obj and col_obj to make them different.
            """
        )

    def get_patches_from_matrix(self):
        """Get all interacting patches from a binary matrix

        Args:
            matrix (np.ndarray): binary matrix
            row_obj (str): chain 1 identifier
            col_obj (str): chain 2 identifier

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

        assert np.unique(self.matrix).tolist() == [0, 1]; (
        f"Matrix must be binary and non-empty, got {np.unique(self.matrix).tolist()}"
        )

        row_sets = self.get_one_sets_from_matrix(self.matrix, axis=0)
        col_sets = self.get_one_sets_from_matrix(self.matrix, axis=1)

        split_row_sets = self.extend_one_sets_by_subsets(row_sets)
        split_col_sets = self.extend_one_sets_by_subsets(col_sets)

        df_row = self.one_sets_to_df(
            split_row_sets, [self.row_obj, self.col_obj]
        )
        df_col = self.one_sets_to_df(
            split_col_sets, [self.col_obj, self.row_obj]
        )

        df_row = self.aggregate_df_rows(df_row, self.col_obj, self.row_obj)
        df_col = self.aggregate_df_rows(df_col, self.row_obj, self.col_obj)

        combined_df = self.combine_dfs(
            df_row, df_col, self.row_obj, self.col_obj
        )

        for col in [self.row_obj, self.col_obj]:
            combined_df[col] = combined_df[col].apply(
                get_res_range_from_key, return_type="set"
            )

        combined_df = self.remove_subset_rows(
            combined_df, self.row_obj, self.col_obj
        )

        return combined_df

    @staticmethod
    def get_one_sets_from_matrix(matrix: np.ndarray, axis: int = 0):
        """Get the indices of 1s in a binary matrix rowwise or columnwise

        Args:
            matrix (np.ndarray): binary matrix
            axis (int, optional): 0 for rowwise, 1 for columnwise.
                Defaults to 0.

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

        assert np.unique(matrix).tolist() == [0, 1]
        "Matrix must be binary"

        one_sets = {}

        if axis == 0:  # row_sets
            for i in range(matrix.shape[0]):
                one_sets[i] = set(np.where(matrix[i] == 1)[0])

        elif axis == 1:  # col_sets
            for j in range(matrix.shape[1]):
                one_sets[j] = set(np.where(matrix[:, j] == 1)[0])

        return one_sets

    @classmethod
    def extend_one_sets_by_subsets(cls, one_sets: dict) -> dict:
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

        split_sets = MatrixPatches.split_one_sets(one_sets)

        new_one_sets = defaultdict(list)
        list_of_sets = []  # unique sets from split_sets
        list_of_sets = [
            set(x)
            for xs in split_sets.values()
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

    @classmethod
    def split_one_sets(cls, one_sets: dict) -> dict:
        """Split the sets in `one_sets` into sub-sets such that
        each subset only contains consecutive indices

        Args:
            one_sets (dict): {k:v} where v is a set of indices of 1s for key k

        Returns:
            new_one_sets (dict): dictionary of lists of lists where each list
                contains the indices of 1s

        Example:
        one_sets = {0: {0, 1, 2, 3, 5, 6}, 1: {1}, 2: {0, 1}} \n
        split_sets(one_sets) -> {0: [[0, 1, 2, 3], [5, 6]], 1: [[1]], 2: [[0, 1]]}
        """

        new_one_sets = {}

        for i, one_set in one_sets.items():

            if not isinstance(one_set, set):
                raise TypeError("one_set must be a set")

            sub_sets = MatrixPatches.split_one_set(one_set)
            new_one_sets[i] = sub_sets

        return new_one_sets

    @staticmethod
    def split_one_set(one_set: set | list) -> list:
        """Split a set of indices into sub-sets such that
        each subset only contains consecutive indices

        Args:
            one_set (set | list): set of indices of 1s

        Returns:
            sub_sets (list): list of lists where each list
                contains the indices of 1s

        Example:
        one_set = {0, 1, 2, 3, 5, 6} \n
        split_one_set(one_set) -> [[0, 1, 2, 3], [5, 6]]
        """

        assert isinstance(
            one_set, set | list
        ), "one_set must be a set or a list"

        sub_sets = []

        if isinstance(one_set, list):
            # need to remove duplicates if any
            one_set = set(one_set)

        one_set = sorted(list(one_set))

        for idx, one_pos in enumerate(one_set):

            curr_pos = one_pos
            prev_pos = one_set[idx - 1] if idx > 0 else None

            if idx == 0:
                # If it's the first position, create a new sub-set
                sub_sets.append([curr_pos])

            elif curr_pos - prev_pos == 1:
                # If the current position is consecutive to the previous one
                # add it to the last sub-set
                sub_sets[-1].append(one_pos)

            else:
                # If the current position is not consecutive to the previous one
                # create a new sub-set
                sub_sets.append([curr_pos])

        return sub_sets

    @staticmethod
    def one_sets_to_df(one_sets: dict, columns: list) -> pd.DataFrame:
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

        if all([isinstance(val, list) for val in one_sets.values()]):

            df_rows = []

            for k, v in one_sets.items():
                for val in v:
                    df_rows.append([str(k), val])

            df = pd.DataFrame(df_rows, columns=columns)

        else:
            raise ValueError("All values in the dictionary must be lists.")

        return df

    @staticmethod
    def aggregate_df_rows(
        df: pd.DataFrame, groupby_col: str, agg_col: str
    ) -> pd.DataFrame:
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

        df_group = (
            df.groupby(df[groupby_col].map(tuple))[agg_col]
            .apply(",".join)
            .reset_index()
        )
        for idx, row in df_group.iterrows():
            one_set = row[agg_col].split(",")
            one_set = [int(x) for x in one_set]
            one_set = sorted(one_set)
            df_group.at[idx, agg_col] = set(one_set)

        df_group[groupby_col] = df_group[groupby_col].apply(lambda x: set(x))

        return df_group

    @staticmethod
    def combine_dfs(
        df1: pd.DataFrame, df2: pd.DataFrame, colname_1: str, colname_2: str
    ) -> pd.DataFrame:
        """Combine two DataFrames with columns colname_1 and colname_2
        into a new DataFrame with interacting residues ranges without duplicates

        Args:
            df1 (pd.DataFrame): DataFrame 1
            df2 (pd.DataFrame): DataFrame 2
            colname_1 (str): column name 1
            colname_2 (str): column name 2

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

        new_df = pd.DataFrame(columns=[colname_1, colname_2])
        df_rows = []

        for _, row in combined_df.iterrows():
            if isinstance(row[colname_1], set) and isinstance(
                row[colname_2], set
            ):
                for res_range1 in get_key_from_res_range(
                    row[colname_1], as_list=True
                ):
                    for res_range2 in get_key_from_res_range(
                        row[colname_2], as_list=True
                    ):
                        df_rows.append([res_range1, res_range2])

        new_df = pd.DataFrame(df_rows, columns=[colname_1, colname_2])
        new_df.drop_duplicates(inplace=True, keep="first")
        new_df.reset_index(drop=True, inplace=True)

        return new_df

    def remove_subset_rows(
        self, df: pd.DataFrame, colname_1: str, colname_2: str
    ):
        """Remove rows that are subsets of other rows
        (from chatgpt)

        Args:
            df (pd.DataFrame): DataFrame with columns colname_1 and colname_2
            colname_1 (str): column name 1
            colname_2 (str): column name 2

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

        rows_to_keep = []

        for i, row in df.iterrows():

            if not any(
                self.is_subset(row, df.iloc[j], colname_1, colname_2)
                for j in range(len(df))
                if i != j
            ):
                rows_to_keep.append(i)

        filtered_df = df.loc[rows_to_keep].reset_index(drop=True)

        return filtered_df

    @staticmethod
    def is_subset(
        row: pd.Series,
        other_row: pd.Series,
        colname_1: str,
        colname_2: str,
    ):
        """Check if row is a subset of other_row for two specified columns

        Args:
            row (pd.Series): _description_
            other_row (pd.Series): _description_
            colname_1 (str): _description_
            colname_2 (str): _description_

        Returns:
            bool: True if row is a subset of other_row, False otherwise
        """

        return row[colname_1].issubset(other_row[colname_1]) and row[
            colname_2
        ].issubset(other_row[colname_2])


def pairwise_alignment_map(
    pairwise_alignment_file:str,
    include_aligned_seq: bool=False,
) -> dict:
    """ Given a pairwise alignment, return one-one map of aligned residues.

    Args:
        pairwise_alignment_file (str):
            Path to the pairwise alignment file

        include_aligned_seq (bool, optional):
            Whether to return the aligned sequences

    Returns:
        dict: mapping of codon numbers to residue numbers
    """

    pairwise_alignment = read_fasta(pairwise_alignment_file)

    if len(pairwise_alignment) != 2:
        raise ValueError(
            f"""
            Expected 2 sequences in the pairwise alignment file,
            found {len(pairwise_alignment)}.
            """
        )

    pairwise_alignment_dict = {}
    q_count, s_count = 0, 0
    qseq, sseq = list(pairwise_alignment.values())

    increment_dict = {
        (True, True): lambda q, s: (q + 1, s + 1),
        (True, False): lambda q, s: (q + 1, s),
        (False, True): lambda q, s: (q, s + 1),
        (False, False): lambda q, s: (q, s),
    }

    init_case = {
        True: lambda d: d + 1,
        False: lambda d: d,
    }

    for idx in range(len(qseq)):

        inc_case = (sseq[idx] != "-", qseq[idx] != "-")
        q_count, s_count = increment_dict[inc_case](q_count, s_count)

        if sseq[idx] != "-" and qseq[idx] != "-":
            pairwise_alignment_dict[init_case[q_count==0](q_count)] = (
                init_case[s_count==0](s_count)
            )

    if include_aligned_seq:
        return pairwise_alignment_dict, qseq, sseq
    else:
        return pairwise_alignment_dict

def handle_pairwise_alignment(
    p_name: str,
    sseq: str,
    qseq: str,
    pairwise_alignment_file: str,
    alignment_program: str = "stretcher",
    ignore_warnings: bool = False,
    include_identical: bool = False,
    include_aligned_seq: bool = False,
):
    """ Handle pairwise alignment between two sequences and return the mapping.

    Args:
        p_name (str): Protein name
        sseq (str): Subject sequence
        qseq (str): Query sequence
        pairwise_alignments_dir (str): Directory to save pairwise alignments

    Returns:
        dict: Mapping of codon to residue number
    """

    # no need to align if sequences are identical
    if sseq == qseq and not include_identical:
        psa_map = {}

    else:
        os.makedirs(os.path.dirname(pairwise_alignment_file),exist_ok=True)

        if not os.path.exists(pairwise_alignment_file):

            pairwise_alignment = pairwise_seq_aln.align(
                program=alignment_program,
                moltype="prot",
                qseq=qseq,
                sseq=sseq,
            )

            with open(pairwise_alignment_file, "w") as f:
                f.write(pairwise_alignment.fasta())
            print(f"Pairwise alignment saved to {pairwise_alignment_file}")

        psa_output = pairwise_alignment_map(
            pairwise_alignment_file,
            include_aligned_seq=include_aligned_seq,
        )

        aligned_qseq = None
        aligned_sseq = None

        if include_aligned_seq:
            psa_map, aligned_qseq, aligned_sseq = psa_output
        else:
            psa_map = psa_output

        # warn if sequences not identical
        if not ignore_warnings:
            warnings.warn(
                f"""
                Sequences not identical for {p_name}. Cross-check the pairwise alignment.\n
                Pairwise alignment for {p_name}:\n{pairwise_alignment.raw}\n
                """
            )
            # print(f"Codon to residue map: {psa_map}")

    if include_aligned_seq:
        return psa_map, aligned_qseq, aligned_sseq
    else:
        return psa_map