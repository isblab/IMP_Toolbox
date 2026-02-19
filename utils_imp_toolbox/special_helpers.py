from collections import defaultdict
import os
import psa
import warnings
import Bio
import Bio.PDB
import Bio.PDB.Structure
import numpy as np
import pandas as pd
from Bio.PDB import Select
from ruamel.yaml import YAML
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import (
    get_key_from_res_range,
    get_res_range_from_key,
)
from IMP_Toolbox.utils_imp_toolbox.file_helpers import read_fasta

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
        df_group[agg_col] = df_group[agg_col].astype(object)
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

        sets_1 = df[colname_1].tolist()
        sets_2 = df[colname_2].tolist()
        idxs = df.index.tolist()

        sets_ = list(zip(sets_1, sets_2))
        common_subset_idxs = [
            i for i in idxs for j in idxs
            if (i != j and (
                sets_[i][0].issubset(sets_[j][0]) and
                sets_[i][1].issubset(sets_[j][1])
            ))
        ]

        rows_to_keep = [i for i in idxs if i not in common_subset_idxs]

        filtered_df = df.loc[rows_to_keep].reset_index(drop=True)

        return filtered_df

def get_closest_mapped_residue(
    psa_map: dict,
    codon_number: int,
    which: str = "lower",
):
    """ Given a pairwise alignment map and a codon number, return the
    closest mapped residue number.

    Args:

        psa_map (dict):
            Pairwise alignment map of codon number to residue number

        codon_number (int):
            Codon number to map
        which (str, optional):
            Whether to return the closest lower or higher residue number.
            Defaults to "lower".
    Returns:
        int:
            Closest mapped residue number
    """

    mapped_residues = sorted(psa_map.keys())

    if which == "lower":
        lower_residues = [
            res for res in mapped_residues if res <= codon_number
        ]
        if lower_residues:
            closest_residue = psa_map[lower_residues[-1]]
        else:
            warnings.warn(
                f"""
                No lower residue found for codon number {codon_number}
                in the pairwise alignment map.
                Returning higher residue instead.
                """
            )
            closest_residue = get_closest_mapped_residue(
                psa_map, codon_number, which="higher"
            )

    elif which == "higher":
        higher_residues = [
            res for res in mapped_residues if res >= codon_number
        ]
        if higher_residues:
            closest_residue = psa_map[higher_residues[0]]
        else:
            warnings.warn(
                f"""
                No higher residue found for codon number {codon_number}
                in the pairwise alignment map.
                Returning lower residue instead.
                """
            )
            closest_residue = get_closest_mapped_residue(
                psa_map, codon_number, which="lower"
            )

    else:
        raise ValueError("which must be 'lower' or 'higher'")

    return closest_residue

def pairwise_alignment_map(
    pairwise_alignment_file: str,
) -> dict:
    """ Given a pairwise alignment, return one-one map of aligned residues.

    Args:
        pairwise_alignment_file (str):
            Path to the pairwise alignment file

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

    for q_res, s_res in zip(qseq, sseq):

        inc_case = (s_res != "-", q_res != "-")
        q_count, s_count = increment_dict[inc_case](q_count, s_count)

        if s_res != "-" and q_res != "-":
            pairwise_alignment_dict[init_case[q_count==0](q_count)] = (
                init_case[s_count==0](s_count)
            )

    return pairwise_alignment_dict

def handle_pairwise_alignment(
    sseq: str,
    qseq: str,
    pairwise_alignment_file: str,
    alignment_program: str = "stretcher",
    xtra_attributes: list = [],
    overwrite: bool = True,
):
    """ Handle pairwise alignment between two sequences and return the mapping.

    Args:
        p_name (str): Protein name
        sseq (str): Subject sequence
        qseq (str): Query sequence
        pairwise_alignments_file (str): Path to the pairwise alignment file
        alignment_program (str, optional): Alignment program to use.
        xtra_attributes (list, optional): Extra attributes to extract from
            the pairwise alignment object. Defaults to [].
        overwrite (bool, optional): Whether to overwrite existing alignment file.

    Returns:
        dict: Mapping of codon to residue number
    """

    os.makedirs(os.path.dirname(pairwise_alignment_file),exist_ok=True)

    if not os.path.exists(pairwise_alignment_file) or overwrite:

        pairwise_alignment = psa.align(
            program=alignment_program,
            moltype="prot",
            qseq=qseq,
            sseq=sseq,
        )

        with open(pairwise_alignment_file, "w") as f:
            f.write(pairwise_alignment.fasta())

        print(f"Pairwise alignment saved to {pairwise_alignment_file}")

    psa_map = pairwise_alignment_map(pairwise_alignment_file)

    xtra = {}
    for attr in xtra_attributes:
        xtra[attr] = getattr(pairwise_alignment, attr, None)

    return psa_map, xtra

def get_mapped_residue(
    psa_map: dict,
    codon_number: int,
    p_name: str = "",
) -> tuple[int | None, str]:
    """ Given a pairwise alignment map and a codon number, return the
    mapped residue number.

    Args:

        psa_map (dict):
            Pairwise alignment map of codon number to residue number

        codon_number (int):
            Codon number to map

        p_name (str, optional):
            Protein name for warning messages

    Returns:

        tuple[int | None, str]:
            Mapped residue number and warning message (if any)
    """

    warn_msg = ""

    if len(psa_map) == 0:
        res_num_mapped = codon_number
    else:
        try:
            res_num_mapped = psa_map[codon_number]
        except KeyError:
            warn_msg += (
                f"""
                Residue number {codon_number} not found in
                pairwise alignment map for protein {p_name}.
                Skipping...
                """
            )
            res_num_mapped = None

    return res_num_mapped, warn_msg

def get_seq_identity(
    qseq: str,
    sseq: str,
    start: int,
    end: int,
    as_percentage: bool=True
):
    """ Calculate sequence identity between two aligned sequences in a given range.

    Args:
        qseq (str): Query sequence (aligned)
        sseq (str): Subject sequence (aligned)
        start (int): Start position (1-based)
        end (int): End position (1-based)

    Returns:
        float: Sequence identity percentage
    """

    seq_count = 0
    start_idx = end_idx = None
    for i, a in enumerate(qseq):
        if a != '-':
            seq_count += 1
        if seq_count == start and start_idx is None:
            start_idx = i
        if seq_count == end:
            end_idx = i
            break

    if start_idx is None or end_idx is None:
        return 0.0

    aligned_qseq = qseq[start_idx:end_idx+1]
    aligned_sseq = sseq[start_idx:end_idx+1]

    aligned_qseq = qseq[start-1:end]
    aligned_sseq = sseq[start-1:end]

    matches = sum(
        1 for a, b in zip(aligned_qseq, aligned_sseq)
        if (a == b and a != '-' and b != '-')
    )
    length = len(aligned_qseq)

    if length == 0:
        return 0.0

    if as_percentage:
        identity = (matches / length) * 100
    else:
        identity = matches

    return identity

def get_gap(
    qseq: str,
    sseq: str,
    start: int,
    end: int,
    as_percentage: bool=True
):
    """ Calculate gap percentage between two aligned sequences in a given range.

    Args:
        qseq (str): Query sequence (aligned)
        sseq (str): Subject sequence (aligned)
        start (int): Start position (1-based)
        end (int): End position (1-based)

    Returns:
        float: Gap percentage
    """

    seq_count = 0
    start_idx = end_idx = None
    for i, a in enumerate(qseq):
        if a != '-':
            seq_count += 1
        if seq_count == start and start_idx is None:
            start_idx = i
        if seq_count == end:
            end_idx = i
            break

    if start_idx is None or end_idx is None:
        return 100.0

    aligned_qseq = qseq[start_idx:end_idx+1]
    aligned_sseq = sseq[start_idx:end_idx+1]

    aligned_qseq = qseq[start-1:end]
    aligned_sseq = sseq[start-1:end]

    gaps = sum(
        1 for a, b in zip(aligned_qseq, aligned_sseq)
        if (a == '-' or b == '-')
    )
    length = len(aligned_qseq)

    if length == 0:
        return 100.0

    if as_percentage:
        gap_percentage = (gaps / length) * 100
    else:
        gap_percentage = gaps

    return gap_percentage

def fasta_str_to_dict(fasta_str):
    """ Convert a FASTA string to a dictionary of sequences.

    Args:
        fasta_str (str): FASTA string

    Returns:
        dict: Dictionary of sequences.

    Example:
    >>> fasta_str = '''>seq1
    ... ABCD
    ... >seq2
    ... ABCD'''

    >>> fasta_str_to_dict(fasta_str)
    {'seq1': 'ABCD', 'seq2': 'ABCD'}
    """

    fasta_str = fasta_str.splitlines()
    fasta_dict = {}
    for i, line in enumerate(fasta_str):
        if line.startswith(">"):
            seq_name = line[1:].strip()
            fasta_dict[seq_name] = ""
        else:
            fasta_dict[seq_name] += line.strip()
    return fasta_dict