import pandas as pd
from collections import Counter
import numpy as np

def symmetrize_matrix(matrix: np.ndarray) -> np.ndarray:
    """ Symmetrize a matrix by averaging it with its transpose.

    ## Arguments:

    - **matrix (np.ndarray)**:<br />
        Input matrix to be symmetrized.

    ## Returns:

    - **np.ndarray**:<br />
        Symmetrized matrix.

    ## Examples:

    >>> matrix = np.array([[1, 2], [3, 4]])
    >>> symmetrize_matrix(matrix)
    array([[1., 2.5],
           [2.5, 4.]])
    """

    assert isinstance(matrix, np.ndarray), "Input must be a numpy array"
    assert matrix.ndim == 2, "Input must be a 2D array"

    sym_matrix = (matrix + matrix.T) / 2

    return sym_matrix

def fill_up_the_blanks(li: list) -> list:
    """ Fill up the blanks in a list by adding missing integers between the minimum and maximum values.

    ## Arguments:

    - **li (list)**:<br />
        List of integers with missing values.

    ## Returns:

    - **list**:<br />
        List of integers with missing values filled in.

    ## Examples:

    >>> fill_up_the_blanks([1, 2, 4, 5])
    [1, 2, 3, 4, 5]
    """

    min_li_val = min(li)
    max_li_val = max(li)

    new_li = [x for x in range(min_li_val, max_li_val+1)]

    return new_li

def get_key_from_res_range(
    res_range: list,
    as_list=False
) -> str | list:
    """ Return a residue range string from a list of residue numbers.

    ## Arguments:

    - **res_range (list)**:<br />
        List of residue numbers.

    - **as_list (bool, optional):**:<br />
        Whether to return a list of strings or single string.

    ## Returns:

    - **str | list**:<br />
        Residue range string or list of residue range strings.

    ## Examples:

    >>> get_key_from_res_range([1, 2, 3, 5, 6, 7])
    '1-3,5-7'
    >>> get_key_from_res_range([1, 2, 3, 5, 6, 7], as_list=True)
    ['1-3', '5-7']
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

def get_res_range_from_key(
    res_range: str,
    return_type: str = "list"
) -> list | set:
    """ Convert a residue range string to alist of residue numbers.

    ## Arguments:

    - **res_range (str)**:<br />
        Residue range string, e.g. "1-3,5-7".

    - **return_type (str, optional):**:<br />
        Whether to return a list or a set of residue numbers. Defaults to "list".

    ## Returns:

    - **list | set**:<br />
        List or set of residue numbers.

    ## Examples:

    >>> get_res_range_from_key("1-3,5-7")
    [1, 2, 3, 5, 6, 7]
    >>> get_res_range_from_key("1-3,5-7", return_type="set")
    {1, 2, 3, 5, 6, 7}
    """

    res_range_list = []

    for res_range in res_range.split(","):

        if "-" in res_range:
            start, end = map(int, res_range.split("-"))
            res_range_list.extend(list(range(start, end+1)))

        else:
            res_range_list.append(int(res_range))

    if return_type == "set":
        return set(res_range_list)

    return res_range_list

def convert_false_to_true(arr: np.ndarray | list, threshold:int=5) -> np.ndarray:
    """ Convert False values to True in a binary True-False array if
    consecutive False values (a patch) are less than or equal to a specified threshold.

    ## Arguments:

    - **arr (np.ndarray | list)**:<br />
        Binary array (list or numpy array) containing True and False values.

    - **threshold (int, optional):**:<br />
        Maximum length of consecutive False values (patch) to be converted to True. Defaults to 5.

    ## Returns:

    - **np.ndarray**:<br />
        Binary array with False values converted to True if they are part of a patch of length <= threshold.

    ## Examples:

    >>> arr = [True, False, False, True, False, False, False, True]
    >>> convert_false_to_true(arr, threshold=2)
    array([ True,  True,  True,  True, False, False, False,  True])
    """

    if isinstance(arr, list):
        arr = np.array(arr)

    where_false = list(np.argwhere(arr == False).flatten())

    false_patches = get_key_from_res_range(where_false, as_list=True)

    for patch in false_patches:

        patch = patch.split("-")
        if len(patch) == 1:
            patch.append(patch[0])

        start = int(patch[0])
        end = int(patch[1])

        if end - start + 1 <= threshold:
            arr[start:end + 1] = True

    return arr

def get_duplicate_indices(
    my_li: list,
    return_type: str = "list",
    keep_which: None | int = 0,
) -> list | dict:
    """ Get the indices of duplicate elements in a list.

    ## Arguments:

    - **my_li (list)**:<br />
        List of elements to check for duplicates.

    - **return_type (str, optional):**:<br />
        Output type, either "list" or "dict". Defaults to "list".

    - **keep_which (None | int, optional):**:<br />
        If specified, the index of the element to keep in case of duplicates.
        If None, all duplicates are returned. Defaults to 0.

    ## Returns:

    - **list | dict**:<br />
        List of indices of duplicate elements if return_type is "list".
        Dictionary with residue IDs as keys and a tuple of first and last
        indices as values if return_type is "dict".

    ## Examples:

    >>> get_duplicate_indices(['A', 'B', 'C', 'A', 'D'], return_type='list')
    [0, 3]
    >>> get_duplicate_indices(['A', 'B', 'C', 'A', 'D'], return_type='dict')
    {'A': (0, 3)}
    """

    token_counts = Counter(my_li)
    # Get the repeated residue IDs, these are pTMs or small molecules
    repeated_tokens = [
        token_id
        for token_id, count
        in token_counts.items()
        if count > 1
    ]

    duplicate_indices = []

    for token_id in repeated_tokens:
        indices = [i for i, x in enumerate(my_li) if x == token_id]
        if keep_which is not None:
            indices.pop(keep_which)
        duplicate_indices.extend(indices)

    if return_type == "list":
        return duplicate_indices

    elif return_type == "dict":

        duplicate_indices = {}

        for res in repeated_tokens:
            indices = [i for i, x in enumerate(my_li) if x == res]
            duplicate_indices[res] = (indices[0], indices[-1])

        return duplicate_indices

# def fasta_str_to_dict(fasta_str):
#     """ Convert a FASTA string to a dictionary of sequences.

#     Args:
#         fasta_str (str): FASTA string

#     Returns:
#         dict: Dictionary of sequences.

#     Example:
#     >>> fasta_str = '''>seq1
#     ... ABCD
#     ... >seq2
#     ... ABCD'''

#     >>> fasta_str_to_dict(fasta_str)
#     {'seq1': 'ABCD', 'seq2': 'ABCD'}
#     """

#     fasta_str = fasta_str.splitlines()
#     fasta_dict = {}
#     for i, line in enumerate(fasta_str):
#         if line.startswith(">"):
#             seq_name = line[1:].strip()
#             fasta_dict[seq_name] = ""
#         else:
#             fasta_dict[seq_name] += line.strip()
#     return fasta_dict