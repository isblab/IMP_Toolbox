import pandas as pd
from collections import Counter
import numpy as np

def symmetrize_matrix(matrix: np.ndarray) -> np.ndarray:
    """Symmetrize a matrix by averaging it with its transpose

    Args:
        matrix (np.ndarray): input matrix

    Returns:
        sym_matrix (np.ndarray): symmetrized matrix

    Example:
    matrix = np.array([[1, 2], [3, 4]]) \n
    symmetrize_matrix(matrix) -> np.array([[1, 2.5], [2.5, 4]])
    """

    assert isinstance(matrix, np.ndarray), "Input must be a numpy array"
    assert matrix.ndim == 2, "Input must be a 2D array"

    sym_matrix = (matrix + matrix.T) / 2

    return sym_matrix


def fill_up_the_blanks(li: list) -> list:
    """Fill up the blanks in a list

    Example:
        fill_up_the_blanks([1, 2, 4, 5]) -> [1, 2, 3, 4, 5]

    Args:
        li (list): list with missing numbers

    Returns:
        new_li (list): list with all the missing numbers filled 
            up between the minimum and maximum values
    """

    min_li_val = min(li)
    max_li_val = max(li)

    new_li = [x for x in range(min_li_val, max_li_val+1)]

    return new_li


def get_key_from_res_range(
    res_range: list,
    as_list=False
) -> str | list:
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


def get_res_range_from_key(res_range: str, return_type: str = "list") -> list | set:
    """Convert a residue range string to a list of residue numbers

    Args:
        res_range (str): residue range string

    Returns:
        res_range_list (list): list of residue numbers

    Example:
    get_res_range_from_key("1-3,5-7") -> [1, 2, 3, 5, 6, 7]
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

def convert_false_to_true(arr: np.ndarray | list, threshold:int=5):
    """ Convert False values in a binary array to True if the patch length is less than or equal to a threshold
        a patch is defined as a sequence of consecutive False values

    Args:
        arr (list): binary array with False values
        threshold (int, optional): _description_. Defaults to 5.

    Returns:
        _type_: _description_
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
):
    """ Get the indices of duplicate elements in a list.

    Args:
        my_li (list): List of elements to check for duplicates.
        return_type (str, optional): Output type, either "list" or "dict".
        keep_which (None | int, optional): If specified, the index of the element to keep in case of duplicates.
            If None, all duplicates are returned. Defaults to 0.

    Returns:
        list: List of indices of duplicate elements if return_type is "list".
        dict: Dictionary with residue IDs as keys and a tuple of first and last indices as values if return_type is "dict".

    Example:
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