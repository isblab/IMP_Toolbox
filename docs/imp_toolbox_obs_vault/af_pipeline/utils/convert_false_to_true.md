```python
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
```