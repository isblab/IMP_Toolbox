```python
def get_res_range_from_key(res_range: str):
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

    return res_range_list
```