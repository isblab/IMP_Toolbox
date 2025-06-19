```python
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
```