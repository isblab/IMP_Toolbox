```python
def fill_up_the_blanks(li: list):
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
```