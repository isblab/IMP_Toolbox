```python
def get_duplicate_indices(
    my_li: list,
    return_type: str = "list",
    keep_which: None | int = 0,
):
    from collections import Counter

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
```