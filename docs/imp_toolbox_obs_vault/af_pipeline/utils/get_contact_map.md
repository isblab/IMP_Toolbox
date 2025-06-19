```python
def get_contact_map(distance_map: np.ndarray, contact_threshold: float):
    """
    Given the distance map, create a binary contact map by thresholding distances.

    Returns a binary matrix, where 1 indicates a contact and 0 indicates no contact.
    """

    contact_map = np.where(
        distance_map <= contact_threshold, 1, 0
    )

    return contact_map
```