```python
@staticmethod
def get_distance_map(
	coords1: np.ndarray,
	coords2: np.ndarray
):
	"""
	Create an all-v-all distance map.

	Returns a matrix of distances between all pairs of atoms/residues in the two sets of coordinates.
	"""

	from scipy.spatial import distance_matrix

	distance_map = distance_matrix(coords1, coords2)

	return distance_map
```
