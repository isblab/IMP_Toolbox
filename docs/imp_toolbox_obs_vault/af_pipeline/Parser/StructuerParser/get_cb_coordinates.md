```python
def get_cb_coordinates(self):
	"""Get the coordinates of representative atoms for all residues.

	Returns:

		coords_list (list):
			List containing the coordinates for each residue index.
	"""

	coords_list = []

	for residue, _chain_id in self.get_residues():

		coords = self.extract_perresidue_quantity(
			residue=residue, quantity="coords"
		)

		coords_list.append(coords)

	return coords_list
```

### Used in
- [[get_attributes]]

### Uses
- [[get_residues]]
- [[extract_perresidue_quantity]]

### Tags
#method 