```python
def get_num_interface_residues(self):
	"""Get the number of interface residues for each chain pair.

	Returns:
		num_interface_residues (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the number of unique interface residues.
		Each key is a tuple of two chain IDs, and the value is the count of unique residues that interact between those chains.
	"""

	num_interface_residues = defaultdict(int)

	for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

		unique_interface_residues = np.unique(np.array(interacting_res_pairs).flatten())
		num_interface_residues[chain_pair] = len(unique_interface_residues)

	return num_interface_residues
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 