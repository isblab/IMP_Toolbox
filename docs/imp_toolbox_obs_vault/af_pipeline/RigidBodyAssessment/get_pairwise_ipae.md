```python
def get_pairwise_ipae(self, pae):
	""" Get pairwise iPAE values for each chain pair.

	Args:
		pae (np.ndarray): A 2D numpy array representing the predicted aligned error (PAE) matrix.

	Returns:
		pairwise_ipae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing iPAE values for residue pairs.
	"""

	pairwise_ipae = defaultdict(dict)

	for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

		pairwise_ipae[chain_pair] = {
			"ij" : {
				(res1_idx, res2_idx): pae[res1_idx, res2_idx]
				for res1_idx, res2_idx in interacting_res_pairs
			},
			"ji" : {
				(res2_idx, res1_idx): pae[res2_idx, res1_idx]
				for res1_idx, res2_idx in interacting_res_pairs
			}
		}

	return pairwise_ipae
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 