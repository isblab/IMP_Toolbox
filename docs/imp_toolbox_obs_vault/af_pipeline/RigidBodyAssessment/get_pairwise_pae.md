```python
def get_pairwise_pae(self, pae):
	""" Get pairwise PAE values for each chain pair.

	Args:
		pae (np.ndarray): A 2D numpy array representing the predicted aligned error (PAE) matrix.

	Returns:
		pairwise_pae (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing PAE values for residue pairs.
	"""

	pairwise_pae = defaultdict(np.ndarray)

	for chain_pair in self.chain_pairs:

		rb_chain_pair_res = self.rb_res_pairs[chain_pair]

		rb_pae_vals_ij = [
			pae[res1_idx, res2_idx] for res1_idx, res2_idx in rb_chain_pair_res
		]

		rb_pae_vals_ji = [
			pae[res2_idx, res1_idx] for res1_idx, res2_idx in rb_chain_pair_res
		]

		if len(rb_pae_vals_ij) > 0:
			pairwise_pae[chain_pair] = {
				"ij": rb_pae_vals_ij,
				"ji": rb_pae_vals_ji,
			}

	return pairwise_pae
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 