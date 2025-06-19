```python
def get_per_chain_plddt(self, plddt_list):
	""" Get per-chain pLDDT scores from a list of pLDDT scores.

	Args:
		plddt_list (list): A list of pLDDT scores for all residues in the structure.

	Returns:
		per_chain_plddt (defaultdict): A dictionary where keys are chain IDs and values are numpy arrays of pLDDT scores for residues in that chain.
	"""

	per_chain_plddt = defaultdict(np.ndarray)

	for chain_id, res_list in self.rb_dict.items():

		res_idxs = [
			self.num_to_idx[chain_id][res_num] for res_num in res_list
		]

		plddt_scores = np.array(plddt_list)[res_idxs]
		per_chain_plddt[chain_id] = plddt_scores

	return per_chain_plddt
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 