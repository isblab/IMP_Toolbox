```python
def get_per_chain_iplddt(self, plddt_list):
	""" Get per-chain ipLDDT scores from a list of pLDDT scores.

	Args:
		plddt_list (list): A list of pLDDT scores for all residues in the structure.

	Returns:
		per_chain_iplddt (defaultdict): A dictionary where keys are chain IDs and values are dictionaries mapping residue indices to their pLDDT scores.
	"""

	per_chain_iplddt = defaultdict(dict)

	for chain_id, interface_res_idxs in self.per_chain_interface_residues.items():

		for res_idx in interface_res_idxs:

			per_chain_iplddt[chain_id][res_idx] = plddt_list[res_idx]

	return per_chain_iplddt
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 