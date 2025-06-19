```python
def get_per_chain_interface_residues(self):
	"""Get interface residues for each chain.

	Returns:
		per_chain_interface_residues (defaultdict): A dictionary where keys are chain IDs and values are lists of residue indices.
		Each list contains the indices of residues that are part of any of the interfaces that the chain is involved in.
	"""

	per_chain_interface_residues = defaultdict(list)

	for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

		chain1, chain2 = chain_pair

		for res1_idx, res2_idx in interacting_res_pairs:

			per_chain_interface_residues[chain1].append(
				res1_idx
			) if res1_idx not in per_chain_interface_residues[chain1] else None

			per_chain_interface_residues[chain2].append(
				res2_idx
			) if res2_idx not in per_chain_interface_residues[chain2] else None

	return per_chain_interface_residues
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 