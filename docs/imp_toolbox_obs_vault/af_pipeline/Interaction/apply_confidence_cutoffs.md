```python
def apply_confidence_cutoffs(
	self,
	plddt1: dict,
	plddt2: dict,
	pae: np.array
):
	"""
	mask low-confidence interactions.

	Args:
		plddt1 (dict): pLDDT values for chain 1
		plddt2 (dict): pLDDT values for chain 2

	Returns:
		plddt_matrix (np.array): binary matrix for plddt values >= plddt_cutoff
		pae (np.array): binary matrix for pae values <= pae_cutoff
	"""

	chain1, chain2 = next(iter(plddt1)), next(iter(plddt2))
	plddt1, plddt2 = plddt1[chain1], plddt2[chain2]
	plddt1, plddt2 = plddt1.reshape(-1, 1), plddt2.reshape(-1,1)

	ch1_cutoff = ch2_cutoff = self.plddt_cutoff
	if chain1 in self.idr_chains:
		ch1_cutoff = self.idr_plddt_cutoff
	if chain2 in self.idr_chains:
		ch2_cutoff = self.idr_plddt_cutoff

	plddt1 = np.where(plddt1 >= ch1_cutoff, 1, 0)
	plddt2 = np.where(plddt2 >= ch2_cutoff, 1, 0)
	plddt_matrix = plddt1 * plddt2.T

	pae = np.where(pae <= self.pae_cutoff, 1, 0)

	return plddt_matrix, pae
```

### Used in
- [[get_confident_interaction_map]]

### Uses


### Tags
#method 