```python
def get_confident_interaction_map(self, region_of_interest: Dict):
	"""
	For the specified regions in the predicted structure, obtain all confident interacting residue pairs.

	Returns:
		confident_interactions (np.array): binary map of confident interacting residues
	"""

	interaction_map, plddt1, plddt2, avg_pae = self.get_interaction_data(
		region_of_interest=region_of_interest
	)

	plddt_matrix, pae_matrix = self.apply_confidence_cutoffs(
		plddt1=plddt1, plddt2=plddt2, pae=avg_pae
	)

	confident_interactions = interaction_map * plddt_matrix * pae_matrix

	return confident_interactions
```

### Used in
- [[save_ppair_interaction]]

### Uses
- [[get_interaction_data]]
- [[apply_confidence_cutoffs]]

### Tags
#method 