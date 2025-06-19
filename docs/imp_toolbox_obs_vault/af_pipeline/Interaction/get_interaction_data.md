```python
def get_interaction_data(self, region_of_interest: Dict):
	"""
	Get the interaction amp, pLDDT, and PAE for the region of interest.

	Args:
		region_of_interest (Dict): Dictionary containing the chain IDs and the residue numbers for the region of interest.

	Returns:
		interaction_map (np.array): binary contact map or distance map
		plddt1 (np.array): plddt values for chain 1
		plddt2 (np.array): plddt values for chain 2
		pae (np.array): PAE matrix for the region of interest
	"""

	chain1, chain2 = list(region_of_interest.keys())
	p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

	# chain1 start and end indices.
	start_idx1 = self.num_to_idx[chain1][p1_region[0]]
	end_idx1 = self.num_to_idx[chain1][p1_region[1]]

	# chain2 start and end indices.
	start_idx2 = self.num_to_idx[chain2][p2_region[0]]
	end_idx2 = self.num_to_idx[chain2][p2_region[1]]

	avg_pae = self.avg_pae[start_idx1:end_idx1+1, start_idx2:end_idx2+1]

	coords1 = np.array(self.coords_list[start_idx1:end_idx1+1])
	coords2 = np.array(self.coords_list[start_idx2:end_idx2+1])

	plddt1 = {
		chain1: np.array(self.plddt_list[start_idx1:end_idx1+1])
	}
	plddt2 = {
		chain2: np.array(self.plddt_list[start_idx2:end_idx2+1])
	}

	# Create a contact map or distance map as specified.
	interaction_map = get_interaction_map(
		coords1=coords1.reshape(-1, 3),
		coords2=coords2.reshape(-1, 3),
		contact_threshold=self.contact_threshold,
		map_type=self.interaction_map_type
	)

	return interaction_map, plddt1, plddt2, avg_pae
```

### Used in
- [[get_confident_interaction_map]]

### Uses
- [[get_interaction_map]]

### Tags
#method 