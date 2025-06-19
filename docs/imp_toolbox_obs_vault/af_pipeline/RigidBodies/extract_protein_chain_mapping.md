```python
def extract_protein_chain_mapping(self, protein_chain_mapping: dict):

	protein_chain_map = {}

	if protein_chain_mapping is None:
		return protein_chain_map

	for p_c_maps in protein_chain_mapping:
		protein_name, chain_ids = p_c_maps.split(":")
		chain_ids = chain_ids.split(",")
		for chain_id in chain_ids:
			if chain_id not in protein_chain_map:
				protein_chain_map[chain_id] = protein_name

	return protein_chain_map
```

### Used in
- [[save_rigid_bodies]]

### Uses


### Tags
#method 