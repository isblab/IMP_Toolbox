```python
def get_glycans(self)-> List[Dict[str, Any]]:
	"""Get the glycans of the protein chains
	- If glycans are provided, use those else return an empty list
		- For proteinChain, get the glycans from the entity_info dictionary
	"""

	glycans = []

	if self.entity_type == "proteinChain" and "glycans" in self.entity_info:
		glycans = self.entity_info["glycans"]
		glycans = [
			{
				"residues": glycan[0],
				"position": glycan[1] - self.start + 1,
			}
			for glycan in glycans
		]

	return glycans
```

### Used in 
- [[fill_up_entity]]

### Uses


### Tags
#method