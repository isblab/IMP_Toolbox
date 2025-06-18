```python
def get_modifications(self)-> List[Dict[str, Any]]:
	"""Get the modifications of the entity

	- If modifications are provided, use those else empty list

		- For proteinChain, get the modifications from the
			entity_info dictionary (ptmType, ptmPosition)

		- For dnaSequence and rnaSequence, get the modifications from the
			entity_info dictionary (modificationType, basePosition)
	"""

	modifications = self.entity_info.get("modifications", [])

	if "modifications" in self.entity_info:

		if self.entity_type == "proteinChain":
			modifications = [
				{
					"ptmType": mod[0],
					"ptmPosition": mod[1] - self.start + 1,
				}
				for mod in modifications
			]

		elif (
			self.entity_type == "dnaSequence"
			or self.entity_type == "rnaSequence"
		):
			modifications = [
				{
					"modificationType": mod[0],
					"basePosition": mod[1] - self.start + 1,
				}
				for mod in modifications
			]

		else:
			raise Exception(
				"Modifications are not supported for this entity type"
			)

	return modifications
```

### Used in
- [[fill_up_entity]]

### Uses


### Tags
#method