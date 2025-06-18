```python
def sanity_check_glycans(self):
	"""Sanity check the glycans
	- check if the glycosylation position is valid 
		(should be within the provided sequence)
	- glycans are only supported for proteinChain, raise exception otherwise
	"""

	if self.entity_type == "proteinChain" and len(self.glycans) > 0:

		# check if the glycosylation position is valid
		for glycan in self.glycans:
			glyc_pos = glycan["position"]

			if glyc_pos < 1 or glyc_pos > len(self.real_sequence):
				raise Exception(
					f"Invalid glycan position at {glyc_pos} \
					in {self.entity_name}"
				)

	if self.entity_type != "proteinChain" and len(self.glycans) > 0:
		raise Exception("Glycosylation is not supported for this entity type")
```

### Used in 
- [[Entity]]

### Uses


### Tags
#method