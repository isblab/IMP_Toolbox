```python
def sanity_check_modifications(self):
	"""Sanity check the modifications
	- check if the modification type is valid
		(should be in the allowed modifications)
	- check if the modification position is valid
		(should be within the provided sequence)
	- modifications are only supported for
		`proteinChain`, `dnaSequence`, `rnaSequence`; raise exception otherwise
	"""

	if (
		self.entity_type not in ["proteinChain", "dnaSequence", "rnaSequence"]
		and len(self.modifications) > 0
	):
		raise Exception("Modifications are not supported for this entity type")

	# if (
	#     self.entity_type in ["proteinChain", "dnaSequence", "rnaSequence"]
	#     and len(self.modifications) > 0
	# ):

	# check if the modification type is valid
	if self.entity_type == "proteinChain":
		if not all([mod["ptmType"] in PTM for mod in self.modifications]):
			raise Exception("Invalid modification type")

	elif self.entity_type == "dnaSequence":
		if not all(
			[mod["modificationType"] in DNA_MOD for mod in self.modifications]
		):
			raise Exception("Invalid modification type")

	elif self.entity_type == "rnaSequence":
		if not all(
			[mod["modificationType"] in RNA_MOD for mod in self.modifications]
		):
			raise Exception("Invalid modification type")

	# check if the modification position is valid
	for mod in self.modifications:
		mod_pos = (
			mod["ptmPosition"]
			if self.entity_type == "proteinChain"
			else mod["basePosition"]
		)

		if mod_pos < 1 or mod_pos > len(self.real_sequence):
			raise Exception(
				f"Invalid modification at {mod_pos} in {self.entity_name}"
			)
```

### Used in 
- [[Entity]]

### Uses


### Tags
#method