```python
def create_af_sequence(self)-> Dict[str, Any]:
	"""Create an AF sequence dictionary

	Returns:
		af_sequence_dict (dict): AF sequence dictionary in the following format:
		- for proteinChain:
			{
				"proteinChain": {
					"sequence": "AAAA",
					"glycans": [... ],
					"modifications": [... ],
					"count": 1
				}
			}
		- for dnaSequence or rnaSequence:
			{
				"dnaSequence"("rnaSequence"): {
					"sequence": "ATCG",
					"modifications": [... ],
					"count": 1
			}
		- for ligand or ion:
			{
				"ligand"("ion"): {
					"ligand": "ATP",
					"count": 1
				}
			}
	"""

	if self.type == "proteinChain":
		af_sequence_dict = {
			self.type: {
				"sequence": self.real_sequence,
				"glycans": self.glycans,
				"modifications": self.modifications,
				"count": self.count,
			}
		}

		af_sequence_dict[self.type].update(self.template_settings)

	elif self.type in ["dnaSequence", "rnaSequence"]:
		af_sequence_dict = {
			self.type: {
				"sequence": self.real_sequence,
				"modifications": self.get_modifications(),
				"count": self.count,
			}
		}

	elif self.type in ["ligand", "ion"]:
		af_sequence_dict = {
			self.type: {self.type: self.name, "count": self.count}
		}

	return af_sequence_dict
```

### Used in
- [[update_af_sequences]]

### Uses
- [[get_modifications]]

### Tags
#method