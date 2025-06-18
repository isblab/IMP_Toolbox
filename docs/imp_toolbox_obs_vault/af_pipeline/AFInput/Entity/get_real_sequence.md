```python
def get_real_sequence(self)-> str:
	"""Get the real sequence of the entity
	- For proteinChain, get the sequence from the protein_sequences
	- For dnaSequence and rnaSequence, get the sequence from 
		the nucleic_acid_sequences

	Raises:
		Exception: Could not find the entity sequence

	Returns:
		str: amino acid or nucleic acid sequence of the entity
	"""

	real_sequence = None

	if self.entity_type == "proteinChain":

		try:
			uniprot_id = self.entities_map[self.entity_name]
			real_sequence = self.protein_sequences[uniprot_id]

		except KeyError:
			try:
				real_sequence = self.protein_sequences[self.entity_name]
			except KeyError:
				raise Exception(
					f"Could not find the entity sequence for {self.entity_name}"
				)

	elif self.entity_type in ["dnaSequence", "rnaSequence"]:

		try:
			nucleic_acid_id = self.entities_map[self.entity_name]
			real_sequence = self.nucleic_acid_sequences[nucleic_acid_id]

		except KeyError:
			try:
				real_sequence = self.nucleic_acid_sequences[self.entity_name]
			except KeyError:
				raise Exception(
					f"Could not find the entity sequence for {self.entity_name}"
				)

	return real_sequence
```

### Used in
- [[fill_up_entity]]

### Uses


### Tags
#method