```python
def update_af_sequences(self):
	"""Update the AF sequences
	- For each entity, create an AFSequence object
	- Get the name fragment for each entity 
		(used in job name if job name is not provided)
	"""

	# add af_sequence for each entity
	for entity_info in self.job_info["entities"]: 
		af_sequence = AFSequence(
			entity_info=entity_info,
			protein_sequences=self.protein_sequences,
			nucleic_acid_sequences=self.nucleic_acid_sequences,
			entities_map=self.entities_map,
		)
		af_sequence_dict = af_sequence.create_af_sequence()
		self.af_sequences.append(af_sequence_dict)

		self.name_fragments.append(af_sequence.get_name_fragment())
```

### Used in
- [[create_job]]

### Uses
- [[AFSequence]]
- [[create_af_sequence]]
- [[get_name_fragment]]

### Tags
#method