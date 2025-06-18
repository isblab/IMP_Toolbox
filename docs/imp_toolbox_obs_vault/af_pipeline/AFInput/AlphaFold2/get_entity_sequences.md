```python
def get_entity_sequences(
	self,
	ranges: List[Tuple[int, int]],
	headers: List[str],
) -> List[str]:
	"""Get the entity sequences

	First try to get the sequence from the protein_sequences dictionary. \n
	If not found, try to get the sequence from the proteins dictionary. \n
	If not found, raise an exception.

	If a range is provided, get the sequence within the range.

	Args:
		ranges (list): [start, end] of the entities
		headers (list): fasta headers

	Returns:
		sequences (list): list of entity sequences
	"""

	sequences = []

	for header in headers:
		try:
			sequences.append(self.protein_sequences[header])
		except KeyError:
			try:
				sequences.append(
					self.protein_sequences[self.entities_map[header]]
				)
			except KeyError:
				raise Exception(
					f"Could not find the entity sequence for {header}"
				)

	for i, range_ in enumerate(ranges):
		if range_:
			start, end = range_
			sequences[i] = sequences[i][start - 1 : end]

	return sequences
```

### Used in
- [[generate_job_entities]]

### Uses


### Tags
#method