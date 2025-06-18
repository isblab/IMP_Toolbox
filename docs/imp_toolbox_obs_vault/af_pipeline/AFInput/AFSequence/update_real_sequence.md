```python
def update_real_sequence(self)-> str:
	"""Update the real sequence of the entity

	real sequence is:
	- amino acid sequence for proteinChain
	- and nucleic acid sequence for dnaSequence and rnaSequence

	Returns:
		real_sequence (str): amino acid or nucleic acid sequence of the entity
	"""

	real_sequence = self.real_sequence
	start, end = self.start, self.end

	if self.type in ["proteinChain", "dnaSequence", "rnaSequence"]:
		real_sequence = real_sequence[start - 1 : end]

	return real_sequence
```

### Used in
- [[AFSequence]]

### Uses


### Tags
#method