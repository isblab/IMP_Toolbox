```python
def get_entity_range(self)-> Tuple[int, int]:
	"""Get the range of the entity
	what part of the sequence to use? (defined by start and end)
	- If range is provided, use that
	- If no range is provided, use the full sequence
	- If no sequence is found (e.g. ligand or ion), use a range of [1, 1]

	Returns:
		tuple: start and end of the entity
	"""

	if "range" in self.entity_info:

		assert (
			len(self.entity_info["range"]) == 2
		), "Invalid range; must be a list of two integers (start and end)"

		start, end = self.entity_info["range"]

	else:
		start, end = 1, 1

		if self.real_sequence:
			start, end = 1, len(self.real_sequence)

	return start, end
```


### Used in
- [[fill_up_entity]]

### Uses


### Tags
#method