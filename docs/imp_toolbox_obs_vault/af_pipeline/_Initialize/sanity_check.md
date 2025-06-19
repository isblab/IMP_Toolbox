```python
def sanity_check(self):
	"""
	Perform sanity checks on the input data. \n
	If structure file path is not provided, the input data file should be in AF3 format.
	"""

	error_statement = "Input data file needs to be in AF3 format if structure path is not provided."

	if not self.lengths_dict:
		raise Exception(f"No chain lengths found. {error_statement}")
```

### Used in
- [[get_attributes]]

### Uses


### Tags
#method 