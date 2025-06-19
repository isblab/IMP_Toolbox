```python
def get_pae(self, data: Dict):
	"""
	Return the PAE matrix from the data dict. \n

	Args:
		data (Dict): data dict from the data file.

	Returns:
		pae (np.array): PAE matrix.
	"""

	# For AF2.
	if "predicted_aligned_error" in data:
		pae = np.array(data["predicted_aligned_error"])

	# For AF3.
	elif "pae" in data:
		pae = np.array(data["pae"])

	else:
		raise Exception("PAE matrix not found...")

	return pae
```

### Used in
- [[get_attributes]]

### Uses


### Tags
#method 
