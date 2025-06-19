```python
@warnings.deprecated
def get_avg_pae(self, pae: np.ndarray):
	"""
	Return the average PAE matrix. \n

	Args:
		pae (np.array): PAE matrix.

	Returns:
		avg_pae (np.array): average PAE matrix.
	"""

	avg_pae = (pae + pae.T) / 2

	return avg_pae
```


### Used in
- [[get_attributes]]

### Uses


### Tags
#method 