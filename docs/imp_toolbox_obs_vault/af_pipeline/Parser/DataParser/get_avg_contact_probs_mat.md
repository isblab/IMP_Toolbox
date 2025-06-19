```python
@warnings.deprecated
def get_avg_contact_probs_mat(self, contact_probs_mat: np.ndarray):
	"""
	Return the average contact probabilities matrix. \n

	Args:
		contact_probs_mat (np.array): contact probabilities matrix.

	Returns:
		avg_contact_probs_mat (np.array): average contact probabilities matrix.
	"""

	avg_contact_probs_mat = (contact_probs_mat + contact_probs_mat.T) / 2

	return avg_contact_probs_mat
```


### Used in
- [[get_attributes]]

### Uses


### Tags
#method 