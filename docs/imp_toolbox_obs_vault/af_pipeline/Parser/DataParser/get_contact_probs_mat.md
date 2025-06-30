```python
@staticmethod
def get_contact_probs_mat(data: Dict):
	"""Get the contact probabilities from the data dict.

	Args:

		data (Dict):
			data dict from the data file.

	Returns:

		contact_probs_mat (np.array):
			contact probabilities matrix from AlphaFold3 output.
	"""

	if "contact_probs" in data:
		contact_probs_mat = np.array(data["contact_probs"])

	else:
		warnings.warn(
			"Contact probabilities not found, data file might not be AF3."
		)
		contact_probs_mat = None

	return contact_probs_mat
```


### Used in
- [[get_attributes]]

### Uses


### Tags
#method 