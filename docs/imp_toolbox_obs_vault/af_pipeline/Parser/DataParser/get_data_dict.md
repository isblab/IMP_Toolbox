```python
def get_data_dict(self) -> Dict:
	"""Parse the AF2/3 data file.

	Args:

		data_file_path (str):
			path to the data file.

	Returns:

		data (Dict):
			data dict from the data file.
	"""

	ext = os.path.splitext(self.data_file_path)[1]

	# AF2 data file
	if "pkl" in ext:
		with open(self.data_file_path, "rb") as f:
			data = pkl.load(f)

	# AF3 data file
	elif "json" in ext:
		with open(self.data_file_path, "r") as f:
			data = json.load(f)

		if isinstance(data, list):
			data = data[0]

	else:
		raise Exception(
			"Incorrect file format.. Suported .pkl/.json only."
		)

	return data
```

### Used in
- [[get_attributes]]

### Uses


### Tags
#method 