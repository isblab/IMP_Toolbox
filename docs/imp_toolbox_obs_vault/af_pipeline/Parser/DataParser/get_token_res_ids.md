```python
def get_token_res_ids(self, data: Dict) -> list:
	""" Get the tokenized residue IDs from the data dict. \n
	This is specific to AF3: "token_res_ids" key. \n
	If data is AF2, None is returned. \n
	However, similar information can be obtained from the structure file. \n
	see :py:meth:`Parser.AfParser.StructureParser.get_token_chain_res_ids`.

	Args:
		data (Dict): data dict from the data file.

	Returns:
		token_res_ids (list): tokenized residue IDs.
	"""

	if "token_res_ids" in data:
		token_res_ids = data["token_res_ids"]

	else:
		warnings.warn(
			"""
			Residue IDs not found, data file might be AF2.
			Structure file is required for AF2.
			"""
		)
		token_res_ids = None

	return token_res_ids
```

### Used in 
- [[get_attributes]]

### Uses


### Tags
#method 