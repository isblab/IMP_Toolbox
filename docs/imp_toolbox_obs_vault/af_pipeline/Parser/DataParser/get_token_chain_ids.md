```python
@staticmethod
def get_token_chain_ids(data: Dict) -> list:
	"""Get the tokenized chain IDs from the data dict.

	This is specific to AF3: "token_chain_ids" key. \n
	If data is AF2, None is returned. \n
	However, similar information can be obtained from the structure file. \n
	see :py:meth:`Parser.StructureParser.get_token_chain_res_ids`.

	Args:

		data (Dict):
			data dict from the data file.

	Returns:

		token_chain_ids (list):
			tokenized chain IDs.
	"""

	if "token_chain_ids" in data:
		token_chain_ids = data["token_chain_ids"]

	else:
		warnings.warn(
			"""
			Chain IDs not found, data file might be AF2.
			Structure file is required for AF2.
			"""
		)
		token_chain_ids = None

	return token_chain_ids
```

### Used in
- [[get_attributes]]

### Uses


### Tags
#method 