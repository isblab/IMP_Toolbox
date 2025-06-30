```python
@staticmethod
def get_chain_lengths(
	token_chain_ids: list,
) -> Dict:
	"""Get the chain lengths.

	lengths_dict is a dictionary containing the chain lengths. \n
	{chain_id: length} \n
	"total" is the total length of the system. \n
	For example, if the system has 2 chains A and B, \n
	lengths_dict = {"A": 100, "B": 50, "total": 150} \n

	Args:

		token_chain_ids (list):
			tokenized chain IDs.

	Returns:

		lengths_dict (Dict):
			dict containing the chain lengths.

	Example:

		>>> token_chain_ids = ['A', 'A', 'B', 'B', 'B']
		>>> lengths_dict = AfParser.get_chain_lengths(token_chain_ids)
		>>> print(lengths_dict)
		{'total': 5, 'A': 2, 'B': 3}
	"""

	lengths_dict = {}
	lengths_dict["total"] = 0

	for chain_id in token_chain_ids:
		if chain_id not in lengths_dict:
			lengths_dict[chain_id] = 1
		else:
			lengths_dict[chain_id] += 1
		lengths_dict["total"] += 1

	return lengths_dict
```

### Used in
- [[_Initialize]]

### Uses


### Tags
#method 