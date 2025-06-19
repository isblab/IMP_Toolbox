```python
def get_num_contacts(self):
	"""Get the number of contacts for each chain pair.

	Returns:
		num_contacts (defaultdict): A dictionary where keys are chain pairs (tuples) and values are the number of contacts.
		Each key is a tuple of two chain IDs, and the value is the count of contacts between those chains.
	"""

	num_contacts = defaultdict(int)

	for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

		num_contacts[chain_pair] = len(interacting_res_pairs)

	return num_contacts
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 