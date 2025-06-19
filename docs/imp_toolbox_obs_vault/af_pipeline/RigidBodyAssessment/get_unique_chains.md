```python
def get_unique_chains(self):
	"""Get unique chains in the rigid body.

	Returns:
		unique_chains (list): List of unique chain IDs in the rigid body.
	"""

	unique_chains = [
		chain_id
		for chain_id in self.rb_dict.keys()
		if len(self.rb_dict[chain_id]) > 0
	]

	return unique_chains
```

### Used in
[[RigidBodyAssessment]]

### Uses


### Tags
#method 