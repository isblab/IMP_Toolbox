```python
def get_chain_pairs(self):
	"""Get all unique chain pairs in the rigid body.

	Returns:
		chain_pairs (list): List of tuples containing unique chain pairs.
		Each tuple contains two chain IDs.
	"""

	chain_pairs = list(combinations(self.unique_chains, 2))

	return [tuple(pair) for pair in chain_pairs]
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 