```python
def get_per_chain_average_iplddt(self):
	""" Get the average ipLDDT score for each chain.

	Returns:
		per_chain_avg_iplddt (dict): A dictionary where keys are chain IDs and values are the average ipLDDT scores for that chain.
	"""

	return {
		chain_id: np.mean(list(iplddt_scores.values()))
		for chain_id, iplddt_scores in self.per_chain_iplddt.items()
	}
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 