```python
def get_per_chain_avg_plddt(self):
	""" Get the average pLDDT score for each chain.

	Returns:
		per_chain_avg_plddt (dict): A dictionary where keys are chain IDs and values are the average pLDDT scores for that chain.
	"""

	return {
		chain_id: np.mean(plddt_scores)
		for chain_id, plddt_scores in self.per_chain_plddt.items()
	}
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 