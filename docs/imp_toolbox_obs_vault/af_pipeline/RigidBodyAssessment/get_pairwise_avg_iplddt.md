```python
def get_pairwise_avg_iplddt(self):
	""" Get the average ipLDDT for each chain pair.

	Returns:
		pairwise_avg_iplddt (defaultdict): A dictionary where keys are chain pairs (tuples) and values are dictionaries containing average ipLDDT values for each chain in the pair.
	"""

	pairwise_avg_iplddt = defaultdict(dict)

	for chain_pair, interacting_res_pairs in self.interface_res_pairs.items():

		chain1, chain2 = chain_pair

		iplddt1_values = [
			self.per_chain_iplddt[chain1].get(res1_idx, np.nan)
			for res1_idx, res2_idx in interacting_res_pairs
		]

		iplddt2_values = [
			self.per_chain_iplddt[chain2].get(res2_idx, np.nan)
			for res1_idx, res2_idx in interacting_res_pairs
		]

		pairwise_avg_iplddt[chain_pair][chain1] = np.mean(iplddt1_values)
		pairwise_avg_iplddt[chain_pair][chain2] = np.mean(iplddt2_values)

	return pairwise_avg_iplddt
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 