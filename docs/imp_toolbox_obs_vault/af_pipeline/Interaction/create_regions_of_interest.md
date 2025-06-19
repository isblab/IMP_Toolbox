```python
def create_regions_of_interest(self):
	"""
	Create regions of interest for all possible chain pairs.

	Returns:
		regions_of_interest (list): list of regions of interest
	"""

	regions_of_interest = []
	token_chain_ids = self.token_chain_ids
	chain_pairs = set()

	for chain1 in set(token_chain_ids):
		for chain2 in set(token_chain_ids):
			if chain1 != chain2:
				pair = tuple(sorted([chain1, chain2]))
				chain_pairs.add(pair)

	chain_pairs = list(chain_pairs)

	for chain1, chain2 in chain_pairs:

		ch1_start = self.renumber.renumber_chain_res_num(
			chain_res_num=1,
			chain_id=chain1
		)
		ch1_end = self.renumber.renumber_chain_res_num(
			chain_res_num=self.lengths_dict[chain1],
			chain_id=chain1
		)
		ch2_start = self.renumber.renumber_chain_res_num(
			chain_res_num=1,
			chain_id=chain2
		)
		ch2_end = self.renumber.renumber_chain_res_num(
			chain_res_num=self.lengths_dict[chain2],
			chain_id=chain2
		)

		region_of_interest = {
			chain1: [ch1_start, ch1_end],
			chain2: [ch2_start, ch2_end],
		}

		# region_of_interest = self.renumber.renumber_region_of_interest(
		#     region_of_interest=region_of_interest,
		# )

		regions_of_interest.append(region_of_interest)

	return regions_of_interest
```

### Used in


### Uses
- [[renumber_chain_res_num]]

### Tags
#method #user_function 