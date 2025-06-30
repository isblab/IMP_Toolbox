```python
def renumber_region_of_interest(
	self,
	region_of_interest: Dict,
):
	"""
	Offset the interacting region to the AF2/3 numbering.

	Interacting region defined by user is as per the original numbering 
	(UniProt in case of proteins). \n
	However, if the prediction is done on a fragment of the protein, the 
	numbering will be different. \n
	This function offsets the interacting region to the numbering of the 
	predicted structure. \n
	By default, the offset is assumed to be 0.

	Args:

		region_of_interest (Dict):
			Dictionary containing the region of interest for each chain.

	Returns:

		renumbered_region_of_interest (Dict):
			Dictionary containing the renumbered region of interest 
			for each chain.
	"""

	renumbered_region_of_interest = {}

	for chain_id in region_of_interest:

		start, end = region_of_interest[chain_id]

		if self.af_offset and chain_id in self.af_offset:

			start = start - (self.af_offset[chain_id][0] - 1)
			end = end - (self.af_offset[chain_id][0] - 1)

		renumbered_region_of_interest[chain_id] = [start, end]

	return renumbered_region_of_interest
```

### Used in


### Uses


### Tags
#method 