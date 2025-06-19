```python
def renumber_region_of_interest(
	self,
	region_of_interest: Dict,
):
	"""
	Offset the interacting region to the AF2/3 numbering. \n
	Interacting region defined by user is as per the original numbering (UniProt in case of proteins). \n
	However, if the prediction is done on a fragment of the protein, the numbering will be different. \n
	This function offsets the interacting region to the numbering of the predicted structure. \n
	By default, the offset is assumed to be 0.

	Args:
		region_of_interest (Dict): dict containing the region of interest for each chain.

	Returns:
		renumbered_region_of_interest (Dict): dict containing the renumbered region of interest for each chain.

	Example:
		consider a prediction involving proteins A (100 aa) and B (50 aa). \n
		prediction is done on a fragment of A (30-100) and B (10-50). \n
		so, user defines - \n
		af_offset = {'A': [30, 100], 'B': [10, 50]} \n
		region_of_interest = {'A': (30, 50), 'B': (20, 40)}

		renumbered_region_of_interest = {'A': (1, 21), 'B': (11, 31)} \n
		i.e. within the predicted structure, the functions in the Interaction class will look for
		interactions in the region of: 1-21 resdiue of A and 11-31 residues of B.
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