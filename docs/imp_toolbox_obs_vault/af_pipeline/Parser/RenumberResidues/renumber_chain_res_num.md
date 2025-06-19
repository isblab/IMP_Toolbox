```python
def renumber_chain_res_num(
	self,
	chain_res_num: int,
	chain_id: str,
):
	"""
	Renumber the residue number based on the offset.

	Args:
		chain_res_num (int): residue number
		af_offset (dict): offset dictionary

	Returns:
		chain_res_num (int): renumbered residue number
	"""

	if self.af_offset and chain_id in self.af_offset:
		chain_res_num += self.af_offset[chain_id][0] - 1

	return chain_res_num
```

### Used in
- [[renumber_structure]]
- [[residue_map]]

### Uses


### Tags
#method 