```python
def residue_map(
	self,
	token_chain_ids: list,
	token_res_ids: list
):
	"""Create a map of residue indices to residue numbers and vice-versa.

	res_idx is essentially token index. \n
	res_num is the residue number. \n
	res_num = res_idx + 1 if af_offset is not provided. \n
	res_num = res_idx + af_offset if af_offset is provided. \n
	af_offset informs what is the starting residue number for each chain.

	Args:

		token_chain_ids (list):
			Tokenized chain IDs.

		token_res_ids (list):
			Tokenized residue IDs.

	Returns:

		tuple (idx_to_num, num_to_idx):
			idx_to_num (Dict):
				Dictionary mapping residue indices to residue numbers.

			num_to_idx (Dict):
				Dictionary mapping residue numbers to residue indices.
	"""

	idx_to_num = {}
	num_to_idx = defaultdict(dict)

	for res_idx, (chain_id, res_num) in enumerate(
		zip(token_chain_ids, token_res_ids)
	):

		res_num = self.renumber_chain_res_num(
			chain_res_num=res_num,
			chain_id=chain_id,
		)

		idx_to_num[res_idx] = {
			"chain_id": chain_id,
			"res_num": res_num,
		}

		num_to_idx[chain_id][res_num] = res_idx

	return idx_to_num, num_to_idx
```

### Used in
- [[get_attributes]]

### Uses
- [[renumber_chain_res_num]]

### Tags
#method 