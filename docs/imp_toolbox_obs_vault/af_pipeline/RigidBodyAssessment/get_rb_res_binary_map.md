```python
def get_rb_res_binary_map(self, lengths_dict):
	"""Get a binary map of residues in the rigid body.

	Returns:
		rb_res_binary_map (np.ndarray): A binary map of residues in the rigid body.
		The shape is (total_length, total_length) where total_length is the sum of lengths of all chains.
		The value is 1 if the residue is part of the rigid body, 0 otherwise.
	"""

	total_len = lengths_dict.get("total", 0)
	rb_res_binary_map = np.zeros((total_len, total_len), dtype=int)
	all_rb_interface_res_idxs = []

	for chain_id, res_list in self.rb_dict.items():

		res_idxs = [self.num_to_idx[chain_id][res_num] for res_num in res_list]
		all_rb_interface_res_idxs.extend(res_idxs)

	all_rb_interface_res_idxs = np.unique(all_rb_interface_res_idxs)

	rb_res_binary_map[
		np.ix_(all_rb_interface_res_idxs, all_rb_interface_res_idxs)
	] = 1

	return rb_res_binary_map
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 