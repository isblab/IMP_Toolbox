```python
def get_interface_res_pairs(
	self,
	contact_map: np.ndarray,
):
	""" Get interface residue pairs from the contact map.

	Args:
		contact_map (np.ndarray): A binary contact map where 1 indicates a contact between residues and 0 indicates no contact.

	Returns:
		interface_res_pairs (defaultdict): A dictionary where keys are chain pairs (tuples) and values are lists of residue index pairs.
	"""

	interface_res_pairs = defaultdict(list)
	contacting_res_indices = np.argwhere(contact_map == 1)

	for chain1, chain2 in tqdm(self.chain_pairs):

		res1_list = self.rb_dict[chain1]
		res2_list = self.rb_dict[chain2]

		res1_idxs = [self.num_to_idx[chain1][res_num] for res_num in res1_list]
		res2_idxs = [self.num_to_idx[chain2][res_num] for res_num in res2_list]

		mask1 = np.isin(contacting_res_indices[:, 0], res1_idxs)
		mask2 = np.isin(contacting_res_indices[:, 1], res2_idxs)

		mask = mask1 & mask2
		contacting_res_pairs = set(
			map(tuple, contacting_res_indices[mask])
		)

		if len(contacting_res_pairs) > 0:
			interface_res_pairs[(chain1, chain2)].extend(
				list(contacting_res_pairs)
			)

	return interface_res_pairs
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 