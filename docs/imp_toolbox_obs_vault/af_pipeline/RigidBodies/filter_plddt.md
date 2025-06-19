```python
def filter_plddt(
	self,
	rb_dict: dict,
	patch_threshold: int = 0,
):
	"""Filter the residues in the rigid bodies based on the pLDDT cutoff.
	- If the pLDDT score of a residue is less than the cutoff, it is removed from the rigid body.
	Args:
		rb_dict (dict): dictionary of rigid bodies
		patch_threshold (int): minimum number of contiguous residues for which the pLDDT score is above the cutoff

	Returns:
		rb_dict (dict): dictionary of rigid bodies with residues filtered based on the pLDDT cutoff
	"""

	# Filter the residues in each chain in the rigid body based on the pLDDT cutoff
	for chain_id, rb_res_num_list in rb_dict.items():

		confident_residues = []

		# sorted list of residue numbers in the rigid body
		rb_res_num_arr = np.array(sorted(rb_res_num_list))

		# sorted list of residue indices in the rigid body
		plddt_res_num_arr = np.array([self.num_to_idx[chain_id][res_num] for res_num in rb_res_num_list])

		# True/False array based on the pLDDT cutoff
		# for e.g. plddt_arr = [70, 78, 90, 65, 65, 80, 90]
		# tf_plddt_filtered = [True, True, True, False, False, True, True] for cutoff = 70
		if chain_id in self.idr_chains:
			tf_plddt_filtered = np.array(self.plddt_list)[plddt_res_num_arr] >= self.plddt_cutoff_idr
		else:
			tf_plddt_filtered = np.array(self.plddt_list)[plddt_res_num_arr] >= self.plddt_cutoff

		# Convert the pLDDT scores to True/False based on the threshold
		# for e.g. if arr = [True, False, False, False, True, True] and threshold = 3
		# the output will be [True, True, True, True, True, True]
		tf_plddt_filtered = convert_false_to_true(
			arr=tf_plddt_filtered,
			threshold=patch_threshold,
		)

		# Get the residue numbers of the confident residues
		confident_residues = rb_res_num_arr[tf_plddt_filtered]
		rb_dict[chain_id] = confident_residues.tolist()

	# Remove chains which have no confident residues
	empty_chains = []

	for chain_id, confident_residues in rb_dict.items():
		if len(confident_residues) == 0:
			empty_chains.append(chain_id)

	for chain_id in empty_chains:
		del rb_dict[chain_id]

	return rb_dict
```

### Used in
- [[predict_domains]]

### Uses
- [[convert_false_to_true]]

### Tags
#method 