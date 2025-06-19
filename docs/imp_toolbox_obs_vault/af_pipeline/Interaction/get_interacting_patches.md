```python
def get_interacting_patches(
	self,
	contact_map: np.array,
	region_of_interest: dict,
):
	"""This is a dirty implementation to get the interacting patches. \n
	This is a temporary solution until we find a better way to get interacting
	patches for the given contact map.

	Args:
		contact_map (np.array): binary contact map.
		region_of_interest (dict): region of interest for the protein pair.

	Returns:
		patches (dict): interacting patches for the given region of interest of the protein pair.
	"""

	patches = {}

	chain1, chain2 = region_of_interest.keys()
	p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

	if np.unique(contact_map).tolist() == [0]: # No interactions found.
		warnings.warn(
			f"No interacting patches found for {chain1}:{p1_region} and {chain2}:{p2_region}."
		)
		return patches

	patches_df = get_patches_from_matrix(
		matrix=contact_map,
		chain1=chain1,
		chain2=chain2
	)

	for patch_idx, patch in patches_df.iterrows():

		ch1_patch = patch[chain1]
		ch2_patch = patch[chain2]

		ch1_patch = sorted([int(x) for x in ch1_patch])
		ch2_patch = sorted([int(x) for x in ch2_patch])

		ch1_patch = np.array(ch1_patch) + region_of_interest[chain1][0]
		ch2_patch = np.array(ch2_patch) + region_of_interest[chain2][0]

		# patches[patch_idx] = {
		#     chain1: (
		#         f"{ch1_patch[0]}-{ch1_patch[-1]}"
		#         if len(ch1_patch) > 1
		#         else str(ch1_patch[0])
		#     ),
		#     chain2: (
		#         f"{ch2_patch[0]}-{ch2_patch[-1]}"
		#         if len(ch2_patch) > 1
		#         else str(ch2_patch[0])
		#     ),
		# }

		patches[patch_idx] = {
			chain1: np.array(ch1_patch),
			chain2: np.array(ch2_patch),
		}

	return patches
```

### Used in
- [[save_ppair_interaction]]

### Uses
- [[get_patches_from_matrix]]

### Tags
#method 