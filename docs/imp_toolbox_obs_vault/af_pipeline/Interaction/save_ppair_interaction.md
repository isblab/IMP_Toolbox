```python
def save_ppair_interaction(
	self,
	region_of_interest: Dict,
	save_plot: bool = False,
	plot_type: str = "static",
	p1_name: str | None = None,
	p2_name: str | None = None,
	concat_residues: bool = True,
	contact_probability: bool = True,
):
	"""Save the interacting patches for the given region of interest of the protein pair.

	Args:
		region_of_interest (Dict): Dictionary containing the chain IDs and the residue indices for the region of interest.
		save_plot (bool, optional): Outputs the plot if True. Defaults to False.
		plot_type (str, optional): Type of plot to be saved. Defaults to "static"; options: ["static", "interactive", "both"].
		p1_name (str, optional): Name of the first protein. Defaults to None.
		p2_name (str, optional): Name of the second protein. Defaults to None.
		concat_residues (bool, optional): Whether to concatenate the residues into residue ranges. Defaults to True.
		contact_probability (bool, optional): Whether to add contact probability column to the output. Defaults to True.
	"""

	chain1, chain2 = list(region_of_interest.keys())
	p1_region, p2_region = region_of_interest[chain1], region_of_interest[chain2]

	contact_map = self.get_confident_interaction_map(
		region_of_interest=region_of_interest
	)

	interacting_patches = self.get_interacting_patches(
		contact_map=contact_map,
		region_of_interest=region_of_interest,
	)

	if p1_name and p2_name:
		p_names = {
			chain1: p1_name,
			chain2: p2_name,
		}
		dir_name_to_replace = "_".join([p1_name, p2_name])
	else:
		p_names = {
			chain1: chain1,
			chain2: chain2,
		}
		dir_name_to_replace = None


	if len(interacting_patches) > 0:

		file_name = "_".join([
			f"{p_names[k]}_{k}:{v[0]}-{v[1]}" for k, v in region_of_interest.items()
		])

		if dir_name_to_replace:
			dir_name = os.path.basename(self.struct_file_path).split(".")[0]
			# self.output_dir = self.output_dir.replace(dir_name, dir_name_to_replace)

		os.makedirs(self.output_dir, exist_ok=True)

		save_map(
			contact_map=contact_map,
			avg_contact_probs_mat=self.avg_contact_probs_mat,
			patches=interacting_patches,
			chain1=chain1,
			chain2=chain2,
			p1_name=p_names[chain1],
			p2_name=p_names[chain2],
			p1_region=p1_region,
			p2_region=p2_region,
			out_file=os.path.join(self.output_dir, f"patches_{file_name}.html"),
			save_plot=save_plot,
			plot_type=plot_type,
			concat_residues=concat_residues,
			contact_probability=contact_probability,
			num_to_idx=self.num_to_idx,
			idx_to_num=self.idx_to_num,
		)
```

### Used in


### Uses
- [[get_confident_interaction_map]]
- [[get_interacting_patches]]
- [[save_map]]

### Tags
#user_function #method 