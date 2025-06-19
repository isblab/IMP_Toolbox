```python
def save_rigid_bodies(
	self,
	domains: list,
	output_dir: str,
	output_format: str = "txt",
	save_structure: bool = True,
	structure_file_type: str = "pdb",
	no_plddt_filter_for_structure: bool = False,
	pae_plot: bool = False,
	rb_assessment: dict | None = None,
	protein_chain_map: dict | None = None,
):
	""" Save the rigid bodies to a file and/or save the structure of the rigid bodies and assess the rigid bodies.
	- The rigid bodies are saved in a plain text format with the chain IDs and residue numbers.
	- The structure of the rigid bodies can be saved in PDB or CIF format. For rigid bodies with modifications, it is recommended to use PDB format.
	- The PAE plot can be saved to visualize the rigid bodies in the PAE matrix.
	- The rigid bodies can be assessed based on the interface residues, number of contacts, interface PAE and pLDDT, average PAE and plDDT and minimum PAE.
	- The assessment is saved in an Excel file.

	parameters for rigid body assessment:\n
	- `as_average`: \n
	whether to report only the average of assessment metric to the output file. Defaults to False. \n
	- `symmetric_pae`: \n
	whether to report a single average PAE value or assymetric PAE value for PAE assessment metrics. Defaults to False. \n

	Args:
		domains (list): list of rigid bodies, where each rigid body is a dictionary with chain IDs as keys and residue numbers as values.
		output_dir (str): Directory to save the output files.
		output_format (str, optional): Defaults to "txt". ("txt" or "csv")
		save_structure (bool, optional): Whether to save the structure of the rigid bodies. Defaults to True.
		structure_file_type (str, optional): File type to save the structure. Defaults to "pdb". ("pdb" or "cif")
		no_plddt_filter_for_structure (bool, optional): Whether to save the structure without filtering based on pLDDT. Defaults to False.
		pae_plot (bool, optional): Whether to save the PAE plot for the rigid bodies. Defaults to False.
		rb_assessment (dict | None, optional): Dictionary containing parameters for rigid body assessment.
	"""

	dir_name = os.path.basename(self.struct_file_path).split(".")[0]
	output_dir = os.path.join(output_dir, dir_name)

	os.makedirs(output_dir, exist_ok=True)

	file_name = (
		os.path.basename(self.struct_file_path).split(".")[0] + "_rigid_bodies"
	)

	protein_chain_map = self.extract_protein_chain_mapping(protein_chain_mapping=protein_chain_map)

	##################################################

	# txt or csv output format
	if output_format == "txt":
		file_name += ".txt"
		output_path = os.path.join(output_dir, file_name)

		with open(output_path, "w") as f:

			for idx, rb_dict in enumerate(domains):
				f.write(f"Rigid Body {idx}\n")

				for chain_id, res_list in rb_dict.items():

					protein_name = protein_chain_map.get(chain_id, None)

					if len(res_list) > 0:
						if protein_name:
							f.write(
								f"{protein_name}_{chain_id}: {get_key_from_res_range(res_range=res_list)}\n"
							)
						else:
							f.write(
								f"{chain_id}:{get_key_from_res_range(res_range=res_list)}\n"
							)

				f.write("\n")

	elif output_format == "csv":
		file_name += ".csv"
		output_path = os.path.join(output_dir, file_name)

		rows = []
		for idx, rb_dict in enumerate(domains):
			for chain_id, res_list in rb_dict.items():
				if len(res_list) > 0:
					protein_name = protein_chain_map.get(chain_id, None)

					if protein_name:
						rows.append({
							"Rigid Body": idx,
							"Chain ID": chain_id,
							"Protein Name": protein_name,
							"Residues": get_key_from_res_range(res_range=res_list),
						})
					else:
						rows.append({
							"Rigid Body": idx,
							"Chain ID": chain_id,
							"Residues": get_key_from_res_range(res_range=res_list),
						})

		df = pd.DataFrame(rows)
		df.to_csv(output_path, index=False)

	else:
		raise ValueError(
			f"Invalid output format: {output_format}. Use 'txt' or 'csv'."
		)

	##################################################
	# Save the structure of the rigid bodies
	if save_structure:

		if structure_file_type == "cif":
			warnings.warn(
				"""
				Protein or nucleotide modifications are stored as HETATM for which sequence connectivity
				is lost in CIF format. \n
				Please use PDB format to save the structure with modifications.
				"""
			)

		# Renumber the structure to match the actual sequence numbering if af_offset is provided
		structure = self.renumber.renumber_structure(
			structure=self.structureparser.structure,
		)

		for idx, rb_dict in enumerate(domains):

			# In the following case, the txt or csv ouput will have pLDDT filtered residues
			# but, the structure file will ignore this filter
			# use this flag when you don't want missing residues in the structure file
			if no_plddt_filter_for_structure:
				for chain_id, res_list in rb_dict.items():
					if len(res_list) > 0:
						res_list = fill_up_the_blanks(res_list)
						rb_dict[chain_id] = res_list

			output_path = os.path.join(output_dir, f"rigid_body_{idx}.{structure_file_type}")

			save_structure_obj(
				structure=structure,
				out_file=output_path,
				res_select_obj=ResidueSelect(rb_dict),
				save_type=structure_file_type,
				preserve_header_footer=False,
			)

	##################################################
	# Save the PAE plot for the rigid bodies
	# the region of the PAE matrix corresponding to the rigid bodies will be highlighted
	if pae_plot:
		for rb_idx, rb_dict in enumerate(domains):

			# patches are the highlighted rectangles in the PAE matrix
			patches = []

			for chain_id1, res_list1 in rb_dict.items():

				for chain_id2, res_list2 in rb_dict.items():

					res_idxs_1 = [
						self.num_to_idx[chain_id1][res_num] for res_num in res_list1
					]
					res_idxs_2 = [
						self.num_to_idx[chain_id2][res_num] for res_num in res_list2
					]
					res_idx_range_1 = get_key_from_res_range(res_range=res_idxs_1, as_list=True)
					res_idx_range_2 = get_key_from_res_range(res_range=res_idxs_2, as_list=True)

					for res_idx_1 in res_idx_range_1:

						for res_idx_2 in res_idx_range_2:

							if "-" in res_idx_1 and "-" in res_idx_2:
								xy_ = (int(res_idx_2.split("-")[0]), int(res_idx_1.split("-")[0])) # xy (0,0) coordinates for the rectangle
								h_ = int(res_idx_1.split("-")[1]) - int(res_idx_1.split("-")[0]) + 1 # patch height
								w_ = int(res_idx_2.split("-")[1]) - int(res_idx_2.split("-")[0]) + 1 # patch width

								if h_ > 0 and w_ > 0:
									patches.append(
										[xy_, h_, w_]
									)

			fig = plt.figure(figsize=(20, 20))
			plt.rcParams['font.size'] = 16
			plt.rcParams['axes.titlesize'] = 28
			plt.rcParams['axes.labelsize'] = 22
			plt.rcParams['xtick.labelsize'] = 13
			plt.rcParams['ytick.labelsize'] = 13
			plt.imshow(
				self.pae,
				# cmap="Greens_r",
				cmap="Greys_r",
				vmax=31.75,
				vmin=0,
				interpolation="nearest",
				)

			for xy, h, w in patches:
				rect = matplotlib.patches.Rectangle(
					xy,
					w,
					h,
					linewidth=0,
					# edgecolor="green",
					facecolor="lime",
					alpha=0.5,
				)
				plt.gca().add_patch(rect)

			cumu_len = 0
			ticks = []
			ticks_labels = []

			for chain_id, p_length in self.lengths_dict.items():
				if chain_id != "total":
					cumu_len += p_length

					if cumu_len != self.pae.shape[1]:
						plt.axhline(y=cumu_len, color='red', linestyle='--', linewidth=0.75)
						plt.axvline(x=cumu_len, color='red', linestyle='--', linewidth=0.75)

					if self.af_offset is not None:

						ticks_labels.extend(["\n" + f"{self.af_offset[chain_id][0]}" , f"{self.af_offset[chain_id][1]}" + "\n" ])
						ticks.extend([cumu_len-p_length, cumu_len]) if cumu_len-p_length not in ticks else ticks.extend([cumu_len-p_length+1, cumu_len])

					else:
						ticks_labels.extend(["\n" + "1", f"{self.lengths_dict[chain_id]}" + "\n"])
						ticks.extend([cumu_len-p_length, cumu_len]) if cumu_len-p_length not in ticks else ticks.extend([cumu_len-p_length+1, cumu_len])

			plt.xlim(0, self.pae.shape[0])
			plt.ylim(0, self.pae.shape[1])

			plt.gca().invert_yaxis()
			plt.yticks(ticks, ticks_labels)

			plt.xticks(ticks, ticks_labels, rotation=90, ha='center')
			plt.title(f"Predicted aligned error (PAE)", pad=20)
			plt.xlabel("Scored residue")
			plt.ylabel("Aligned residue")

			ax = plt.gca()

			divider = make_axes_locatable(ax)
			cax = divider.append_axes("bottom", size="5%", pad=1.2)
			plt.colorbar(
				label="Predicted Alignment Error (PAE)",
				orientation="horizontal",
				cax=cax,
			)

			plt.savefig(os.path.join(output_dir, f"rigid_body_{rb_idx}.png"), transparent=True)
			plt.close(fig)

	##################################################
	# Save the assessment of rigid bodies
	if rb_assessment:

		_start = time.time()
		print("Assessing rigid bodies...")

		assessment_file_name = (
			os.path.basename(self.struct_file_path).split(".")[0] + "_rb_assessment.xlsx"
		)
		save_path = os.path.join(output_dir, assessment_file_name)

		coords = np.array(self.coords_list)

		contact_map = get_interaction_map(
			coords1=coords,
			coords2=coords,
			contact_threshold=8,
			map_type="contact",
		)

		for rb_idx, rb_dict in enumerate(domains):

			rb_save_path = save_path.replace(
				".xlsx", f"_rb_{rb_idx}.xlsx"
			)

			rb_assess = RigidBodyAssessment(
				rb_dict=rb_dict,
				num_to_idx=self.num_to_idx,
				idx_to_num=self.idx_to_num,
				contact_map=contact_map,
				plddt_list=self.plddt_list,
				pae=self.pae,
				lengths_dict=self.lengths_dict,
				save_path=rb_save_path,
				symmetric_pae=rb_assessment.get("symmetric_pae", False),
				as_average=rb_assessment.get("as_average", False),
				idr_chains=self.idr_chains,
				protein_chain_map=protein_chain_map,
			)

			rb_assess.save_rb_assessment()

		print(f"Time taken to save rigid body assessment: {time.time() - _start:.2f} seconds")
```

### Used in


### Uses
- [[extract_protein_chain_mapping]]
- [[get_key_from_res_range]]
- [[renumber_structure]]
- [[fill_up_the_blanks]]
- [[save_structure_obj]]
- [[get_interaction_map]]
- [[save_rb_assessment]]

### Tags
#method #user_function 