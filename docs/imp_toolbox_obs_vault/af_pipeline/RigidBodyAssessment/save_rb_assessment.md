```python
def save_rb_assessment(self):
	""" Save the assessment of the rigid bodies to an Excel file.

	The assessment includes:
	- Per chain assessment: Average pLDDT, Average iLDDT, Number of interface residues, Chain type (IDR or R)
	- Per chain pair assessment: Number of interface residues, Number of contacts, Average PAE, Average iPAE, Minimum PAE, Average iLDDT for each chain, Chain type (IDR or R) for each chain
	- Overall assessment: Average pLDDT, Average iLDDT, Number of interface residues, Chain type (IDR or R)

	The assessment is saved in an Excel file with three sheets:
	- "Chain Wise Assessment": Contains per chain assessment data.
	- "Chain Pairwise Assessment": Contains per chain pair assessment data.
	- "Overall Assessment": Contains overall assessment data.
	"""

	chain_wise_assessment_rows = []
	chain_pairwise_assessment_rows = []
	overall_assessment_rows = []

	for chain_id in self.unique_chains:
		chain_wise_assessment_rows.append({
			"Chain ID": chain_id,
			"Protein Name": self.protein_chain_map.get(chain_id, None),
			"Average pLDDT": self.per_chain_avg_plddt[chain_id],
			"Average ipLDDT": self.per_chain_avg_iplddt.get(chain_id, np.nan),
			"Number of Interface Residues": len(self.per_chain_interface_residues[chain_id]),
			"Chain Type": "IDR" if chain_id in self.idr_chains else "R",
		})

	for chain_pair in self.chain_pairs:
		chain1, chain2 = chain_pair

		if self.as_average:

			if self.symmetric_pae:
				chain_pairwise_assessment_rows.append({
					"Chain Pair": f"{chain1}-{chain2}",
					"Protein Name 1": self.protein_chain_map.get(chain1, "Unknown"),
					"Protein Name 2": self.protein_chain_map.get(chain2, "Unknown"),
					"Number of Interface Residues": self.num_interface_residues[chain_pair],
					"Number of Contacts": self.num_contacts[chain_pair],
					"Average PAE": self.pairwise_avg_pae[chain_pair],
					"Average iPAE": self.pairwise_avg_ipae[chain_pair],
					"Minimum PAE": self.pairwise_min_pae[chain_pair],
					"Average ipLDDT chain1": self.pairwise_avg_iplddt[chain_pair].get(chain1, np.nan),
					"Average ipLDDT chain2": self.pairwise_avg_iplddt[chain_pair].get(chain2, np.nan),
					"Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
					"Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
				})
			else:
				chain_pairwise_assessment_rows.append({
					"Chain Pair": f"{chain1}-{chain2}",
					"Protein Name 1": self.protein_chain_map.get(chain1, "Unknown"),
					"Protein Name 2": self.protein_chain_map.get(chain2, "Unknown"),
					"Number of Interface Residues": self.num_interface_residues[chain_pair],
					"Number of Contacts": self.num_contacts[chain_pair],
					"Average PAE ij": self.pairwise_avg_pae[chain_pair]["ij"],
					"Average PAE ji": self.pairwise_avg_pae[chain_pair]["ji"],
					"Average iPAE ij": self.pairwise_avg_ipae[chain_pair]["ij"] if chain_pair in self.pairwise_avg_ipae else np.nan,
					"Average iPAE ji": self.pairwise_avg_ipae[chain_pair]["ji"] if chain_pair in self.pairwise_avg_ipae else np.nan,
					"Minimum PAE ij": self.pairwise_min_pae[chain_pair]["ij"] if chain_pair in self.pairwise_min_pae else np.nan,
					"Minimum PAE ji": self.pairwise_min_pae[chain_pair]["ji"] if chain_pair in self.pairwise_min_pae else np.nan,
					"Average ipLDDT chain1": self.pairwise_avg_iplddt[chain_pair].get(chain1, np.nan),
					"Average ipLDDT chain2": self.pairwise_avg_iplddt[chain_pair].get(chain2, np.nan),
					"Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
					"Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
				})

		else:

			if self.symmetric_pae:
				for res1_idx, res2_idx in self.interface_res_pairs[chain_pair]:
					ipae_val = (
						self.pairwise_ipae[chain_pair]["ij"].get((res1_idx, res2_idx), np.nan) +
						self.pairwise_ipae[chain_pair]["ji"].get((res2_idx, res1_idx), np.nan)
					) / 2
					chain_pairwise_assessment_rows.append({
						"Chain Pair": f"{chain1}-{chain2}",
						"Protein Name 1": self.protein_chain_map.get(chain1, "Unknown"),
						"Protein Name 2": self.protein_chain_map.get(chain2, "Unknown"),
						"Residue Pair": f"{res1_idx}-{res2_idx}",
						"iPAE": ipae_val,
						"ipLDDT res1": self.per_chain_iplddt[chain1].get(res1_idx, np.nan),
						"ipLDDT res2": self.per_chain_iplddt[chain2].get(res2_idx, np.nan),
						"Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
						"Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
					})
			else:
				for res1_idx, res2_idx in self.interface_res_pairs[chain_pair]:
					ipae_ij = self.pairwise_ipae[chain_pair]["ij"].get((res1_idx, res2_idx), np.nan)
					ipae_ji = self.pairwise_ipae[chain_pair]["ji"].get((res2_idx, res1_idx), np.nan)
					chain_pairwise_assessment_rows.append({
						"Chain Pair": f"{chain1}-{chain2}",
						"Protein Name 1": self.protein_chain_map.get(chain1, "Unknown"),
						"Protein Name 2": self.protein_chain_map.get(chain2, "Unknown"),
						"Residue Pair": f"{res1_idx}-{res2_idx}",
						"iPAE ij": ipae_ij,
						"iPAE ji": ipae_ji,
						"ipLDDT res1": self.per_chain_iplddt[chain1].get(res1_idx, np.nan),
						"ipLDDT res2": self.per_chain_iplddt[chain2].get(res2_idx, np.nan),
						"Chain Type 1": "IDR" if chain1 in self.idr_chains else "R",
						"Chain Type 2": "IDR" if chain2 in self.idr_chains else "R",
					})

	overall_assessment_keys = {
		"Number of Chains": "num_chains",
		"Number of Interacting Chain Pairs": "num_interacting_chain_pairs",
		"Number of Interface Residues": "num_interface_residues",
		"Number of Contacts": "num_contacts",
		"Average ipLDDT": "avg_iplddt",
		"Average IDR ipLDDT": "avg_idr_iplddt",
		"Average iPAE ij": "avg_ipae_ij",
		"Average iPAE ji": "avg_ipae_ji",
	}

	for col_head, key in overall_assessment_keys.items():
		if self.overall_assessment.get(key, np.nan) is not np.nan:
			overall_assessment_rows.append({
				"Key": col_head,
				"Value": self.overall_assessment.get(key)
			})

	chain_pairwise_assessment_df = pd.DataFrame(chain_pairwise_assessment_rows)
	chainwise_assessment_df = pd.DataFrame(chain_wise_assessment_rows)
	overall_assessment_df = pd.DataFrame(overall_assessment_rows)

	df_dict = {
		"chain_pairwise_assessment": chain_pairwise_assessment_df,
		"chainwise_assessment": chainwise_assessment_df,
		"overall_assessment": overall_assessment_df,
	}

	for k, df_ in df_dict.items():
		df_dict[k] = df_.fillna(np.nan)
		df_dict[k] = df_.map(lambda x: round(x, 2) if isinstance(x, (int, float)) else x)

	with pd.ExcelWriter(self.save_path, engine='openpyxl', mode='w') as writer:
		for sheet_name, df in df_dict.items():

			if df.empty:
				warnings.warn(f"Skipping empty DataFrame for sheet: {sheet_name}")
				continue

			df.to_excel(
				writer,
				sheet_name=sheet_name,
				index=False,
			)
```

### Used in
- [[save_rigid_bodies]]

### Uses


### Tags
#method 