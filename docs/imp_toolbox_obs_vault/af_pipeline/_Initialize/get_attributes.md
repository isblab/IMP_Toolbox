```python
def get_attributes(self):
	"""
	Extract the following from the input files:
		1. Residue positions of all residues for each chain.
		2. Ca or representative atom coordinates.
		3. Ca or representative atom pLDDT.
		4. Tokenized chain IDs.
		5. Tokenized residue IDs.
		6. Chain lengths.
		7. PAE matrix.
		8. Average PAE matrix.
		9. Min PAE for each residue.
	"""

	data = self.dataparser.get_data_dict()
	self.pae = self.dataparser.get_pae(data=data)
	self.avg_pae = self.dataparser.get_avg_pae(pae=self.pae)
	self.contact_probs_mat = self.dataparser.get_contact_probs_mat(data=data)

	if self.contact_probs_mat is not None:
		self.avg_contact_probs_mat = self.dataparser.get_avg_contact_probs_mat(
			contact_probs_mat=self.contact_probs_mat
		)

	self.token_chain_ids = self.dataparser.get_token_chain_ids(data=data)
	self.token_res_ids = self.dataparser.get_token_res_ids(data=data)

	if self.token_chain_ids is None or self.token_res_ids is None:
		self.token_chain_ids, self.token_res_ids = self.structureparser.get_token_chain_res_ids()

	# Ca-coords of all residues for each chain.
	self.coords_list = self.structureparser.get_ca_coordinates()
	# Ca-plddt of all residues for each chain.
	self.plddt_list = self.structureparser.get_ca_plddt()

	self.pae = self.update_pae(
		pae=self.pae,
		token_res_ids=self.token_res_ids,
		token_chain_ids=self.token_chain_ids,
		average_atom_pae=self.average_atom_pae,
	)
	self.avg_pae = self.update_pae(
		pae=self.avg_pae,
		token_res_ids=self.token_res_ids,
		token_chain_ids=self.token_chain_ids,
		average_atom_pae=self.average_atom_pae,
	)

	if self.contact_probs_mat is not None:

		self.contact_probs_mat = self.update_contact_probs(
			contact_probs_mat=self.contact_probs_mat,
			token_chain_ids=self.token_chain_ids,
			token_res_ids=self.token_res_ids,
			average_atom_pae=self.average_atom_pae,
		)
		self.avg_contact_probs_mat = self.update_contact_probs(
			contact_probs_mat=self.avg_contact_probs_mat,
			token_chain_ids=self.token_chain_ids,
			token_res_ids=self.token_res_ids,
			average_atom_pae=self.average_atom_pae,
		)

	self.token_chain_ids, self.token_res_ids = self.update_token_ids(
		token_chain_ids=self.token_chain_ids,
		token_res_ids=self.token_res_ids,
		average_atom_pae=self.average_atom_pae,
	)

	self.lengths_dict = self.get_chain_lengths(
		token_chain_ids=self.token_chain_ids,
	)
	self.sanity_check()
	self.idx_to_num, self.num_to_idx = self.renumber.residue_map(
		token_chain_ids=self.token_chain_ids,
		token_res_ids=self.token_res_ids,
	)
```

### Used in
- [[_Initialize]]

### Uses
- [[get_data_dict]]
- [[get_pae]]
- [[get_contact_probs_mat]]
- [[get_token_chain_ids]]
- [[get_token_res_ids]]
- [[get_token_chain_res_ids]]
- [[get_cb_coordinates]]
- [[get_cb_plddt]]
- [[update_pae]]
- [[update_contact_probs]]
- [[update_token_ids]]
- [[get_chain_lengths]]
- [[sanity_check]]
- [[residue_map]]

### Tags
#method 