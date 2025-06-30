```python
@staticmethod
def extract_perresidue_quantity(
	residue: Bio.PDB.Residue.Residue,
	quantity: str
):
	"""Extract per-residue quantities from the Biopython residue object.

	Given the Biopython residue object, return the specified quantity: \n
		1. residue or nucleotide or ion position \n
		2. Cb-coordinate or representative atom coordinate \n
		3. Cb-pLDDT or representative atom pLDDT

	Args:

		residue (Bio.PDB.Residue.Residue):
			Biopython residue object.

		quantity (str):
			Quantity to extract.

	Returns:

		extracted_quantity (int | np.ndarray | float):
			The extracted quantity based on the specified `quantity`.
			- If `quantity` is "res_pos", returns residue position (int).
			- If `quantity` is "coords", returns coordinates of the 
				representative atom (np.ndarray).
			- If `quantity` is "plddt", returns pLDDT value of the 
				representative atom (float).
	"""

	# Using representative atoms as specified by AF3.
	# https://github.com/google-deepmind/alphafold3/blob/main/src/alphafold3/model/features.py#L1317

	symbol = residue.get_resname()
	rep_atom = residue.child_list[0].get_name()

	if residue.xtra.get("entityType") == "proteinChain":

		if (
			"CB" in residue.child_dict and symbol in PROTEIN_ENTITIES
		):  # this includes modifications
			rep_atom = "CB"

		elif (
			"CB" not in residue.child_dict and symbol in ONLY_CA_RESIDUES
		):  # this includes modifications
			rep_atom = "CA"

		else:
			raise Exception(
				f"""
				Are you sure this is a protein chain?
				residue {symbol} in chain {residue.parent.id}
				does not have a Cb-atom or a Ca-atom.
				"""
			)

	elif residue.xtra.get("entityType") in ["dnaSequence", "rnaSequence"]:

		if symbol in PURINES:  # this includes modifications
			rep_atom = "C4"

		elif symbol in PYRIMIDINES:  # this includes modifications
			rep_atom = "C2"

	elif residue.xtra.get("entityType") == "ion" and symbol in ION:
		rep_atom = symbol

	elif residue.xtra.get("entityType") == "ligand":
		rep_atom = residue.child_list[0].get_name()
		warnings.warn(
			f"""
			Can not determine representative atom for ligand {symbol}
			in chain {residue.parent.get_id()}
			Setting representative atom to {rep_atom}.
			"""
		)

	else:
		rep_atom = residue.child_list[0].get_name()
		warnings.warn(
			f"""
			Unknown entity type for residue {symbol}
			in chain {residue.parent.id}.
			It could be a glycan modification.
			Setting representative atom to {rep_atom}.
			"""
		)

	if quantity == "res_pos":
		return residue.id[1]

	elif quantity == "coords":
		coords = residue[rep_atom].coord
		return coords

	elif quantity == "plddt":
		plddt = residue[rep_atom].bfactor
		return plddt

	else:
		raise Exception(
			f"Specified quantity: {quantity} does not exist for {symbol}"
		)

```

### Used in
- [[get_token_chain_res_ids]]
- [[get_cb_coordinates]]
- [[get_cb_plddt]]

### Uses


### Tags
#method 