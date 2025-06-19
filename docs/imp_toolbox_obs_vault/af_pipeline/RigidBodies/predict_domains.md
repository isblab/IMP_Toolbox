```python
def predict_domains(
	self,
	num_res: int = 1,
	num_proteins: int = 1,
	plddt_filter: bool = True
):
	"""Predict domains from a PAE file.
	- Three implementations are available:
		1. igraph based
		2. networkx based
		3. label_propagation based

	(1) is significantly faster than (2)

	Args:
		num_res (int): Minimum number of residues in a rigid body
		num_proteins (int): Minimum number of proteins in a rigid body
		plddt_filter (bool): Filter the residues based on the pLDDT cutoff

	Raises:
		ValueError: Invalid library specified. Use 'igraph' or 'networkx'

	Returns:
		domains (list): List of domains in which each domain is a rigid body dictionary

	A rigid body dictionary is of the form:
	- {
		ch1: [res_num, ...],
		ch2: [res_num, ...],
		...
	}
	"""

	print("Predicting domains...")
	start_time = time.time()

	pae_matrix = self.pae

	if self.library == "igraph":
		f = domains_from_pae_matrix_igraph

	elif self.library == "networkx":
		f = domains_from_pae_matrix_networkx

	elif self.library == "label_propagation":
		f = domains_from_pae_matrix_label_propagation

	else:
		raise ValueError("Invalid library specified. Use 'igraph' or 'networkx")

	if f == domains_from_pae_matrix_igraph or f == domains_from_pae_matrix_networkx:
		domains = f(
			pae_matrix,
			pae_power=self.pae_power,
			pae_cutoff=self.pae_cutoff,
			graph_resolution=self.resolution,
		)
	elif f == domains_from_pae_matrix_label_propagation:
		domains = f(
			pae_matrix,
			pae_power=self.pae_power,
			pae_cutoff=self.pae_cutoff,
			random_seed=self.random_seed,
		)

	# domains is a list of frozensets or lists
	# each frozenset/list contains residue indices for residues in a domain
	for idx, domain in enumerate(domains):

		if isinstance(domain, frozenset):
			domain = list(domain)

		# rb_dict is a dictionary of rigid bodies
		# each rigid body is represented as a dictionary with chain_id as the key and a list of residue numbers as the value
		rb_dict = self.domain_to_rb_dict(domain=domain)

		# removing residues with pLDDT score below the cutoff
		if plddt_filter:
			rb_dict = self.filter_plddt(
				rb_dict=rb_dict,
				patch_threshold=self.patch_threshold,
			)

		domains[idx] = rb_dict

	# Remove domains with number of proteins less than `num_proteins`
	domains = [
		rb_dict
		for rb_dict in domains
		if len(rb_dict) >= num_proteins
	]

	# Remove domains with number of residues less than `num_res`
	domains = [
		rb_dict
		for rb_dict in domains
		if sum([len(res_list) for res_list in rb_dict.values()]) >= num_res
	]

	end_time = time.time()
	print(f"Done predicting pseudo-rigid domains in {end_time - start_time:.2f} seconds")

	return domains
```

### Used in


### Uses
- [[domains_from_pae_matrix_igraph]]
- [[domains_from_pae_matrix_networkx]]
- [[domains_from_pae_matrix_label_propagation]]
- [[domain_to_rb_dict]]
- [[filter_plddt]]

### Tags
#method #user_function 