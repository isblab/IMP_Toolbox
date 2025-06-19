```python
def get_overall_assessment(self):
	""" Get overall assessment of the rigid body.

	Returns:
		overall_assessment (dict): A dictionary containing overall statistics about the rigid body.
		It includes the number of chains, number of interacting chain pairs, number of interface residues,
		number of contacts, average ipLDDT, average IDR ipLDDT, average iPAE ij, and average iPAE ji.
	"""

	overall_assessment = {}

	overall_assessment["num_chains"] = len(self.unique_chains)

	overall_assessment["num_interacting_chain_pairs"] = len(self.interface_res_pairs)

	overall_assessment["num_interface_residues"] = sum(
		len(res_list)
		for res_list in self.per_chain_interface_residues.values()
	)

	overall_assessment["num_contacts"] = sum(
		len(contact_pairs)
		for contact_pairs in self.interface_res_pairs.values()
	)

	global_iplddt_scores = [
		iplddt
		for iplddt_scores in self.per_chain_iplddt.values()
		for iplddt in iplddt_scores.values()
	]

	overall_assessment["avg_iplddt"] = (
		np.mean(global_iplddt_scores) if global_iplddt_scores else np.nan
	)

	global_idr_iplddt_scores = [
		iplddt
		for chain_id, iplddt_scores in self.per_chain_iplddt.items()
		for iplddt in iplddt_scores.values()
		if chain_id in self.idr_chains
	]

	overall_assessment["avg_idr_iplddt"] = (
		np.mean(global_idr_iplddt_scores) if global_idr_iplddt_scores else np.nan
	)

	global_ipae_ij_scores = [
		ipae
		for ipae_dict in self.pairwise_ipae.values()
		for ipae in ipae_dict["ij"].values()
	]

	global_ipae_ji_scores = [
		ipae
		for ipae_dict in self.pairwise_ipae.values()
		for ipae in ipae_dict["ji"].values()
	]

	overall_assessment["avg_ipae_ij"] = (
		np.mean(global_ipae_ij_scores) if global_ipae_ij_scores else np.nan
	)

	overall_assessment["avg_ipae_ji"] = (
		np.mean(global_ipae_ji_scores) if global_ipae_ji_scores else np.nan
	)

	return overall_assessment
```


### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 