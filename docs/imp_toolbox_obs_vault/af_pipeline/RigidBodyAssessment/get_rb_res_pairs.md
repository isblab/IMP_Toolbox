```python
def get_rb_res_pairs(self):
	"""Get all unique residue pairs in the rigid body.

	Returns:
		rb_res_pairs (defaultdict): A dictionary where keys are chain pairs (tuples) and values are lists of residue index pairs.
		Each residue index pair is a tuple of indices from the two chains in the rigid body.
	"""

	rb_res_pairs = defaultdict(list)

	for chain_pair in self.chain_pairs:

		chain1, chain2 = chain_pair

		res1_list = self.rb_dict[chain1]
		res2_list = self.rb_dict[chain2]

		res1_idxs = [self.num_to_idx[chain1][res_num] for res_num in res1_list]
		res2_idxs = [self.num_to_idx[chain2][res_num] for res_num in res2_list]

		# Create pairs of residues from the two chains
		pairs = list(product(res1_idxs, res2_idxs))
		rb_res_pairs[chain_pair].extend(pairs)

	return rb_res_pairs
```

### Used in
- [[RigidBodyAssessment]]

### Uses


### Tags
#method 