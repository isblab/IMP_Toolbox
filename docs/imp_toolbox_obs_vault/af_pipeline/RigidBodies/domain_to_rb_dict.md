```python
def domain_to_rb_dict(self, domain: list):
	"""Convert the domain list to a dictionary of rigid bodies.
	- The rigid bodies are represented as a dictionary with chain_id as the key and
		a list of residue numbers as the value.

	Args:
		domain (list): list of residue indices in the domain

	Returns:
		rb_dict (dict): pseudo-rigid body in the form of a dictionary

	Example:
		if predicted structure has chains: A (20 aa), B (30 aa), C (50 aa) \n
		such that, actual residue numbers are
		A: [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39] \n
		B: [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49] \n
		C: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49] \n
		and detected domain is [0, 1, 2, 3, 4, 5, 20, 21, 22, 23, 54, 55, 56, 57, 58] \n
		rb_dict = {
			'A': [20, 21, 22, 23, 24, 25],
			'B': [20, 21, 22, 23, 24],
			'C': [5, 6, 7, 8, 9]
		}
	"""

	rb_dict = defaultdict(list)

	for res_idx in domain:

		res_num = self.idx_to_num[res_idx].get("res_num")
		chain_id = self.idx_to_num[res_idx].get("chain_id")

		rb_dict[chain_id].append(res_num)

	return rb_dict
```

### Used in
- [[predict_domains]]

### Uses


### Tags
#method 