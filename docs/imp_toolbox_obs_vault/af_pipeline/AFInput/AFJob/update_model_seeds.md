```python
def update_model_seeds(self):
	"""Update the model seeds
	- If modelSeeds is an integer, generate that many seeds
	- If modelSeeds is a list, use those seeds
	- If modelSeeds is not provided, empty list (auto seed by AF3)

	Raises:
		Exception: modelSeeds must be an integer or a list
	"""

	model_seeds = self.job_info.get("modelSeeds")

	if "modelSeeds" in self.job_info:

		if isinstance(model_seeds, int):
			self.model_seeds = self.generate_seeds(num_seeds=model_seeds)

		elif isinstance(model_seeds, list):
			self.model_seeds = model_seeds

		else:
			raise Exception("modelSeeds must be an integer or a list")
```

### Used in 
- [[create_job]]

### Uses
- [[generate_seeds]]

### Tags
#method