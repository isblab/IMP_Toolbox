```python
def generate_seeds(self, num_seeds: int) -> List[int]:
	"""Generate model seeds"""

	model_seeds = random.sample(range(1, SEED_MULTIPLIER * num_seeds), num_seeds)

	return model_seeds
```

### Used in
- [[update_model_seeds]]

### Uses

### Tags
#method