```python
@staticmethod
def sanity_check_small_molecule(entity_type, entity_name):
	"""Sanity check the small molecules"""

	if (entity_type == "ligand" and entity_name not in LIGAND) or (
		entity_type == "ion" and entity_name not in ION
	):
		raise Exception(f"Invalid small molecule {entity_name}")
```

### Used in
- [[Entity]]

### Uses


### Tags
#staticmethod 