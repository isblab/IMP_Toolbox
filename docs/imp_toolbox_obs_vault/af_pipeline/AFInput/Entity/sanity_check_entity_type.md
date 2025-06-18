```python
@staticmethod
def sanity_check_entity_type(entity_type):
	"""Sanity check the entity
	allowed entity types: proteinChain, dnaSequence, rnaSequence, ligand, ion
	"""

	if entity_type not in ENTITY_TYPES:
		raise Exception(f"Invalid entity type {entity_type}")
```

### Used in 
- [[Entity]]

### Uses


### Tags
#staticmethod