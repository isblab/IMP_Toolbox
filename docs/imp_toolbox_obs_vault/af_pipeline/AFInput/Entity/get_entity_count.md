```python
def get_entity_count(self)-> int:
	"""Get the count of the entity

	Returns:
		int: count or copy number of the entity (default: 1)
	"""

	entity_count = self.entity_info.get("count", 1)

	return entity_count
```

### Used in
- [[fill_up_entity]]

### Uses


### Tags
#method