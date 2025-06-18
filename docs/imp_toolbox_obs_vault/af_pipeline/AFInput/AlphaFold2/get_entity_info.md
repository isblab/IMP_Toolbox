```python
def get_entity_info(
	self,
	job_info: Dict[str, Any],
	info_type: str,
	default_val: Any
) -> List[Dict[str, Any]]:
	"""Get the entity information

	Get the required information for each entity in the job

	Args:
		job_info (dict): job information (name, range, count, type)
		info_type (str): type of information to get (name, range, count, type)
		default_val (Any): default value if not found

	Returns:
		List[Dict[str, Any]]: list of entity information for the given type
	"""

	return [
		entity.get(info_type, default_val)
		for entity in job_info["entities"]
		if entity["type"] == "proteinChain"
	]
```

### Used in
- [[generate_job_entities]]

### Uses


### Tags
#method