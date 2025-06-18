```python
@staticmethod
def warning_not_protien(
	job_info: Dict[str, Any],
	job_name: str
):
	"""Warn if entity is not a protein

	AF2/ ColabFold only supports proteinChain entities. \n
	Will skip the entities which are not proteins. \n

	Args:
		job_info (dict): job information
		job_name (str): job name
	"""

	if any(
		[
			entity_type != "proteinChain"
			for entity_type in [
				entity["type"] for entity in job_info["entities"]
			]
		]
	):
		warnings.warn(
			f"""
			AF2/ ColabFold only supports proteinChain entities.
			Will skip the entities which are not proteins.
			{job_name} will be created with only proteinChain entities.
			"""
		)
```

### Used in
- [[generate_job_entities]]

### Uses


### Tags
#method #staticmethod 