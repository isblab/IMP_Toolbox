```python
def generate_job_name(
	self,
	job_dict: Dict[str, Any],
) -> str:
	"""Generate job name (if not provided)

	see :py:mod:`AFInput.generate_job_entities` for the job dictionary format.

	Args:
		job_dict (dict): job dictionary

	Returns:
		job_name (str): job name
	"""

	job_name = ""

	fragments = defaultdict(list)

	for entity in job_dict["entities"]:
		header = entity["header"]
		start, end = entity["range"]
		count = entity["count"]

		fragments[f"{header}_{start}to{end}"].append(count)

	fragments = {k: max(v) for k, v in fragments.items()}

	for header, count in fragments.items():
		header_, range_ = header.split("_")
		job_name += f"{header_}_{count}_{range_}_"

	job_name = job_name[:-1] if job_name[-1] == "_" else job_name

	return job_name
```

### Used in
- [[generate_job_entities]]

### Uses


### Tags
#method