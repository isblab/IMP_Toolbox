```python
def generate_job_entities(
	self,
	job_info: Dict[str, Any],
) -> Tuple[Dict[str, str], str]:
	"""Generate job entities

	job entities are the collection of entities within a job. \n
	Each entity is a proteinChain with a header and sequence. \n

	Args:
		job_info (dict): job information (name, range, count, type)

	Returns:
		Tuple[Dict[str, str], str]: sequences_to_add, job_name
	"""

	# get the job name if provided
	job_name = job_info.get("name", None)

	# get the information for each proteinChain
	headers = self.get_entity_info(job_info, "name", None)
	ranges = self.get_entity_info(job_info, "range", None)
	counts = self.get_entity_info(job_info, "count", 1)

	sequences = self.get_entity_sequences(ranges=ranges, headers=headers)

	job_dict = {
		"job_name": job_name,
		"entities": [],
	}

	for entity_count, (header, sequence, range_, count_) in enumerate(
		zip(headers, sequences, ranges, counts)
	):
		for count_ in range(1, count_ + 1):
			job_dict["entities"].append(
				{
					"header": header,
					"sequence": sequence,
					"range": range_ if range_ else [1, len(sequence)],
					"count": count_,
				}
			)

	# generate job name if not provided
	if not job_name:
		job_name = self.generate_job_name(job_dict)

	# create fasta dictionary for each job {header: sequence}
	sequences_to_add = {}

	for entity in job_dict["entities"]:
		for entity_count in range(1, entity["count"] + 1):
			header = entity["header"]
			sequence = entity["sequence"]
			start, end = entity["range"]

			sequences_to_add[
				f"{header}_{entity_count}_{start}to{end}"
			] = sequence

	# warn if any of the entities is not a proteinChain
	self.warning_not_protien(job_info, job_name)

	return (sequences_to_add, job_name)
```

### Used in 
- [[create_af2_job_cycles]]
- [[create_colabfold_job_cycles]]

### Uses
- [[af_pipeline/AFInput/AlphaFold2/generate_job_name|generate_job_name]]
- [[get_entity_info]]
- [[get_entity_sequences]]
- [[warning_not_protein]]

### Tags
#method