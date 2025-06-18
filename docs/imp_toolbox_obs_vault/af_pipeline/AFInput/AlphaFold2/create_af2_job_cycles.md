```python
def create_af2_job_cycles(self) -> Dict[str, List[Tuple[Dict[str, str], str]]]:
	"""Create job cycles for AlphaFold2

	each job cycle is a list of jobs. \n
	each job is a tuple of `sequences_to_add` and `job_name`. \n
	`sequences_to_add` is a dictionary of fasta sequences {header: sequence} \n

	Returns:
		job_cycles (dict): dictionary of job cycles {(}job_cycle: job_list}
	"""

	job_cycles = {}

	for job_cycle, jobs_info in self.input_yml.items():

		job_list = []

		for job_info in jobs_info:
			sequences_to_add, job_name = self.generate_job_entities(
				job_info=job_info
			)
			job_list.append((sequences_to_add, job_name))

		job_cycles[job_cycle] = job_list

	return job_cycles
```

### Used in


### Uses
- [[generate_job_entities]]

### Tags
#user_function 

### Tags
#method