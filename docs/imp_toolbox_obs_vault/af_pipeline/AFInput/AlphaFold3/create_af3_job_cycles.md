```python
def create_af3_job_cycles(self) -> Dict[str, List[Dict[str, Any]]]:
	"""Create job cycles for AlphaFold3

	each job cycle is a list of jobs with each job being a dictionary
	see :py:mod:`AFCycle.seed_jobs` or `AFJob.create_job` for the job 
	dictionary format

	Returns:
		job_cycles (dict): dictionary of job cycles {job_cycle: job_list}
	"""

	job_cycles = {}

	for job_cycle, jobs_info in self.input_yml.items():

		print("Creating job cycle", job_cycle, "\n")

		af_cycle = AFCycle(
			jobs_info=jobs_info,
			protein_sequences=self.protein_sequences,
			nucleic_acid_sequences=self.nucleic_acid_sequences,
			entities_map=self.entities_map,
		)
		af_cycle.update_cycle()
		job_cycles[job_cycle] = af_cycle.job_list

	return job_cycles
```

### Used in
- 

### Uses
- [[AFCycle]]
- [[update_cycle]]

### Tags
#user_function #method 