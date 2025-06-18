```python
def write_job_files(
	self,
	job_cycles: Dict[str, List[Dict[str, Any]]],
	output_dir: str = "./output/af_input",
	num_jobs_per_file: int = 20,
):
	"""Write job files to the output directory

	Args:
		job_cycles (dict): dictionary of job cycles {job_cycle: job_list}
		output_dir (str, optional): Directory to save the job_files. Defaults to "./output/af_input".
		num_jobs_per_file (int, optional): Number of jobs per file. Defaults to 20.
	"""

	assert 101 > num_jobs_per_file > 0; "Number of jobs per file must be within 1 and 100"

	for job_cycle, job_list in job_cycles.items():

		sets_of_n_jobs = [
			job_list[i : i + num_jobs_per_file]
			for i in range(0, len(job_list), num_jobs_per_file)
		]
		os.makedirs(output_dir, exist_ok=True)

		self.write_to_json(
			sets_of_n_jobs=sets_of_n_jobs,
			file_name=job_cycle,
			output_dir=os.path.join(output_dir, job_cycle),
		)
```

### Used in


### Uses
- [[write_to_json]]

### Tags
#user_function #method 