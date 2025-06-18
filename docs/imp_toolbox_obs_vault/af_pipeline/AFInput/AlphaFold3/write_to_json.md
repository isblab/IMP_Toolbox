```python
def write_to_json(
	self,
	sets_of_n_jobs: List[List[Dict[str, Any]]],
	file_name: str,
	output_dir: str = "./output/af_input",
):
	"""Write the sets of n jobs to json files

	Args:
		sets_of_n_jobs (list): list of lists, each list containing n jobs
		file_name (str): name of the file
		output_dir (str, optional): Directory to save the job_files.
			Defaults to "./output/af_input".
	"""

	os.makedirs(output_dir, exist_ok=True)
	for i, job_set in enumerate(sets_of_n_jobs):

		save_path = os.path.join(output_dir, f"{file_name}_set_{i}.json")

		with open(save_path, "w") as f:
			json.dump(job_set, f, indent=4)

		print(f"{len(job_set)} jobs written for {file_name}_set_{i}")
```

### Used in
- [[af_pipeline/AFInput/AlphaFold3/write_job_files|write_job_files]]

### Uses

### Tags
#method