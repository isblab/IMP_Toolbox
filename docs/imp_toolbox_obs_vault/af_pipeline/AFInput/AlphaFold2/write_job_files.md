```python
def write_job_files(
	self,
	job_cycles: Dict[str, List[Tuple[Dict[str, str], str]]],
	output_dir: str = "./output/af_input",
):
	"""Write job files to the output directory

	Args:
		job_cycles (dict): dictionary of job cycles {job_cycle: job_list}
		output_dir (str, optional): Defaults to "./output/af_input".
	"""

	for job_cycle, job_list in job_cycles.items():

		os.makedirs(os.path.join(output_dir, job_cycle), exist_ok=True)

		for fasta_dict, job_name in job_list:

			self.write_to_fasta(
				fasta_dict=fasta_dict,
				file_name=job_name,
				output_dir=os.path.join(output_dir, job_cycle),
			)

	print("\nAll job files written to", output_dir)
```

### Used in


### Uses
- [[write_to_fasta]]

### Tags
#user_function #method 