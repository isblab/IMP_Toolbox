```python
def write_to_fasta(
	self,
	fasta_dict: Dict[str, str],
	file_name: str,
	output_dir: str = "./output/af_input",
):
	"""Write the fasta sequences to a file

	Args:
		fasta_dict (dict): dictionary of fasta sequences {header: sequence}
		file_name (str): name of the file
		output_dir (str, optional): Directory to save the fasta file. Defaults to "./output/af_input".
	"""

	os.makedirs(output_dir, exist_ok=True)
	save_path = os.path.join(output_dir, f"{file_name}.fasta")

	with open(save_path, "w") as f:
		for header, sequence in fasta_dict.items():
			f.write(f">{header}\n{sequence}\n")

	print(f"\nFasta file written to {save_path}")
```

### Used in
- [[af_pipeline/AFInput/AlphaFold2/write_job_files|write_job_files]]

### Uses


### Tags
#method