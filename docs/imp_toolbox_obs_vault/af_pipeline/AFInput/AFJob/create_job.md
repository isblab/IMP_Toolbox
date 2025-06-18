```python
def create_job(self) -> Dict[str, Any]:
	"""Create a job from the job info

	Returns:
		job_dict (dict): job dictionary in the following format:
			{
				"name": "job_name",
				"modelSeeds": [1, 2],
				"sequences": [... ]
			}
	"""

	self.update_job_name()
	self.update_model_seeds()
	self.update_af_sequences()

	if self.job_name is None:
		self.generate_job_name()

	job_dict = {
		"name": self.job_name,
		"modelSeeds": self.model_seeds,
		"sequences": self.af_sequences,
	}

	return job_dict
```


### Used in
- [[update_cycle]]

### Uses
- [[update_job_name]]
- [[update_model_seeds]]
- [[update_af_sequences]]
- [[af_pipeline/AFInput/AFJob/generate_job_name|generate_job_name]]

### Tags
#method