```python
def update_cycle(self):
	"""Update the cycle with the jobs

	For each job in jobs_info, creates an AFJob instance and uses it to create
	a job dictionary. The job dictionary is then seeded to create multiple jobs
	based on model seeds.
	"""

	for job_info in self.jobs_info:
		af_job = AFJob(
			job_info=job_info,
			protein_sequences=self.protein_sequences,
			nucleic_acid_sequences=self.nucleic_acid_sequences,
			entities_map=self.entities_map,
		)

		job_dict = af_job.create_job()
		self.seed_jobs(job_dict)
```

- Create a AF3-compatible job dictionary for each job given `self.jobs_info`.
### Used in
- [[create_af3_job_cycles]]

### Uses
- [[seed_jobs]]
- [[create_job]]
- [[AFJob]]

### Tags
#method