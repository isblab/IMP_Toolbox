```python
def seed_jobs(
	self,
	job_dict: Dict[str, Any],
):
	"""Create a job for each model seed

	Args:
		job_dict (dict): job dictionary in the following format:
			{
				"name": "job_name",
				"modelSeeds": [1, 2],
				"sequences": [... ]
			}

	will lead to -->
		{
			"name": "job_name",
			"modelSeeds": [1],
			"sequences": [... ]
		},
		{
			"name": "job_name",
			"modelSeeds": [2],
			"sequences": [... ]
		}
	"""

	if len(job_dict["modelSeeds"]) == 0:
		self.job_list.append(job_dict)

	else:
		for seed in job_dict["modelSeeds"]:
			job_copy = job_dict.copy()
			job_copy["modelSeeds"] = [seed]
			job_copy["name"] = f"{job_dict['name']}_{seed}"
			self.job_list.append(job_copy)
```

- Split the jobs in case of multiple seeds and add each job to `self.jobs_list`.

### Used in
- [[update_cycle]]

### Uses

### Tags
#method