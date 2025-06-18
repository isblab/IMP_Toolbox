```python
def generate_job_name(self):
	"""Generate a job name"""

	job_name = "_".join(self.name_fragments)
	self.job_name = job_name
```

### Used in
- [[create_job]]

### Uses


### Tags
#method