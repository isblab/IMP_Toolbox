```python
def get_name_fragment(self)-> str:
	"""Get the name fragments of the entity"""

	return f"{self.name}_{self.count}_{self.start}to{self.end}"
```

### Used in
- [[update_af_sequences]]

### Uses


### Tags
#method