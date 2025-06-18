```python
def fill_up_entity(self):
	"""Fill up the entity with the required information"""

	self.entity_count = self.get_entity_count()
	self.real_sequence = self.get_real_sequence()
	self.start, self.end = self.get_entity_range()
	self.glycans = self.get_glycans()
	self.modifications = self.get_modifications()
	self.template_settings = self.get_template_settings()
```

### Used in
- [[Entity]]

### Uses
- [[get_template_settings]]
- [[get_modifications]]
- [[get_glycans]]
- [[get_entity_range]]
- [[get_real_sequence]]
- [[get_entity_count]]

### Tags
#method