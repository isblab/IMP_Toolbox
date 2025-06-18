```python
def get_template_settings(self):
	"""Get the template settings for the entity

	- For proteinChain, get the template settings from the entity_info dictionary
	- For dnaSequence, rnaSequence, ligand and ion, return an empty dictionary
	- For proteinChain, if useStructureTemplate is not provided, set it to True
	- For proteinChain, if maxTemplateDate is not provided, set it to 2021-09-30
	- For proteinChain, if useStructureTemplate is False, ignore maxTemplateDate
		and raise a warning if maxTemplateDate is provided

	Returns:
		dict: template settings for the entity:
			- maxTemplateDate
			- useStructureTemplate
	"""

	template_dict = {}

	if self.entity_type == "proteinChain":
		if self.entity_info.get("useStructureTemplate", True):
			template_dict = {
				"maxTemplateDate": self.entity_info.get(
					"maxTemplateDate", "2021-09-30"
				),
				"useStructureTemplate": True
			}
		else:
			if "maxTemplateDate" in self.entity_info:
				warnings.warn(
					f"maxTemplateDate is provided for {self.entity_name} \
					but useStructureTemplate is False. \
					Ignoring maxTemplateDate."
				)
			template_dict = {
				"useStructureTemplate": False
			}

	return template_dict
```

### Used in 
- [[fill_up_entity]]

### Uses


### Tags
#method