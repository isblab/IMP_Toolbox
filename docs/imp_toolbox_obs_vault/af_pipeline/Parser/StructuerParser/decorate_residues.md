```python
def decorate_residue(self, residue: Bio.PDB.Residue.Residue):

	symbol = residue.get_resname()

	if symbol in PROTEIN_ENTITIES:
		residue.xtra["entityType"] = "proteinChain"

	elif symbol in DNA_ENTITIES:
		residue.xtra["entityType"] = "dnaSequence"

	elif symbol in RNA_ENTITIES:
		residue.xtra["entityType"] = "rnaSequence"

	elif symbol in ALLOWED_LIGANDS:
		residue.xtra["entitiyType"] = "ligand"

	elif symbol in ION:
		residue.xtra["entityType"] = "ion"

	else:
		warnings.warn(
			f"""
			The residue {symbol} does not belong to any known entity types.
			Setting 'entityType' to None.
			"""
		)
		residue.xtra["entityType"] = None
```


### Used in 
- [[get_structure]]

### Uses


### Tags
#method 