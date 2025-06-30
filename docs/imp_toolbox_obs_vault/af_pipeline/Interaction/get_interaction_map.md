```python
@staticmethod
def get_interaction_map(
	coords1: np.ndarray,
	coords2: np.ndarray,
	contact_threshold: float,
	map_type: str
	):
	"""
	Create an interaction map, given the input coordinates.

	Returns a distance map or a contact map, based on the map_type specified.
	"""

	distance_map = Interaction.get_distance_map(
		coords1 = coords1,
		coords2 = coords2
	)

	if map_type == "distance":
		return distance_map

	elif map_type == "contact":
		contact_map = Interaction.get_contact_map(
			distance_map = distance_map,
			contact_threshold = contact_threshold
		)
		return contact_map

	else:
		raise Exception("Invalid map_type specified...")
```

