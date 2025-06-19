```python
def update_config(
    input_file: str,
    updates: dict = None,
    mode: str = "replace",
):
    """ Update config file with a new field or update an existing field

    Args:
        input_file (str): Path to input config file
        updates (dict, optional): Fields to update in the config file. Defaults to None.
        mode (str, optional): Mode to update the config file. Defaults to "replace". ("append" or "replace")
    """

    yaml = YAML()

    update_fields = list(updates.keys()) if updates else []

    if len(update_fields) == 0:

        print("No fields to update in config")
        return None

    yaml.preserve_quotes = True

    with open(input_file, "r") as f:
        config_yaml = yaml.load(f)

    existing_fields = list(config_yaml.keys())

    for field in update_fields:

        add_field = False

        if field in existing_fields:
            if mode == "replace":
                config_yaml[field] = updates[field]
            elif mode == "append":
                # need to change this, not working as expected
                config_yaml[field].update(updates[field])
            else:
                raise ValueError("Invalid mode. Use 'replace' or 'append")

        else:
            print(f"{field} not found in config")
            print("Adding field to config")
            add_field = True

        if add_field:
            config_yaml[field] = updates[field]
            add_field = False

    with open(input_file, "w") as f:
        yaml.dump(config_yaml, f)

    print(f"Config file updated with {update_fields}")
```