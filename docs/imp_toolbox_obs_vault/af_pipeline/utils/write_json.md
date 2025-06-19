```python
def write_json(file_path, data):
    """Write data to a json file

    Args:
        file_path (str): path to json file
        data (dict): data to write
    """

    with open(file_path, "w") as f:
        json.dump(data, f)
```