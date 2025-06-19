```python
def read_json(file_path):
    """Load a json file

    Args:
        file_path (str): path to json file

    Returns:
        data (dict): data from json file
    """

    with open(file_path, "r") as f:
        data = json.load(f)
    return data
```