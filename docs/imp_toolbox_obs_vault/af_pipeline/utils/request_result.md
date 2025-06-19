```python
def request_result(
    get_request: requests.models.Response,
    uniprot_id: str,
    ignore_error: bool=False
):
    """Get the result of a get request

    Args:
        get_request (requests.models.Response): get request
        uniprot_id (str): valid entrypoint id
        ignore_error (bool, optional): Defaults to False.

    """

    if get_request.status_code == 200:
        try:
            return get_request.json()
        except json.decoder.JSONDecodeError:
            return get_request.content
    else:
        print(f"Error while requesting {uniprot_id}") \
            if not ignore_error else None

        return None
```