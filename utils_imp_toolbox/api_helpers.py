import requests
import requests.adapters
import json

def request_session(max_retries: int=3) -> requests.sessions.Session:
    """Create a request session with set max retries

    Args:
        max_retries (int, optional): Defaults to 3.

    Returns:
        req_sess (requests.sessions.Session): request session
    """

    req_sess = requests.Session()
    req_sess.mount(
        "https://",
        requests.adapters.HTTPAdapter(max_retries=max_retries)
    )
    req_sess.trust_env = False  # Disable environment variables for proxies
    return req_sess


def request_result(
    get_request: requests.models.Response,
    uniprot_id: str,
    ignore_error: bool=False
) -> dict | bytes | None:
    """Get the result of a `get_request`

    Args:
        get_request (requests.models.Response): `get_request`
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