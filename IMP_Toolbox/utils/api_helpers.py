import requests
import requests.adapters
import json

def request_session(max_retries: int=3) -> requests.sessions.Session:
    """ Create a request session with set max retries

    ## Arguments:

    - **max_retries (int, optional):**:<br />
        Number of times to retry a request in case of failure. Defaults to 3.

    ## Returns:

    - **requests.sessions.Session**:<br />
        A request session with the specified max retries.
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
    identifier: str | None=None,
    ignore_error: bool=False
) -> dict | bytes | None:
    """ Get the result of a `get_request` and return it as a dict, bytes or None if an error occurred.

    ## Arguments:

    - **get_request (requests.models.Response)**:<br />
        The response object returned by a requests.get() call.

    - **identifier (str | None, optional):**:<br />
        The identifier of the request, used for error reporting.

    - **ignore_error (bool, optional):**:<br />
        Whether to ignore errors and not print error messages. Defaults to False.

    ## Returns:

    - **dict | bytes | None**:<br />
        The result of the request as a dict if successful, bytes if not JSON, or None if an error occurred.
    """

    if get_request.status_code == 200:
        try:
            return get_request.json()
        except json.decoder.JSONDecodeError:
            return get_request.content
    else:
        print(f"Error while requesting {identifier}") if not ignore_error else None

        return None