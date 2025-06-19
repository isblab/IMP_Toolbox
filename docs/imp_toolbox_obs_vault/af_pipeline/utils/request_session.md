```python
def request_session(max_retries=3):
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
```