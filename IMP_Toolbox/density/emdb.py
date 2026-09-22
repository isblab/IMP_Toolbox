from IMP_Toolbox.utils.api_helpers import request_session
from IMP_Toolbox.constants.imp_toolbox_constants import (
    APIurl,
    MAX_API_RETRIES,
)


def fetch_emdb_map(emdb_id: str, max_retries: int = MAX_API_RETRIES) -> bytes:
    """ Fetch density map from EMDB

    ## Arguments:

    - **emdb_id (str)**:<br />
        EMDB ID for which the density map is to be fetched.
        This should be in the format "EMD-XXXX" where XXXX is a 4 digit number.

    - **max_retries (int, optional):**:<br />
        Maximum number of retries for the API request. Default is 3.

    ## Returns:

    - **bytes**:<br />
        Density map file content in bytes.
    """

    EMDB_MAP_URL = APIurl.emdb_ftp_map.substitute(
        emdb_id_hyphen=emdb_id,
        emdb_id_underscore=emdb_id.lower().replace("-", "_")
    )

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(EMDB_MAP_URL)

    if response.status_code == 200:
        print("Successfully fetched EMDB map for given EMDB id")
        return response.content

    else:
        raise Exception("Error while requesting EMDB map for given EMDB id")

def fetch_emdb_mask(
    emdb_id: str,
    mask_name: str,
    max_retries: int = MAX_API_RETRIES,
) -> bytes:
    """ Fetch mask from EMDB

    ## Arguments:

    - **emdb_id (str)**:<br />
        EMDB ID for which the mask is to be fetched.
        This should be in the format "EMD-XXXX" where XXXX is a 4 digit number.

    - **mask_name (str)**:<br />
        Name of the mask to be fetched. This should be in the format "mask_XX".

    - **max_retries (int, optional):**:<br />
        Maximum number of retries for the API request. Default is 3.

    ## Returns:

    - **bytes**:<br />
        Mask file content in bytes.
    """

    if mask_name is None:
        raise ValueError("Mask name must be provided to fetch EMDB mask")

    EMDB_MASK_URL = APIurl.emdb_ftp_mask.substitute(
        emdb_id_hyphen=emdb_id,
        mask_name=mask_name
    )

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(EMDB_MASK_URL)

    if response.status_code == 200:
        print("Successfully fetched EMDB mask for given EMDB id and mask name")
        return response.content

    else:
        raise Exception(
            "Error while requesting EMDB mask for given EMDB id and mask name"
        )