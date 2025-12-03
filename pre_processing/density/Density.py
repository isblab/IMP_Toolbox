from IMP_Toolbox.utils import request_session

class EMDBMaps:
    """Fetch density maps and masks from EMDB"""

    def __init__(
        self,
        emdb_id: str,
    ):
        self.emdb_id = emdb_id

    def fetch_emdb_map(self, max_retries: int = 3) -> bytes:
        """Fetch density map from EMDB

        Args:
            max_retries (int, optional): Maximum number of retries.

        Returns:
            bytes: Density map file content
        """

        emdb_id_ = self.emdb_id.lower().replace("-", "_")

        EMDB_BASE_URL = f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/"
        EMDB_MAP_URL = f"{EMDB_BASE_URL}{self.emdb_id}/map/{emdb_id_}.map.gz"

        req_sess = request_session(max_retries=max_retries)
        response = req_sess.get(EMDB_MAP_URL)

        if response.status_code == 200:
            print("Successfully fetched EMDB map for given EMDB id")
            return response.content

        else:
            raise Exception("Error while requesting EMDB map for given EMDB id")

    def fetch_emdb_mask(self, mask_name, max_retries=3):

        if mask_name is None:
            raise ValueError("Mask name must be provided to fetch EMDB mask")

        EMDB_BASE_URL = f"https://ftp.ebi.ac.uk/pub/databases/emdb/structures/"
        EMDB_MASK_URL = f"{EMDB_BASE_URL}{self.emdb_id}/masks/{mask_name}.map"

        req_sess = request_session(max_retries=max_retries)
        response = req_sess.get(EMDB_MASK_URL)

        if response.status_code == 200:
            print("Successfully fetched EMDB mask for given EMDB id and mask name")
            return response.content

        else:
            raise Exception("Error while requesting EMDB mask for given EMDB id and mask name")