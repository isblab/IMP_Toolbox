import json
import requests
import platform
import os

def get_dropbox_path():
    system_spec_paths = {
    "Windows": f"C:\\Users\\{os.getlogin()}\\Dropbox",
    "Linux": f"/home/{os.getlogin()}/Dropbox"
}

    if platform.system() in system_spec_paths:
        path_to_dropbox = system_spec_paths[platform.system()]
    else:
        path_to_dropbox = input("Enter the path to your Dropbox folder: ")
    return path_to_dropbox

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
    return req_sess

def request_result(get_request, uniprot_id, ignore_error=False):
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
        print(f"Error while requesting {uniprot_id}") if not ignore_error else None
        return None

def write_json(file_path, data):
    """Write data to a json file

    Args:
        file_path (str): path to json file
        data (dict): data to write
    """
    with open(file_path, "w") as f:
        json.dump(data, f)

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


def read_fasta(fasta_file):
    """
    Read a fasta file and return a dictionary of sequences
    """

    all_sequences = {}

    with open(fasta_file, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith(">"):
            seq_id = line[1:].strip()
        else:
            seq = line.strip()
            all_sequences[seq_id] = seq if seq_id not in all_sequences else all_sequences[seq_id] + seq

    return all_sequences

def get_key_from_res_range(res_range):
    """Returns a residue range string from a list of residue numbers.

    Args:
        res_range (list): List of residue numbers, e.g., [1, 2, 3, 5, 6, 7]

    Returns:
        str: Residue range string, e.g., "1-3,5-7"
    """

    if not res_range:
        return ""
    res_range = sorted(res_range)
    ranges = []
    start = prev = res_range[0]
    for num in res_range[1:]:
        if num == prev + 1:
            prev = num
        else:
            ranges.append(f"{start}-{prev}") if start != prev else ranges.append(str(start))
            start = prev = num
    if start == prev:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{prev}")
    return ",".join(ranges)

def save_df(df, file_path, index=False, header=True, sep=","):
    """Save a pandas DataFrame to a file

    Args:
        df (pandas.DataFrame): DataFrame to save
        file_path (str): path to save the DataFrame
    """
    df.to_csv(file_path, index=index, header=header, sep=sep)

def generate_cmap(n):
    """
    Generate a custom colormap with n colors.
    """

    import random
    import time
    colors = []

    start = time.time()

    while len(colors) < n:

        color = "%06x" % random.randint(0, 0xFFFFFF)

        if f"#{color}" not in colors:
            colors.append(f"#{color}")

        if time.time() - start > 10:
            break

    return colors

def str_join( prots: list, sep: str ):
	"""
	Given a list of strings, join them by the provided separator (sep).
	"""
	return f"{sep}".join( prots )


def str_split( prot_pair: str, sep: str ):
	"""
	Split the given string by the  provided separator (sep).
	"""
	return prot_pair.split( sep )

