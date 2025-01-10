import json
import requests
import platform
import os

def get_dropbox_path():
    """
    Determine the path to the Dropbox folder based on the current operating system.
    
    This function provides a cross-platform method to locate the Dropbox folder by using predefined paths for Windows and Linux. 
    If the operating system is not recognized, it prompts the user to manually enter the Dropbox folder path.
    
    Returns:
        str: The full path to the Dropbox folder on the current system.
    
    Raises:
        No explicit exceptions are raised, but user input is required for unrecognized systems.
    
    Examples:
        >>> dropbox_path = get_dropbox_path()
        >>> print(dropbox_path)
        # On Windows: C:\Users\username\Dropbox
        # On Linux: /home/username/Dropbox
    """
    system_spec_paths = {
    "Windows": f"C:\\Users\\{os.getlogin()}\\Dropbox",
    "Linux": f"/home/{os.getlogin()}/Dropbox"
}

    if platform.system() in system_spec_paths.keys():
        path_to_dropbox = system_spec_paths[platform.system()]
    else:
        path_to_dropbox = input("Enter the path to your Dropbox folder: ")
    return path_to_dropbox

def request_session(max_retries=3):
    """
    Create a requests session with configurable retry mechanism for HTTPS requests.
    
    This function initializes a requests Session object and configures it with an HTTPAdapter
    to automatically retry failed HTTPS requests. The number of retry attempts can be customized.
    
    Parameters:
        max_retries (int, optional): Maximum number of retry attempts for failed HTTPS requests.
            Defaults to 3. Determines how many times a request will be retried before giving up.
    
    Returns:
        requests.sessions.Session: A configured requests session with retry capabilities for HTTPS endpoints.
    
    Example:
        session = request_session(max_retries=5)  # Create a session with 5 retry attempts
        response = session.get('https://example.com')
    """
    req_sess = requests.Session()
    req_sess.mount(
        "https://",
        requests.adapters.HTTPAdapter(max_retries=max_retries)
    )
    return req_sess

def request_result(get_request, uniprot_id, ignore_error=False):
    """
    Process the result of an HTTP GET request for a UniProt entry.
    
    Parameters:
        get_request (requests.models.Response): HTTP GET request response object
        uniprot_id (str): UniProt identifier for the requested entry
        ignore_error (bool, optional): Flag to suppress error messages. Defaults to False.
    
    Returns:
        dict or bytes or None: JSON response if successful, raw content if JSON decoding fails,
        or None if request was unsuccessful
    
    Raises:
        json.decoder.JSONDecodeError: If response cannot be decoded as JSON
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
    """
    Write data to a JSON file.
    
    Args:
        file_path (str): Path to the output JSON file where data will be written.
        data (dict): Dictionary containing the data to be serialized and saved.
    
    Raises:
        IOError: If the file cannot be opened or written.
        TypeError: If the data cannot be serialized to JSON.
    """
    with open(file_path, "w") as f:
        json.dump(data, f)

def read_json(file_path):
    """
    Read a JSON file and return its contents as a dictionary.
    
    Args:
        file_path (str): Path to the JSON file to be read.
    
    Returns:
        dict: Parsed JSON data from the file.
    
    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the file contains invalid JSON data.
    """
    with open(file_path, "r") as f:
        data = json.load(f)
    return data


def read_fasta(fasta_file):
    """
    Read a FASTA file and return a dictionary of sequences.
    
    Parameters:
        fasta_file (str): Path to the FASTA file to be read.
    
    Returns:
        dict: A dictionary where keys are sequence identifiers and values are their corresponding sequences.
              Supports multi-line sequences for the same identifier.
    
    Raises:
        FileNotFoundError: If the specified FASTA file cannot be found.
        IOError: If there are issues reading the file.
    
    Example:
        sequences = read_fasta('proteins.fasta')
        # Returns: {'seq1': 'MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTDGLLYGSQTPNEECLFLERLEENHYNTYISKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV', ...}
    """

    all_sequences = {}

    with open(fasta_file, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith(">"):
            seq_id = line[1:].strip()
        else:
            seq = line.strip()
            all_sequences[seq_id] = seq if seq_id not in all_sequences else all_sequences[seq_id] + seq

    return all_sequences


def get_key_from_res_range(res_range):
    """
    Convert a list of residue numbers into a compact range string representation.
    
    Parameters:
        res_range (list): A list of sequential residue numbers to be converted.
    
    Returns:
        str: A compact string representation of residue ranges, such as "1-3,5-7".
             Returns an empty string if the input list is empty.
    
    Examples:
        >>> get_key_from_res_range([1, 2, 3, 5, 6, 7])
        '1-3,5-7'
        >>> get_key_from_res_range([10, 11, 12, 15, 16])
        '10-12,15-16'
        >>> get_key_from_res_range([])
        ''
    """
    if len(res_range) == 0:
        return ""
    key = str(res_range[0])
    for i, res in enumerate(res_range):
        if i == 0:
            continue
        elif res_range[i-1] != res - 1:
            key += f",{res}"
        if i+1 < len(res_range):
            if res_range[i+1] != res + 1:
                key += f"-{res_range[i]}"
        if i+1 == len(res_range):
            key += f"-{res}"
    return key

def save_df(df, file_path, index=False, header=True, sep=","):
    """
    Save a pandas DataFrame to a CSV file with configurable export options.
    
    Parameters:
        df (pandas.DataFrame): The DataFrame to be saved to a file
        file_path (str): The destination file path for the CSV output
        index (bool, optional): Whether to include the DataFrame index in the output. Defaults to False.
        header (bool, optional): Whether to include column names in the output. Defaults to True.
        sep (str, optional): Delimiter to use between values. Defaults to ",".
    
    Example:
        save_df(my_dataframe, "/path/to/output.csv")
        save_df(my_dataframe, "/path/to/output.csv", index=True, sep="\t")
    """
    df.to_csv(file_path, index=index, header=header, sep=sep)