import warnings
import re
import xmltodict
import os
from tqdm import tqdm
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
from IMP_Toolbox.utils_imp_toolbox.api_helpers import (
    request_session,
    request_result,
)
from IMP_Toolbox.utils_imp_toolbox.file_helpers import (
    write_json,
    read_json,
)
from IMP_Toolbox.pre_processing.mutations.mutation_constants import (
    MISSENSE_REGEX,
    TRUNCATION_NOTATIONS,
    API_URLS,
    AMINO_ACID_MAP,
)

def get_af_missense_data(
    uniprot_id: str,
    api_url: str,
    api_parameters: dict = {},
    return_type: str = "csv",
    ignore_error: bool = False,
    max_retries: int = 3,
):
    """ Fetch alpha missense variants for a given Uniprot ID from the
    AlphaMissense database.

    Args:
        uniprot_id (str): Valid Uniprot ID to fetch variants for.
        ignore_error (bool, optional): Defaults to False.
        max_retries (int, optional): Defaults to 3.

    Returns:
        list: List of alpha missense variants for the given Uniprot ID if found
    """

    for key, value in api_parameters.items():
        api_url = api_url.replace(key, value)

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(api_url)

    if response.status_code == 200:

        if return_type == "json":
            return response.json()

        elif return_type == "csv":
            return response.text

        else:
            raise ValueError(f"Invalid return type: {return_type}")

    elif ignore_error:
        return None

    else:
        raise ValueError(
            f"""Error fetching AlphaMissense data for {uniprot_id}
            {response.status_code}: {response.text}
            """
        )

def get_uniprot_variants(
    uniprot_id: str,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
):
    """ Fetch variants for a given Uniprot ID from the UniProt database.

    Args:
        uniprot_id (str): Valid Uniprot ID to fetch variants for.
        ignore_error (bool, optional): Defaults to False.
        max_retries (int, optional): Defaults to 3.

    Returns:
        list: List of variants for the given Uniprot ID if found
    """

    for key, value in api_parameters.items():
        api_url = api_url.replace(key, value)

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(api_url)

    return request_result(response, uniprot_id, ignore_error=ignore_error)

def fetch_nuccore_ref_seq_from_id(
    ref_seq_id: str,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
) -> str:
    """ Fetch the nucleotide sequence for a given RefSeq ID from the NCBI
    nuccore database.

    Args:
        ref_seq_id (str):
            The RefSeq ID to fetch the nucleotide sequence for.

    Returns:
        str:
            The nucleotide sequence for the given RefSeq ID.
    """

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(api_url, params=api_parameters)

    if response.status_code == 200:
        result = response.text

    elif ignore_error:
        return None

    else:
        raise ValueError(
            f"""Error fetching RefSeq ID {ref_seq_id}
            {response.status_code}: {response.text}
            """
        )

    root = ET.fromstring(result)
    p_sequence = root.find('.//*NCBIeaa').text.replace("\n", "").strip()

    return {ref_seq_id: p_sequence}

def get_lovd_variants_for_gene(
    gene_name: str,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
    format="json"
):
    """ Fetch variants for a given gene name from the LOVD.

    Args:
        gene_name (str):
            Valid gene name to fetch variants for.
        ignore_error (bool, optional):
            Defaults to False.
        max_retries (int, optional):
            Defaults to 3.

    Returns:
        list:
            List of variants for the given gene name if found.
    """

    for key, value in api_parameters.items():
        api_url = api_url.replace(key, value)

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(api_url)

    if response.status_code == 200:

        if format == "json":
            return response.json()

        elif format == "txt":
            return response.text

        else:
            raise ValueError(f"Invalid format: {format}")

    elif ignore_error:
        return None

    else:
        raise ValueError(
            f"""Error fetching variants for gene {gene_name} in LOVD.
            {response.status_code}: {response.text}
            """
        )

def get_variant_ids_from_clinvar(
    gene_name:str,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
) -> list:
    """ Fetch variant IDs from ClinVar for a given gene name.

    Args:
        gene_name (str): The gene name to search for.

    Returns:
        list: List of ClinVar variant IDs.
    """

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(api_url, params=api_parameters)
    result = request_result(response, gene_name)

    if result is None:
        # try without returnmode
        api_parameters.pop("retmode")
        response = req_sess.get(api_url, params=api_parameters)

        if response.status_code == 200:
            result = response.text

        elif ignore_error:
            return []

        else:
            raise ValueError(
                f"""Error fetching variant ids for {gene_name}
                {response.status_code}: {response.text}
                """
            )

    if isinstance(result, str):
        # parse xml
        root = ET.fromstring(result)
        id_list = root.find("IdList")

        if id_list is None:
            return []

        ids = [id_elem.text for id_elem in id_list.findall("Id")]

        return ids

    elif isinstance(result, dict):
        return result.get("esearchresult", {}).get("idlist", [])

    return []

def get_variant_details_from_clinvar(
    variant_ids:list,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
) -> dict:
    """ Fetch variant details from ClinVar for a list of variant IDs.
    Note: If you parallelize it, make sure to not exceed API's rate limits.
    We get info in XML format from the API, which we convert and store as JSON.

    Args:
        variant_ids (list): List of ClinVar variant IDs.

    Returns:
        dict: Dictionary with variant IDs as keys and their details as values.
    """

    variant_id_batches = [
        variant_ids[i:i + 100] for i in range(0, len(variant_ids), 100)
    ]

    variant_info_dict = {}

    for idx, batch in enumerate(tqdm(variant_id_batches)):

        api_parameters["id"] = ",".join(batch)

        req_sess = request_session(max_retries=max_retries)
        response = req_sess.get(api_url, params=api_parameters)

        if response.status_code == 200:
            result = response.text
        else:
            if ignore_error:
                return {}
            else:
                raise ValueError(f"Error fetching variant details for batch {idx}")

        root = ET.fromstring(result)
        variation_info = root.findall("VariationArchive")
        variant_ids = [
            var_id_elem.get("VariationID") for var_id_elem in variation_info
        ]

        for idx, variation_id in enumerate(variant_ids):
            variant_info_dict[variation_id] = xmltodict.parse(
                ET.tostring(variation_info[idx])
            )

    return variant_info_dict

def is_missense_mutation(
    p_mutation: str,
    allow_truncation: bool = False,
) -> bool:
    """ Check if a given protein mutation string represents a missense mutation.

    Args:
        p_mutation (str): Protein mutation string in the format "p.(Ala123Val)".

    Returns:
        bool: True if the mutation is a missense mutation, False otherwise.
    """

    pattern = None
    for key, value in MISSENSE_REGEX.items():
        if value["condition"](p_mutation):
            pattern = value["regex"]
            break

    if pattern is None:
        warnings.warn(f"No regex pattern found for mutation: {p_mutation}")
        return False

    match = re.match(pattern, p_mutation)

    if not allow_truncation and match and len(match.groups()) == 3:
        if match.group(3) in TRUNCATION_NOTATIONS:
            return False

    return match is not None and len(match.groups()) == 3

def split_missense_mutation(
    p_mutation: str,
    return_type: str = "all",
) -> tuple | str | int | None:
    """ Extract the residue number from a mutation string.
    Example: Arg150Trp -> 150

    Args:
        p_mutation (str): mutation string (e.g. Arg150Trp)

    Returns:
        int: residue number
    """

    pattern = None
    for key, value in MISSENSE_REGEX.items():
        if value["condition"](p_mutation):
            pattern = value["regex"]
            break

    if pattern is None:
        warnings.warn(f"No regex pattern found for mutation: {p_mutation}")
        return None

    match = re.match(pattern, p_mutation)

    if match:
        if return_type == "all":
            return match.groups()

        elif return_type == "wt_aa":
            return match.group(1)

        elif return_type == "res_num":
            return int(match.group(2))

        elif return_type == "mut_aa":
            return match.group(3)

        else:
            raise ValueError(f"Invalid return type: {return_type}")
    else:
        warnings.warn(f"Could not parse mutation: {p_mutation}")

    return None

def get_ncbi_ref_seq(
    ncbi_ref_seq_id: str,
    ref_seq_file: str
):
    """ Handle fetching and parsing of reference sequence from NCBI.

    Args:
        ncbi_ref_seq_id (str): NCBI RefSeq ID
        ref_seq_file (str): Path to save the RefSeq XML file

    Returns:
        str: Protein sequence
    """

    if not os.path.exists(ref_seq_file):
        p_sequence_dict = fetch_nuccore_ref_seq_from_id(
            ncbi_ref_seq_id,
            api_url=API_URLS["ncbi_efetch"],
            api_parameters={
                "db": "nuccore",
                "id": ncbi_ref_seq_id,
                "rettype": "native",
                "retmode": "xml",
            },
            ignore_error=False,
            max_retries=3,
        )
        write_json(ref_seq_file, p_sequence_dict)

    else:
        p_sequence_dict = read_json(ref_seq_file)

    p_sequence = p_sequence_dict.get(ncbi_ref_seq_id, None)

    return p_sequence

def get_protein_sequence_from_hgmd_cdna_html(hgmd_cdna_html_path:str) -> str:
    """ Extracts the amino acid sequence from an HGMD cDNA HTML file.

    Args:
        hgmd_cdna_html_path (str): Path to the HGMD cDNA HTML file.

    Returns:
        str: Amino acid sequence as a string.
    """

    aa_sequence = []

    with open(hgmd_cdna_html_path, "r") as file:
        html = file.read()

    htmlParse = BeautifulSoup(html, 'html.parser')
    content = htmlParse.find_all("div", class_="content")

    if len(content) != 1:
        raise ValueError("Expected exactly one content div in the HTML file.")

    content = content[0].get_text(separator="|").splitlines()[1:-2]

    for line in content:

        amino_acid = line.split("|")[2]
        amino_acid = AMINO_ACID_MAP.get(amino_acid, None)

        # Stop at termination codon
        if f'{line.split('|')[2]}' == "Ter":
            break

        if amino_acid is None:
            warnings.warn(
                f"Amino acid '{line.split('|')[2]}' not found in mapping."
            )
            continue

        aa_sequence.append(amino_acid)

    return "".join(aa_sequence)