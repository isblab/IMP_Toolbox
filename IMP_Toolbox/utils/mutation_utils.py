import warnings
import re
import os
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
from IMP_Toolbox.utils.api_helpers import (
    request_session,
    request_result,
)
from IMP_Toolbox.utils.file_helpers import (
    write_json,
    read_json,
)
from IMP_Toolbox.constants.mutation_constants import (
    MISSENSE_REGEX,
    TRUNCATION_NOTATIONS,
    API_URLS,
    AMINO_ACID_MAP,
)

def get_uniprot_variants(
    uniprot_id: str,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
) -> dict | None:
    """ Fetch variants for a given UniProt ID from the UniProt database.

    ## Arguments:

    - **uniprot_id (str)**:<br />
        The UniProt ID to fetch variants for.

    - **api_url (str)**:<br />
        #TODO: shouldn't be passed as an argument
        The API URL to fetch variants from.

    - **api_parameters (dict, optional):**:<br />
        A dictionary of parameters to replace in the API URL.

    - **ignore_error (bool, optional):**:<br />
        Whether to ignore errors and return None instead of raising an exception.

    - **max_retries (int, optional):**:<br />
        The maximum number of retries for the API request in case of failure.

    ## Returns:

    - **dict | None**:<br />
        The response from the API request, or None if ignore_error is True and an error occurs.
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
    """ Fetch the nucleotide sequence for a given RefSeq ID from the NCBI nuccore
    database.

    ## Arguments:

    - **ref_seq_id (str)**:<br />
        The RefSeq ID to fetch the nucleotide sequence for.

    - **api_url (str)**:<br />
        #TODO: shouldn't be passed as an argument
        The API URL to fetch the nucleotide sequence from.

    - **api_parameters (dict, optional):**:<br />
        A dictionary of parameters to replace in the API URL.

    - **ignore_error (bool, optional):**:<br />
        Whether to ignore errors and return None instead of raising an exception.

    - **max_retries (int, optional):**:<br />
        The maximum number of retries for the API request in case of failure.

    ## Returns:

    - **str**:<br />
        The nucleotide sequence for the given RefSeq ID, or None if
        `ignore_error` is True and an error occurs.
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
) -> list | str | None:
    """ Fetch variants for a given gene name from LOVD.

    ## Arguments:

    - **gene_name (str)**:<br />
        Valid gene name to fetch variants for.

    - **api_url (str)**:<br />
        The API URL to fetch variants from.

    - **api_parameters (dict, optional):**:<br />
        A dictionary of parameters to replace in the API URL.

    - **ignore_error (bool, optional):**:<br />
        Whether to ignore errors and return None instead of raising an exception.

    - **max_retries (int, optional):**:<br />
        The maximum number of retries for the API request in case of failure.

    - **format (str, optional):**:<br />
        The format to return the variants in. Can be "json" or "txt". Defaults to "json".

    ## Returns:

    - **list | str | None**:<br />
        List of variants for the given gene name if found, or a string if `format`
        is "txt", or None if `ignore_error` is True and an error occurs.
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

def is_missense_mutation(
    p_mutation: str | None,
    allow_truncation: bool = False,
    ignore_warnings: bool = False,
) -> bool:
    """ Check if a protein mutation string represents a missense mutation

    ## Arguments:

    - **p_mutation (str | None)**:<br />
        Protein mutation string in the format "p.(Ala123Val)".

    - **allow_truncation (bool, optional):**:<br />
        Whether to allow truncation mutations (e.g. p.Arg150Ter) to be considered

    ## Returns:

    - **bool**:<br />
        True if the mutation is a missense mutation, False otherwise.
    """

    if p_mutation is None:
        return False

    pattern = None
    for key, value in MISSENSE_REGEX.items():
        if value["condition"](p_mutation):
            pattern = value["regex"]
            break

    if pattern is None and ignore_warnings is False:
        warnings.warn(f"No regex pattern found for mutation: {p_mutation}")
        return False

    match = re.match(pattern, p_mutation)

    if not allow_truncation and match and len(match.groups()) == 3:
        if match.group(3) in TRUNCATION_NOTATIONS:
            return False

    #NOTE: the following will miss missense mutations where
    # more than one residue is mutated
    # e.g. p.Gly1094_His1095delinsValAsn in DSG2
    return_val = match is not None and len(match.groups()) == 3
    if return_val is False:
        warnings.warn(
            f"""
            Protein mutation not found or not missense for mutation:
            {p_mutation=}
            Molecular consequence might be "missense", but more than
            one residue might be mutated. Skipping...
            """
        )

    return return_val

def split_missense_mutation(
    p_mutation: str,
    return_type: str = "all",
    ignore_warnings: bool = False,
) -> tuple | str | int | None:
    """ Split a missense mutation string into its components.

    Assuming the mutation string is in one of the recognized formats,
    this function extracts the wild-type amino acid, residue number,
    and mutant amino acid.

    ## Arguments:

    - **p_mutation (str)**:<br />
        Mutation string (e.g. Arg150Trp).

    - **return_type (str, optional):**:<br />
        Type of information to return:
        - "all": return (wt_aa, res_num, mut_aa)
        - "wt_aa": return wild-type amino acid
        - "res_num": return residue number
        - "mut_aa": return mutant amino acid
        Defaults to "all".

    - **ignore_warnings (bool, optional):**:<br />
        Whether to ignore warnings.

    ## Returns:

    - **tuple | str | int | None**:<br />
        Depending on `return_type`, returns:
        - (wt_aa, res_num, mut_aa) if `return_type` is "all"
        - wild-type amino acid if `return_type` is "wt_aa"
        - residue number (int) if `return_type` is "res_num"
        - mutant amino acid if `return_type` is "mut_aa"
        - None if the mutation string could not be parsed
    """

    pattern = None
    for _key, value in MISSENSE_REGEX.items():
        if value["condition"](p_mutation):
            pattern = value["regex"]
            break

    if pattern is None:
        if ignore_warnings is False:
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

    elif ignore_warnings is False:
        warnings.warn(f"Could not parse mutation: {p_mutation}")

    return None

def get_ncbi_ref_seq(
    ncbi_ref_seq_id: str,
    ref_seq_file: str
) -> str:
    """ Handle fetching and parsing of reference sequence from NCBI.

    ## Arguments:

    - **ncbi_ref_seq_id (str)**:<br />
        NCBI RefSeq ID to fetch the reference sequence for.

    - **ref_seq_file (str)**:<br />
        Path to save the RefSeq XML file. If the file already exists, it will be
        read instead of fetching from NCBI again.

    ## Returns:

    - **str**:<br />
        Protein sequence as a string.
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

    ## Arguments:

    - **hgmd_cdna_html_path (str)**:<br />
        Path to the HGMD cDNA HTML file.

    ## Returns:

    - **str**:<br />
        Amino acid sequence as a string.
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