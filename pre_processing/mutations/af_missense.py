import os
from typing import Any
import yaml
import warnings
import argparse
import pandas as pd
from io import StringIO
from pprint import pprint
from IMP_Toolbox.pre_processing.sequence.Sequence import FetchSequences
from IMP_Toolbox.utils_imp_toolbox.file_helpers import read_fasta
from IMP_Toolbox.utils_imp_toolbox.special_helpers import get_mapped_residue
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import fasta_str_to_dict
from IMP_Toolbox.pre_processing.mutations.utils_mutation import (
    split_missense_mutation,
)
from IMP_Toolbox.utils_imp_toolbox.api_helpers import (
    request_session,
)
from IMP_Toolbox.pre_processing.mutations.mutation_constants import (
    AMINO_ACID_MAP,
    AFM_PATHOGENICITY,
    API_URLS,
    AF_MISSENSE_CSV_SUFFIX,
    AF_MISSENSE_PAIR_ALN_SUFFIX,
    AF_MISSENSE_AA_SUBSTITUTIONS_TSV,
)

def get_af_missense_data_online(
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

def get_af_missense_data_offline(
    uniprot_ids: list,
    tsv_path: str,
):
    try:
        import duckdb
    except ImportError:
        raise ImportError(
            "duckdb is required for offline AlphaMissense data retrieval. "
            "Please install it via"
            "'pip install duckdb' or 'conda install -c conda-forge duckdb'."
        )

    assert os.path.exists(tsv_path), (
        f"AlphaMissense aa substitutions TSV file not found at {tsv_path}."
    )

    uniprot_ids = list(set(uniprot_ids))
    uids_df = pd.DataFrame({"uniprot_id": uniprot_ids})
    query = f"""
    SELECT t.*
    FROM '{tsv_path}' t
    JOIN uids_df u USING (uniprot_id)
    """

    result_df = duckdb.sql(query).fetchdf()

    # split into multiple dataframes per uniprot id
    result_dfs = {}
    for uid in uniprot_ids:
        result_dfs[uid] = result_df[result_df["uniprot_id"] == uid].copy()

    return result_dfs

def fetch_af_missense_data(
    alpha_missense_dir: str,
    uniprot_bases: list,
    mode: str = "online",
    overwrite: bool = False,
    **kwargs,
):
    """ Fetch AlphaMissense data for a list of UniProt IDs.

    Args:

        alpha_missense_dir (str):
            Directory to save AlphaMissense CSV files.

        uniprot_bases (list):
            List of UniProt base IDs (without isoform suffix).

        mode (str, optional):
            Mode to fetch data: 'online' or 'offline'.

        overwrite (bool, optional):
            Whether to overwrite existing files. Defaults to False.

    Kwargs:

        af_missense_tsv (str, optional):
            Path to AlphaMissense aa substitutions TSV file for offline mode.

    Yields:
        tuple:
            pd.DataFrame: AlphaMissense dataframe for the UniProt ID.
            str: UniProt base ID.
    """

    os.makedirs(alpha_missense_dir, exist_ok=True)

    remainder_bases = get_remainder_uniprot_bases(
        alpha_missense_dir,
        uniprot_bases,
        overwrite
    )

    if len(remainder_bases) == 0:
        print("All AlphaMissense data files already exist. Nothing to fetch.")
        return

    if mode == "online":

        for uniprot_base in sorted(uniprot_bases):

            if uniprot_base not in remainder_bases:
                af_missense_file = os.path.join(
                    alpha_missense_dir,
                    f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
                )
                af_missense_df = pd.read_csv(af_missense_file)
                yield af_missense_df, uniprot_base
                continue

            try:
                af_missense_csv = get_af_missense_data_online(
                    uniprot_base,
                    api_url=API_URLS["af_missense_csv"],
                    api_parameters={"uniprot_id": uniprot_base},
                    return_type="csv",
                    ignore_error=False,
                )
                af_missense_df = pd.read_csv(StringIO(af_missense_csv))
                remainder_bases.remove(uniprot_base)

                yield af_missense_df, uniprot_base

            except Exception as e:
                warnings.warn(
                    f"""Error fetching AlphaMissense data for {uniprot_base}:
                    {e}"""
                )
                yield pd.DataFrame(), uniprot_base

    elif mode == "offline":

        af_missense_dfs = get_af_missense_data_offline(
            uniprot_ids=remainder_bases,
            tsv_path=kwargs.get(
                "af_missense_tsv",
                AF_MISSENSE_AA_SUBSTITUTIONS_TSV,
            ),
        )

        # for uniprot_base in remainder_bases:
        for uniprot_base in sorted(uniprot_bases):

            if uniprot_base not in remainder_bases:
                af_missense_file = os.path.join(
                    alpha_missense_dir,
                    f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
                )
                af_missense_df = pd.read_csv(af_missense_file)
                yield af_missense_df, uniprot_base
                continue

            af_missense_df = af_missense_dfs.get(uniprot_base, pd.DataFrame())

            yield af_missense_df, uniprot_base

def get_remainder_uniprot_bases(
    alpha_missense_dir: str,
    uniprot_bases: list,
    overwrite: bool=False,
):
    """ Get the list of UniProt bases for which AlphaMissense data

    Args:

        alpha_missense_dir (str):
            Directory to save AlphaMissense CSV files.

        uniprot_bases (list):
            List of UniProt base IDs (without isoform suffix).

        overwrite (bool, optional):
            Whether to overwrite existing files. Defaults to False.

    Returns:
        list:
            List of UniProt bases for which AlphaMissense data needs to be fetched.
    """

    remainder_bases = []

    if overwrite:
        return uniprot_bases

    for uniprot_base in uniprot_bases:
        af_missense_file = os.path.join(
            alpha_missense_dir, f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
        )

        if os.path.exists(af_missense_file):
            print(f"File {af_missense_file} already exists. Skipping...")
        else:
            remainder_bases.append(uniprot_base)

    return remainder_bases

def get_af_missense_attribute(
    af_missense_dict: dict,
    attribute: str,
    p_name: str,
    p_mutation: str,
) -> Any:
    """ Get a specific attribute for a given protein mutation from
    the AlphaMissense variants dictionary.

    Args:

        af_missense_dict (dict):
            AlphaMissense variants dictionary

        attribute (str):
            Attribute to retrieve (e.g., "patho_score", "v_pathogenicity")

        p_name (str):
            Protein name

        p_mutation (str):
            Protein mutation in the format "Arg123Cys"

    Returns:
        Any:
            Value of the requested attribute for the given protein mutation,
            or None if not found.
    """

    if p_name not in af_missense_dict:
        warnings.warn(f"Protein {p_name} not found in AlphaMissense data.")
        return None

    if p_mutation not in af_missense_dict[p_name]:
        warnings.warn(
            f"Mutation {p_mutation} not found for protein {p_name} "
            f"in AlphaMissense data."
        )
        return None

    if attribute not in af_missense_dict[p_name][p_mutation]:
        warnings.warn(
            f"Attribute {attribute} not found for "
            f"protein {p_name} mutation {p_mutation}."
        )
        return None

    return af_missense_dict[p_name][p_mutation][attribute]

def af_missense_df_to_dict(
    p_name: str,
    af_missense_df: pd.DataFrame,
    af_missense_dict: dict = {},
    afm_psa_map: dict = {},
    ignore_warnings: bool = False,
) -> dict:
    """ Convert AlphaMissense dataframe to a dictionary.

    For each UniProt ID, AF-missense data is available in a CSV file.
    This function converts the dataframe to a dictionary for easier access.

    The data is stored in a nested dictionary format in which the outer key
    is the protein name and the inner key is the missense mutation in the
    format "Arg123Cys".

    Additionally, if pairwise alignment map is provided, the residue
    numbers are mapped to the modeled (query) sequence residue numbers.

    Args:
        p_name (str): Protein name
        af_missense_df (pd.DataFrame): AlphaMissense dataframe
        af_missense_dict (dict, optional): AlphaMissense variants dictionary
        afm_psa_map (dict): Residue mapping of query and subject sequences from
            pairwise alignment

    Returns:
        dict: Updated dictionary of AlphaMissense variants

    Example:
        >>> df = pd.DataFrame(
        ... [{"protein_variant": "M1A", "am_pathogenicity": 0.5, "am_class": "LPath"}]
        ... )
        >>> af_missense_dict = af_missense_df_to_dict("PKP2", df)
        >>> pprint(af_missense_dict)
        {'PKP2': {'Met1Ala': {'mut_aa': 'Ala',
                              'patho_score': np.float64(0.5),
                              'res_num': 1,
                              'v_pathogenicity': 'Likely pathogenic',
                              'wt_aa': 'Met'}}}
    """

    if p_name not in af_missense_dict:
        af_missense_dict[p_name] = {}

    protein_variants = af_missense_df["protein_variant"].values
    patho_scores = af_missense_df["am_pathogenicity"].values
    v_pathogenicities = af_missense_df["am_class"].values

    for p_variant, patho_score, v_pathogenicity in zip(
        protein_variants, patho_scores, v_pathogenicities
    ):

        wt_aa, res_num, mut_aa = split_missense_mutation(
            p_variant, return_type="all"
        )

        if any(aa not in AMINO_ACID_MAP for aa in [wt_aa, mut_aa]) and ignore_warnings is False:
            warnings.warn(
                f"""
                Invalid amino acid in {p_variant} for protein {p_name}.
                Got {wt_aa=}, {mut_aa=}.
                Skipping...
                """
            )
            continue

        wt_aa = AMINO_ACID_MAP.get(wt_aa, wt_aa)
        mut_aa = AMINO_ACID_MAP.get(mut_aa, mut_aa)
        res_num_int = int(res_num)

        res_num_mapped, warn_msg = get_mapped_residue(
            psa_map=afm_psa_map,
            codon_number=res_num_int,
        )

        if ignore_warnings is False and len(warn_msg) > 0:
            warnings.warn(f"""{warn_msg} (Protein: {p_name})""")
            continue

        if res_num_mapped is None:
            continue

        p_variant_key = f"{wt_aa}{res_num_mapped}{mut_aa}"

        af_missense_dict[p_name][p_variant_key] = {
            "patho_score": patho_score,
            "v_pathogenicity": AFM_PATHOGENICITY.get(v_pathogenicity),
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "res_num": res_num_mapped,
        }

    return af_missense_dict

def get_fasta_dict_for_af_missense(protein_uniprot_map: dict):
    """ Fetch FASTA sequences for AlphaMissense reference sequences.

    Args:

        protein_uniprot_map (dict):
            Mapping of protein names to UniProt IDs.

    Returns:
        dict:
            Dictionary of FASTA sequences for AlphaMissense reference sequences.
    """

    # For AF-missense, the reference sequence is the Uniprot sequence
    # for non-isoform-specific uniprot ids
    uniprot_ids = list(set((protein_uniprot_map.values())))
    uniprot_ids = list(set([uid.split("-")[0] for uid in uniprot_ids]))
    fetchit = FetchSequences(uniprot_ids=uniprot_ids)
    fasta_str = fetchit.query_uniprot_api_for_sequences()
    fasta_str = fetchit.only_uniprot_id_as_name(fasta_str)
    fasta_dict = fasta_str_to_dict(fasta_str)

    return fasta_dict

def fetch_fasta_dict_for_af_missense(
    savepath: str,
    protein_uniprot_map: dict,
    overwrite: bool = False,
):
    if os.path.exists(savepath) and not overwrite:
        print(f"FASTA file {savepath} already exists. Loading...")
        fasta_dict = read_fasta(savepath)
        return fasta_dict

    os.makedirs(os.path.dirname(savepath), exist_ok=True)
    fasta_dict = get_fasta_dict_for_af_missense(protein_uniprot_map)

    with open(savepath, "w") as f:
        for uid, seq in fasta_dict.items():
            f.write(f">{uid}\n{seq}\n")

    print(f"Saved AlphaMissense FASTA sequences to {savepath}")
    return fasta_dict

def export_af_missense_data(
    alpha_missense_dir: str,
    af_missense_df_gen: iter,
    overwrite: bool = False,
):
    """ Export AlphaMissense dataframes to CSV files.

    Args:

        alpha_missense_dir (str):
            Directory to save AlphaMissense CSV files.

        af_missense_df_gen (iter):
            Generator yielding tuples of AlphaMissense dataframe and
            UniProt base ID.

        overwrite (bool, optional):
            Whether to overwrite existing files. Defaults to False.
    """

    for af_missense_df, uniprot_base in af_missense_df_gen:

        if af_missense_df.empty:
            print(f"No AlphaMissense data found for {uniprot_base}. Skipping...")
            continue

        af_missense_file = os.path.join(
            alpha_missense_dir, f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
        )

        if os.path.exists(af_missense_file) and not overwrite:
            print(f"File {af_missense_file} already exists. Skipping save...")
            continue

        af_missense_df.to_csv(af_missense_file, index=False)

        print(f"Saved AlphaMissense variants to {af_missense_file}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Fetch and save Alpha-Missense variants."
    )
    parser.add_argument(
        "--config_file",
        type=str,
        default="/home/omkar/Projects/cardiac_desmosome/input/config.yaml",
        help="Path to configuration YAML file.",
    )
    parser.add_argument(
        "--alpha_missense_dir",
        type=str,
        default="/home/omkar/Projects/cardiac_desmosome/data/literature_parsing/mutations/alpha_missense",
        help="Directory to save AlphaMissense variant CSV files.",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="offline",
        help="Mode to fetch AlphaMissense data: 'online' or 'offline'.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Whether to overwrite existing AlphaMissense files.",
    )
    parser.add_argument(
        "--af_missense_tsv",
        type=str,
        required=False,
        default=AF_MISSENSE_AA_SUBSTITUTIONS_TSV,
        help="Path to AlphaMissense aa substitutions TSV file for offline mode.",
    )
    # parser.add_argument(
    #     "--pairwise_alignments_dir",
    #     type=str,
    #     default="/home/omkar/Projects/cardiac_desmosome/data/sequence_alignments/pairwise_alignments",
    #     help="Directory to save pairwise alignment files.",
    # )
    # parser.add_argument(
    #     "--odp_sequences_fasta",
    #     type=str,
    #     default="/home/omkar/Projects/cardiac_desmosome/data/sequences/odp_protein_sequences.fasta",
    #     help="FASTA file containing modeled protein sequences.",
    # )

    args = parser.parse_args()

    config_yaml = yaml.load(open(args.config_file, "r"), Loader=yaml.FullLoader)
    protein_uniprot_map = config_yaml["cardiac_odp_protein_uniprot_map"]
    uniprot_bases = [uid.split("-")[0] for uid in protein_uniprot_map.values()]

    af_missense_df_gen = fetch_af_missense_data(
        args.alpha_missense_dir,
        uniprot_bases,
        mode="offline",
        overwrite=args.overwrite,
        af_missense_tsv=args.af_missense_tsv,
    )

    export_af_missense_data(
        alpha_missense_dir=args.alpha_missense_dir,
        af_missense_df_gen=af_missense_df_gen,
        overwrite=args.overwrite,
    )

    # afm_fasta_dict = get_fasta_dict_for_af_missense(protein_uniprot_map)
    # odp_sequences = read_fasta(args.odp_sequences_fasta)
    # af_missense_dict = {}
    # for p_name, uniprot_id in protein_uniprot_map.items():

    #     print(f"\nProcessing {p_name} ({uniprot_id})...")

    #     modeled_seq = odp_sequences.get(protein_uniprot_map[p_name], None)
    #     if modeled_seq is None:
    #         warnings.warn(
    #             f"""
    #             Modeled sequence not found for {p_name} ({uniprot_id}). Got:
    #             {modeled_seq=}

    #             Skipping...
    #             """
    #         )
    #         continue

    ###########################################################################
    # Get AlphaMissense scores for variants in the protein
    ###########################################################################


        # afm_pairwise_alignment_file = os.path.join(
        #     args.pairwise_alignments_dir, f"{p_name}{AF_MISSENSE_PAIR_ALN_SUFFIX}.fasta"
        # )
        # uniprot_base = uniprot_id.split("-")[0]

        # af_missense_ref_seq = afm_fasta_dict.get(uniprot_base, None)

        # if af_missense_ref_seq is None:
        #     warnings.warn(
        #         f"""
        #         Reference sequence not found for {p_name}. Got:
        #         {af_missense_ref_seq=}

        #         Skipping...
        #         """
        #     )
        #     continue

        # won't align if sequences are identical
        # afm_psa_map = handle_pairwise_alignment(
        #     p_name=p_name,
        #     sseq=af_missense_ref_seq,
        #     qseq=modeled_seq,
        #     pairwise_alignment_file=afm_pairwise_alignment_file,
        #     ignore_warnings=True
        # )
