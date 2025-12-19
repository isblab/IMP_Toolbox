import os
import yaml
import warnings
import argparse
import pandas as pd
from io import StringIO
from pprint import pprint
from IMP_Toolbox.pre_processing.sequence.Sequence import FetchSequences
from IMP_Toolbox.utils_imp_toolbox.file_helpers import read_fasta
from IMP_Toolbox.utils_imp_toolbox.special_helpers import handle_pairwise_alignment
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

def get_af_missense_data_offline(
    uniprot_ids: list,
):
    try:
        import duckdb
    except ImportError:
        raise ImportError(
            "duckdb is required for offline AlphaMissense data retrieval. "
            "Please install it via"
            "'pip install duckdb' or 'conda install -c conda-forge duckdb'."
        )

    assert os.path.exists(AF_MISSENSE_AA_SUBSTITUTIONS_TSV), (
        f"AlphaMissense aa substitutions TSV file not found at "
        f"{AF_MISSENSE_AA_SUBSTITUTIONS_TSV}."
    )

    uniprot_ids = list(set(uniprot_ids))
    uids_df = pd.DataFrame({"uniprot_id": uniprot_ids})
    query = f"""
    SELECT t.*
    FROM '{AF_MISSENSE_AA_SUBSTITUTIONS_TSV}' t
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
                af_missense_csv = get_af_missense_data(
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

    for uniprot_base in uniprot_bases:
        af_missense_file = os.path.join(
            alpha_missense_dir, f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
        )

        if os.path.exists(af_missense_file):
            print(f"File {af_missense_file} already exists. Skipping...")
        else:
            remainder_bases.append(uniprot_base)

    remainder_bases = uniprot_bases if overwrite else remainder_bases

    return remainder_bases

def af_missense_df_to_dict(
    p_name: str,
    af_missense_df: pd.DataFrame,
    af_missense_dict: dict = {},
    afm_psa_map: dict = {},
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

        if any(aa not in AMINO_ACID_MAP for aa in [wt_aa, mut_aa]):
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

        res_num_mapped = afm_psa_map.get(res_num_int, res_num_int)
        p_variant_key = f"{wt_aa}{res_num_mapped}{mut_aa}"

        af_missense_dict[p_name][p_variant_key] = {
            "patho_score": patho_score,
            "v_pathogenicity": AFM_PATHOGENICITY.get(v_pathogenicity),
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "res_num": res_num_mapped,
        }

    return af_missense_dict

def get_fasta_dict_for_af_missense(protein_uniprot_map):

    # For AF-missense, the reference sequence is the Uniprot sequence
    # for non-isoform-specific uniprot ids
    uniprot_ids = list(protein_uniprot_map.values())
    uniprot_ids = [uid.split("-")[0] for uid in uniprot_ids]
    fetchit = FetchSequences(uniprot_ids=uniprot_ids)
    fasta_str = fetchit.query_uniprot_api_for_sequences()
    fasta_str = fetchit.only_uniprot_id_as_name(fasta_str)
    fasta_dict = fasta_str_to_dict(fasta_str)

    return fasta_dict

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
    )

    for af_missense_df, uniprot_base in af_missense_df_gen:

        if af_missense_df.empty:
            print(f"No AlphaMissense data found for {uniprot_base}. Skipping...")
            continue

        af_missense_file = os.path.join(
            args.alpha_missense_dir, f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
        )

        if os.path.exists(af_missense_file) and not args.overwrite:
            print(f"File {af_missense_file} already exists. Skipping save...")
            continue

        af_missense_df.to_csv(af_missense_file, index=False)

        print(f"Saved AlphaMissense variants to {af_missense_file}")

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
