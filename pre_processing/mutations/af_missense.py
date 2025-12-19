import os
import yaml
import warnings
import argparse
import pandas as pd
from io import StringIO
from pprint import pprint
from IMP_Toolbox.pre_processing.sequence.Sequence import FetchSequences
# from cardiac_desmosome.utils.where_is_it import WhereIsIt
# from cardiac_desmosome.constants.input_constants import ConfigYaml
from IMP_Toolbox.utils_imp_toolbox.file_helpers import read_fasta
from IMP_Toolbox.utils_imp_toolbox.special_helpers import handle_pairwise_alignment
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import fasta_str_to_dict
from IMP_Toolbox.pre_processing.mutations.utils_mutation import (
    split_missense_mutation,
    get_af_missense_data,
)
from IMP_Toolbox.pre_processing.mutations.mutation_constants import (
    AMINO_ACID_MAP,
    AFM_PATHOGENICITY,
    API_URLS,
)


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

if __name__ == "__main__":

    # where_is_it = WhereIsIt()

    # config_file = where_is_it.config_file
    # odp_sequences_fasta = where_is_it.sequences.odp_sequences_fasta

    # # config_file = "/home/omkar/Projects/wnt2/config.yaml"
    # # odp_sequences_fasta = "/home/omkar/Projects/wnt2/sequences.fasta"

    # pairwise_alignments_dir = where_is_it.alignments.pairwise_alignments_dir
    # alpha_missense_dir = where_is_it.literature_parsing.alpha_missense_dir

    # config_yaml = yaml.load(open(config_file, "r"), Loader=yaml.FullLoader)

    # # protein_uniprot_map = config_yaml["protein_uniprot_map"]
    # protein_uniprot_map = config_yaml[ConfigYaml.cardiac_p_to_u]

    parser = argparse.ArgumentParser(
        description="Fetch and save UniProt variants for cardiac ODP proteins."
    )

    parser.add_argument(
        "--config_file",
        type=str,
        # default=WhereIsIt().config_file,
        help="Path to configuration YAML file.",
    )

    parser.add_argument(
        "--alpha_missense_dir",
        type=str,
        # default=alpha_missense_dir,
        help="Directory to save AlphaMissense variant CSV files.",
    )
    parser.add_argument(
        "--pairwise_alignments_dir",
        type=str,
        # default=pairwise_alignments_dir,
        help="Directory to save pairwise alignment files.",
    )
    parser.add_argument(
        "--odp_sequences_fasta",
        type=str,
        # default=odp_sequences_fasta,
        help="FASTA file containing modeled protein sequences.",
    )

    args = parser.parse_args()

    config_yaml = yaml.load(open(args.config_file, "r"), Loader=yaml.FullLoader)
    protein_uniprot_map = config_yaml["protein_uniprot_map"]

    odp_sequences = read_fasta(args.odp_sequences_fasta)
    # For AF-missense, the reference sequence is the Uniprot sequence
    # for non-isoform-specific uniprot ids
    uniprot_ids = list(protein_uniprot_map.values())
    uniprot_ids = [uid.split("-")[0] for uid in uniprot_ids]
    fetchit = FetchSequences(uniprot_ids=uniprot_ids)
    fasta_str = fetchit.query_uniprot_api_for_sequences()
    fasta_str = fetchit.only_uniprot_id_as_name(fasta_str)
    fasta_dict = fasta_str_to_dict(fasta_str)

    af_missense_dict = {}
    for p_name, uniprot_id in protein_uniprot_map.items():

        print(f"\nProcessing {p_name} ({uniprot_id})...")

        modeled_seq = odp_sequences.get(protein_uniprot_map[p_name], None)
        if modeled_seq is None:
            warnings.warn(
                f"""
                Modeled sequence not found for {p_name} ({uniprot_id}). Got:
                {modeled_seq=}

                Skipping...
                """
            )
            continue

    ###########################################################################
    # Get AlphaMissense scores for variants in the protein
    ###########################################################################

        af_missense_file = os.path.join(
            args.alpha_missense_dir, f"{p_name}_alpha_missense_variants.csv"
        )
        afm_pairwise_alignment_file = os.path.join(
            args.pairwise_alignments_dir, f"{p_name}_afm_vs_modeled.fasta"
        )

        af_missense_ref_seq = fasta_dict.get(uniprot_id.split("-")[0], None)

        if af_missense_ref_seq is None:
            warnings.warn(
                f"""
                Reference sequence not found for {p_name}. Got:
                {af_missense_ref_seq=}

                Skipping...
                """
            )
            continue

        # won't align if sequences are identical
        afm_psa_map = handle_pairwise_alignment(
            p_name=p_name,
            sseq=af_missense_ref_seq,
            qseq=modeled_seq,
            pairwise_alignment_file=afm_pairwise_alignment_file,
            ignore_warnings=True
        )

        if os.path.exists(af_missense_file):
            print(f"File {af_missense_file} already exists. Skipping...")
            af_missense_df = pd.read_csv(af_missense_file)

        else:
            af_missense_csv = get_af_missense_data(
                uniprot_id,
                api_url=API_URLS["af_missense_csv"],
                api_parameters={"uniprot_id": uniprot_id.split("-")[0]},
                return_type="csv",
                ignore_error=False
            )
            af_missense_df = pd.read_csv(StringIO(af_missense_csv))
            af_missense_df.to_csv(af_missense_file, index=False)
            print(f"Saved AlphaMissense variants to {af_missense_file}")
