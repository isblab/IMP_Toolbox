import os
from string import Template
import xmltodict
import yaml
import pprint
import warnings
import argparse
import pandas as pd
from tqdm import tqdm
from io import StringIO
from datetime import datetime
import xml.etree.ElementTree as ET
from IMP_Toolbox.pre_processing.sequence.Sequence import FetchSequences
from IMP_Toolbox.utils_imp_toolbox.api_helpers import (
    request_session,
    request_result,
)
from IMP_Toolbox.utils_imp_toolbox.file_helpers import (
    read_fasta,
    read_json,
    write_json,
)
from IMP_Toolbox.utils_imp_toolbox.special_helpers import handle_pairwise_alignment
from IMP_Toolbox.utils_imp_toolbox.obj_helpers import fasta_str_to_dict
from IMP_Toolbox.pre_processing.mutations.utils_mutation import (
    split_missense_mutation,
    is_missense_mutation,
    get_ncbi_ref_seq,
)
from IMP_Toolbox.pre_processing.mutations.mutation_constants import (
    AF_MISSENSE_CSV_SUFFIX,
    AF_MISSENSE_PAIR_ALN_SUFFIX,
    CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE,
    CLINVAR_TEMPLATE_QUERY_DETAIL,
    CLINVAR_TEMPLATE_QUERY_ID,
    DATE_FORMAT,
    API_URLS,
)
from IMP_Toolbox.pre_processing.mutations.af_missense import (
    af_missense_df_to_dict,
    fetch_fasta_dict_for_af_missense,
    get_af_missense_data,
    get_fasta_dict_for_af_missense,
    fetch_af_missense_data,
)

def get_variant_ids_from_clinvar(
    gene_name:str,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
) -> list:
    """ Fetch variant IDs from ClinVar for a given gene name.

    TODO: Error handling is messed up here, fix it.

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

    print(response.url)

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

def fetch_variant_ids_from_clinvar(
    save_path: str,
    gene_name:str,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
    overwrite: bool = False,
) -> list:

    if os.path.exists(save_path) and not overwrite:
        variant_ids = read_json(save_path)
        return variant_ids

    variant_ids = get_variant_ids_from_clinvar(
        gene_name,
        api_url,
        api_parameters,
        ignore_error,
        max_retries,
    )

    write_json(save_path, variant_ids)

    return variant_ids

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

def fetch_variant_details_from_clinvar(
    save_path: str,
    variant_ids:list,
    api_url: str,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
    overwrite: bool = False,
) -> dict:

    if os.path.exists(save_path) and overwrite is False:
        print(f"Variant details file {save_path} already exists. Loading...")
        variant_info_dict = read_json(save_path)
        return variant_info_dict

    variant_info_dict = get_variant_details_from_clinvar(
        variant_ids,
        api_url,
        api_parameters,
        ignore_error,
        max_retries,
    )

    write_json(save_path, variant_info_dict)

    return variant_info_dict

def get_molecular_consequence_list(
    hgvs_list: dict,
    ncbi_ref_seq_id: str,
    missense_only: bool = True,
) -> list:
    """ Get the molecular consequence list from the HGVS list for a given
    RefSeq ID.

    Args:
        hgvs_list (dict):
            The HGVS list from the ClinVar variant information.
        ncbi_ref_seq_id (str):
            The RefSeq ID to match against.
        missense_only (bool, optional):
            Whether to filter for missense variants only. Defaults to True.

    Returns:
        list:
            List of molecular consequences.
    """

    def check_mc_type(mc: str) -> bool:
        if missense_only:
            return mc.get("@Type", "") == "missense variant"
        return True

    molecular_consequence_list = []

    if isinstance(hgvs_list.get("HGVS", []), dict):
        hgvs_list["HGVS"] = [hgvs_list["HGVS"]]

    for hgvs_entry in hgvs_list.get("HGVS", []):

        # we only want consequence at protein level
        if hgvs_entry.get("@Type", "") != "coding":
            continue

        # depending on the sequence, molecular consequence can differ
        # we only want those that match the RefSeq ID in the title
        # check if the sequence identifier matches
        nucleotide_express = hgvs_entry.get("NucleotideExpression", {})
        seq_accession = nucleotide_express.get("@sequenceAccession", "")

        # can potentially cause issues if the aa sequence changed in
        # newer versions of same base RefSeq ID
        # "@sequenceAccessionVersion" provides the full ID, but it can be
        # different from the one in the title
        if seq_accession != ncbi_ref_seq_id.split(".")[0]:
            continue

        mol_consequence = hgvs_entry.get("MolecularConsequence", {})

        if isinstance(mol_consequence, list):
            for mc in mol_consequence:
                if check_mc_type(mc):
                    molecular_consequence_list.extend([
                        mc.get("@Type", "") for mc in mol_consequence
                    ])
                    break

        elif isinstance(mol_consequence, dict):
            if check_mc_type(mol_consequence):
                molecular_consequence_list.append(
                    mol_consequence.get("@Type", "")
                )

    return list(set(molecular_consequence_list))

def extract_preferred_trait_names(
    trait_dict: dict,
    preferred_trait_names: list
) -> list:
    """ Extract preferred trait names from a trait dictionary.

    Args:
        trait_dict (dict): Trait dictionary from ClinVar variant information.
        preferred_trait_names (list): List to append preferred trait names to.

    Returns:
        list: Updated list of preferred trait names.
    """

    trait_names = trait_dict.get("Name", [])

    if not isinstance(trait_names, list):
        trait_names = [trait_names]

    for trait_name in trait_names:

        if not isinstance(trait_name, dict):
            continue

        if trait_name.get("ElementValue", {}).get("@Type", "") == "Preferred":
            preferred_trait_names.append(
                trait_name.get("ElementValue", {}).get("#text", "")
            )

    return preferred_trait_names

def get_variant_associated_traits(germline_classification: dict) -> list:
    """ Get the associated traits for a given germline classification.

    Args:
        germline_classification (dict):
            Germline classification from ClinVar variant information.

    Returns:
        list:
            List of associated traits.
    """

    traits = []
    trait_set = germline_classification.get("ConditionList", {}).get(
        "TraitSet", []
    )

    if isinstance(trait_set, dict):
        trait_set = [trait_set]

    for trait in trait_set:
        # trait should be a disease and should contribute to agg. classification
        if (
            trait.get("@Type", "") != "Disease"
            or trait.get(
                "@ContributesToAggregateClassification", "false"
            ) == "false"
        ):
            continue

        preferred_trait_names = []
        trait_ = trait.get("Trait", {})

        if isinstance(trait_, dict):
            preferred_trait_names = extract_preferred_trait_names(
                trait_, preferred_trait_names
            )

        elif isinstance(trait_, list):
            for _trait_ in trait_:
                preferred_trait_names = extract_preferred_trait_names(
                    _trait_, preferred_trait_names
                )

        traits.extend(preferred_trait_names)

    return list(set(traits))

def get_all_clinical_assertions(
    classified_record: dict,
    sort_by_date: bool=True,
    date_format:str=DATE_FORMAT,
) -> list:
    """ Get all clinical assertions from a classified record.

    Clinical assertion represents a single submission to ClinVar corresponding
    to the variant. Multiple submissions can be present for a single variant.
    We only look at those that contribute to the agg. clinical significance.

    if sort_by_date is True, the assertions are sorted by the date they were
    last evaluated. If the evaluation date is not available, the date last
    updated is used instead.

    Args:
        classified_record (dict):
            Classified record from ClinVar variant information.
        sort_by_date (bool, optional):
            Whether to sort the assertions by date. Defaults to True.
        date_format (str, optional):
            Date format for parsing dates. Defaults to DATE_FORMAT.

    Returns:
        list: List of clinical assertions.
    """

    clinical_assertions = []

    clinical_assertion_list = classified_record.get(
        "ClinicalAssertionList", {}
    ).get("ClinicalAssertion", [])

    if isinstance(clinical_assertion_list, dict):
        clinical_assertion_list = [clinical_assertion_list]

    for clinical_assertion in clinical_assertion_list:
        # we look at only those that contribute
        contributes_to_agg_class = clinical_assertion.get(
            "@ContributesToAggregateClassification", "false"
        )
        if contributes_to_agg_class == "false":
            continue

        last_updated_date = clinical_assertion.get(
            "@DateLastUpdated", ""
        )
        assertion_classification = clinical_assertion.get(
            "Classification", {}
        )
        asserted_germline_class = assertion_classification.get(
            "GermlineClassification", ""
        )

        assertion_comment = assertion_classification.get("Comment", "")

        if isinstance(assertion_comment, dict):
            assertion_comment = assertion_comment.get("#text", "")

        last_evaluated_date = assertion_classification.get(
            "@DateLastEvaluated", ""
        )

        last_updated_date = datetime.strptime(
            last_updated_date, date_format
        )

        try:
            last_evaluated_date = datetime.strptime(
                last_evaluated_date, date_format
            )
        except ValueError:
            last_evaluated_date = last_updated_date

        clinical_assertions.append({
            "last_updated_date": last_updated_date,
            "asserted_germline_class": asserted_germline_class,
            "assertion_comment": assertion_comment,
            "last_evaluated_date": last_evaluated_date,
        })

    if sort_by_date:
        clinical_assertions = sorted(
            clinical_assertions,
            key=lambda x: x["last_evaluated_date"],
            reverse=True
        )

    return clinical_assertions

def warn_sequence_not_found(
    p_name,
    uniprot_id,
):
    warnings.warn(
        f"""
        Sequence not found for {p_name} ({uniprot_id}).
        Skipping...
        """
    )

class VariantInfo:

    def __init__(
        self,
        variant_id: str,
        variant_info: dict,
    ):
        self.variant_id = variant_id
        self.variant_info = variant_info
        self.set_variant_archive()
        self.set_variant_name()
        self.set_variant_type()
        self.set_classified_record()
        self.set_hgvs_list()
        self.set_germline_classification()
        self.set_agg_significance()
        self.set_trait_set()

    def set_variant_archive(self):
        self.variant_archive = self.variant_info.get("VariationArchive", {})

    def set_variant_name(self):
        self.variant_name = self.variant_archive.get("@VariationName", "")

    def set_variant_type(self):
        self.variant_type = self.variant_archive.get("@VariationType", "")

    def set_classified_record(self) -> dict:
        self.classified_record = self.variant_archive.get("ClassifiedRecord", {})

    def set_hgvs_list(self) -> dict:
        self.hgvs_list = self.classified_record.get("SimpleAllele", {}).get("HGVSlist", {})

    def set_germline_classification(self) -> dict:
        self.germline_classification = self.classified_record.get(
            "Classifications", {}
        ).get("GermlineClassification", {})

    def set_agg_significance(self) -> str:
        ###################################################################
        # get aggregate significance
        # - pathogenic or likely pathogenic only
        # - depending on user input, include VUS
        ###################################################################
        self.agg_significance = self.germline_classification.get(
            "Description", ""
        )

    @staticmethod
    def get_ncbi_ref_seq_id(variant_name):
        # essentially points to which isoform the variant is described for
        # assumes format "NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)"
        # i.e. RefSeqID(GeneName):c.DNAChange (p.ProteinChange)
        ncbi_ref_seq_id = variant_name.split(":")[0].split("(")[0]
        return ncbi_ref_seq_id

    def get_molecular_consequence_list(
        self,
        ncbi_ref_seq_id: str,
        missense_only: bool = True,
    ) -> list:
        """ Get the molecular consequence list from the HGVS list for a given
        RefSeq ID.

        Args:
            hgvs_list (dict):
                The HGVS list from the ClinVar variant information.
            ncbi_ref_seq_id (str):
                The RefSeq ID to match against.
            missense_only (bool, optional):
                Whether to filter for missense variants only. Defaults to True.

        Returns:
            list:
                List of molecular consequences.
        """
        ###################################################################
        # get molecular consequence
        # - missense only
        ###################################################################
        def check_mc_type(mc: str) -> bool:
            if missense_only:
                return mc.get("@Type", "") == "missense variant"
            return True

        molecular_consequence_list = []

        if isinstance(self.hgvs_list.get("HGVS", []), dict):
            self.hgvs_list["HGVS"] = [self.hgvs_list["HGVS"]]

        for hgvs_entry in self.hgvs_list.get("HGVS", []):

            # we only want consequence at protein level
            if hgvs_entry.get("@Type", "") != "coding":
                continue

            # depending on the sequence, molecular consequence can differ
            # we only want those that match the RefSeq ID in the title
            # check if the sequence identifier matches
            nucleotide_express = hgvs_entry.get("NucleotideExpression", {})
            seq_accession = nucleotide_express.get("@sequenceAccession", "")

            # can potentially cause issues if the aa sequence changed in
            # newer versions of same base RefSeq ID
            # "@sequenceAccessionVersion" provides the full ID, but it can be
            # different from the one in the title
            if seq_accession != ncbi_ref_seq_id.split(".")[0]:
                continue

            mol_consequence = hgvs_entry.get("MolecularConsequence", {})

            if isinstance(mol_consequence, list):
                for mc in mol_consequence:
                    if check_mc_type(mc):
                        molecular_consequence_list.extend([
                            mc.get("@Type", "") for mc in mol_consequence
                        ])
                        break

            elif isinstance(mol_consequence, dict):
                if check_mc_type(mol_consequence):
                    molecular_consequence_list.append(
                        mol_consequence.get("@Type", "")
                    )

        return list(set(molecular_consequence_list))

    @staticmethod
    def get_mutation_descs(variant_name: str) -> list:
        mutation_descs = variant_name.split(":")[1].split(" ")
        return mutation_descs

    @staticmethod
    def get_ncbi_g_name(variant_name: str) -> str:
        ncbi_g_name = variant_name.split(":")[0].split("(")[1].replace(")", "")
        return ncbi_g_name

    @staticmethod
    def get_p_mutation(variant_name: str) -> str:
        ###################################################################
        # get protein mutation
        ###################################################################
        # assumes format "NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)"
        # i.e. RefSeqID(GeneName):c.DNAChange (p.ProteinChange)
        mutation_descs = VariantInfo.get_mutation_descs(variant_name)
        p_mutation = None
        for desc in mutation_descs:
            desc = desc.replace("(", "").replace(")", "")
            if desc.startswith("p."):
                p_mutation = desc.replace("p.", "")

        return p_mutation

    def set_trait_set(self):
        trait_set = self.germline_classification.get("ConditionList", {}).get(
            "TraitSet", []
        )
        if isinstance(trait_set, dict):
            self.trait_set = [self.trait_set]
        elif isinstance(trait_set, list):
            self.trait_set = trait_set

    def get_variant_associated_traits(self) -> list:
        """ Get the associated traits for a given germline classification.

        Args:
            germline_classification (dict):
                Germline classification from ClinVar variant information.

        Returns:
            list:
                List of associated traits.
        """

        traits = []
        # trait_set = self.germline_classification.get("ConditionList", {}).get(
        #     "TraitSet", []
        # )

        # if isinstance(self.trait_set, dict):
        #     self.trait_set = [self.trait_set]

        for trait in self.trait_set:
            # trait should be a disease and should contribute to agg. classification
            if (
                trait.get("@Type", "") != "Disease"
                or trait.get(
                    "@ContributesToAggregateClassification", "false"
                ) == "false"
            ):
                continue

            preferred_trait_names = []
            trait_ = trait.get("Trait", {})

            if isinstance(trait_, dict):
                preferred_trait_names = extract_preferred_trait_names(
                    trait_, preferred_trait_names
                )

            elif isinstance(trait_, list):
                for _trait_ in trait_:
                    preferred_trait_names = extract_preferred_trait_names(
                        _trait_, preferred_trait_names
                    )

            traits.extend(preferred_trait_names)

        return list(set(traits))

    def set_clinical_assertion_list(self):
        clinical_assertion_list = self.classified_record.get(
            "ClinicalAssertionList", {}
        ).get("ClinicalAssertion", [])

        if isinstance(clinical_assertion_list, dict):
            self.clinical_assertion_list = [clinical_assertion_list]
        elif isinstance(clinical_assertion_list, list):
            self.clinical_assertion_list = clinical_assertion_list

    def get_all_clinical_assertions(
        self,
        sort_by_date: bool=True,
        date_format:str=DATE_FORMAT,
    ) -> list:
        """ Get all clinical assertions from a classified record.

        Clinical assertion represents a single submission to ClinVar corresponding
        to the variant. Multiple submissions can be present for a single variant.
        We only look at those that contribute to the agg. clinical significance.

        if sort_by_date is True, the assertions are sorted by the date they were
        last evaluated. If the evaluation date is not available, the date last
        updated is used instead.

        Args:
            classified_record (dict):
                Classified record from ClinVar variant information.
            sort_by_date (bool, optional):
                Whether to sort the assertions by date. Defaults to True.
            date_format (str, optional):
                Date format for parsing dates. Defaults to DATE_FORMAT.

        Returns:
            list: List of clinical assertions.
        """

        clinical_assertions = []

        # clinical_assertion_list = classified_record.get(
        #     "ClinicalAssertionList", {}
        # ).get("ClinicalAssertion", [])

        # if isinstance(clinical_assertion_list, dict):
        #     clinical_assertion_list = [clinical_assertion_list]

        for clinical_assertion in self.clinical_assertion_list:
            # we look at only those that contribute
            contributes_to_agg_class = clinical_assertion.get(
                "@ContributesToAggregateClassification", "false"
            )
            if contributes_to_agg_class == "false":
                continue

            last_updated_date = clinical_assertion.get(
                "@DateLastUpdated", ""
            )
            assertion_classification = clinical_assertion.get(
                "Classification", {}
            )
            asserted_germline_class = assertion_classification.get(
                "GermlineClassification", ""
            )

            assertion_comment = assertion_classification.get("Comment", "")

            if isinstance(assertion_comment, dict):
                assertion_comment = assertion_comment.get("#text", "")

            last_evaluated_date = assertion_classification.get(
                "@DateLastEvaluated", ""
            )

            last_updated_date = datetime.strptime(
                last_updated_date, date_format
            )

            try:
                last_evaluated_date = datetime.strptime(
                    last_evaluated_date, date_format
                )
            except ValueError:
                last_evaluated_date = last_updated_date

            clinical_assertions.append({
                "last_updated_date": last_updated_date,
                "asserted_germline_class": asserted_germline_class,
                "assertion_comment": assertion_comment,
                "last_evaluated_date": last_evaluated_date,
            })

        if sort_by_date:
            clinical_assertions = sorted(
                clinical_assertions,
                key=lambda x: x["last_evaluated_date"],
                reverse=True
            )

        return clinical_assertions

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Fetch missense variants from ClinVar."
    )
    parser.add_argument(
        "--config_file",
        type=str,
        default="/home/omkar/Projects/cardiac_desmosome/input/config.yaml",
        # default=config_file,
        help="Path to the configuration YAML file.",
    )
    parser.add_argument(
        "--include_VUS",
        action="store_true",
        default=False,
        help="Include variants of uncertain significance (VUS).",
    )
    parser.add_argument(
        "--clinvar_output_dir",
        type=str,
        default="/home/omkar/Projects/cardiac_desmosome/data/literature_parsing/mutations/clinvar",
        # default=clinvar_output_dir,
        help="Directory to save ClinVar variant information.",
    )
    parser.add_argument(
        "--alpha_missense_dir",
        type=str,
        default="/home/omkar/Projects/cardiac_desmosome/data/literature_parsing/mutations/alpha_missense",
        # default=alpha_missense_dir,
        help="Directory to save AlphaMissense variant information.",
    )
    parser.add_argument(
        "--pairwise_alignments_dir",
        type=str,
        default="/home/omkar/Projects/cardiac_desmosome/data/sequence_alignments/pairwise_alignments",
        # default=pairwise_alignments_dir,
        help="Directory to save pairwise alignments.",
    )
    parser.add_argument(
        "--include_AF_missense",
        action="store_true",
        default=False,
        help="Include AlphaMissense pathogenicity scores.",
    )
    parser.add_argument(
        "--af_missense_mode",
        type=str,
        choices=["online", "offline"],
        default="offline",
        help="Mode to fetch AlphaMissense data.",
    )
    parser.add_argument(
        "--odp_sequences_fasta",
        type=str,
        default="/home/omkar/Projects/cardiac_desmosome/data/sequences/odp_protein_sequences.fasta",
        # default=odp_sequences_fasta,
        help="Fasta file containing modeled protein sequences.",
    )
    args = parser.parse_args()
    print("Doing something")

    config_yaml = yaml.load(open(args.config_file, "r"), Loader=yaml.FullLoader)
    protein_uniprot_map = config_yaml["cardiac_odp_protein_uniprot_map"]
    protein_gene_map = config_yaml["cardiac_odp_protein_gene_map"]
    odp_sequences = read_fasta(args.odp_sequences_fasta)

    if args.include_VUS:
        CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE.append("Uncertain significance")

    os.makedirs(args.clinvar_output_dir, exist_ok=True)

    df_rows = []
    fasta_dict = fetch_fasta_dict_for_af_missense(
        os.path.join(args.alpha_missense_dir, "af_missense_sequences.fasta"),
        protein_uniprot_map,
    )
    uniprot_bases = [uid.split("-")[0] for uid in protein_uniprot_map.values()]

    af_missense_dict = {}

    if args.include_AF_missense:
        _user = os.getlogin()
        af_missense_script_path = f"/home/{_user}/Projects/IMP_Toolbox/pre_processing/mutations/af_missense.py"
        af_missense_command = f"""
        python {af_missense_script_path} \\
            --config_file {os.path.abspath(args.config_file)} \\
            --alpha_missense_dir {os.path.abspath(args.alpha_missense_dir)} \\
            --mode {args.af_missense_mode} \\
        """
        os.system(af_missense_command)

    from collections import defaultdict
    pubmed_ids = defaultdict(set)

    for p_name, uniprot_id in protein_uniprot_map.items():

        uniprot_base = uniprot_id.split("-")[0]

        print(f"Processing {p_name}...")
        g_name = protein_gene_map[p_name]

        af_missense_dict = {}
        modeled_seq = odp_sequences.get(protein_uniprot_map[p_name], None)
        if modeled_seq is None:
            warn_sequence_not_found(p_name, uniprot_id)
            continue

        #######################################################################
        # Get AlphaMissense scores for variants in the protein
        #######################################################################
        if args.include_AF_missense:

            af_missense_file = os.path.join(
                args.alpha_missense_dir, f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
            )
            afm_pairwise_alignment_file = os.path.join(
                args.pairwise_alignments_dir, f"{p_name}{AF_MISSENSE_PAIR_ALN_SUFFIX}.fasta"
            )

            af_missense_ref_seq = fasta_dict.get(uniprot_base, None)

            if af_missense_ref_seq is None:
                warn_sequence_not_found(p_name, uniprot_base)
                continue

            # won't align if sequences are identical
            afm_psa_map = handle_pairwise_alignment(
                p_name=p_name,
                sseq=af_missense_ref_seq,
                qseq=modeled_seq,
                pairwise_alignment_file=afm_pairwise_alignment_file,
                ignore_warnings=True
            )
            # write_json(
            #     os.path.join(
            #         args.pairwise_alignments_dir,
            #         f"{p_name}_afm_pairwise_alignment_map.json"
            #     ),
            #     afm_psa_map
            # )

            if not os.path.exists(af_missense_file):
                warnings.warn(
                    f"AlphaMissense data file not found for {p_name}. Got: {af_missense_file}"
                )
                continue

            af_missense_df = pd.read_csv(af_missense_file)

            af_missense_dict = af_missense_df_to_dict(
                p_name,
                af_missense_df,
                af_missense_dict=af_missense_dict,
                afm_psa_map=afm_psa_map, # may not be appropriate in all cases
            )
            # write_json(
            #     os.path.join(
            #         args.alpha_missense_dir,
            #         f"{uniprot_base}_af_missense_dict.json"
            #     ),
            #     af_missense_dict
            # )

        #######################################################################
        # Get gene specific variants from ClinVar
        #######################################################################
        clinvar_variants_file = os.path.join(
            args.clinvar_output_dir, f"{g_name}_clinvar_variants1.json"
        )

        print(f"Fetching ClinVar variants for {g_name}...")

        api_parameters = CLINVAR_TEMPLATE_QUERY_ID.copy()
        api_parameters["term"] = Template(api_parameters["term"]).substitute(
            gene_name=g_name
        )

        variant_ids = fetch_variant_ids_from_clinvar(
            save_path=os.path.join(
                args.clinvar_output_dir, f"{g_name}_clinvar_variant_ids.json"
            ),
            gene_name=g_name,
            api_url=API_URLS["ncbi_esearch"],
            api_parameters=api_parameters,
            ignore_error=False,
            max_retries=3,
            overwrite=False,
        )

        if len(variant_ids) == 0:
            print(f"No variants found for {g_name}")
            continue

        clinvar_variants = fetch_variant_details_from_clinvar(
            save_path=clinvar_variants_file,
            variant_ids=variant_ids,
            api_url=API_URLS["ncbi_efetch"],
            api_parameters=CLINVAR_TEMPLATE_QUERY_DETAIL,
            ignore_error=False,
            max_retries=3,
            overwrite=False,
        )

        ncbi_ref_seq_ids = set()

        for variant_id, variant_info in clinvar_variants.items():

            vi = VariantInfo(variant_id, variant_info)

            if vi.agg_significance not in CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE:
                continue

            ncbi_ref_seq_id = VariantInfo.get_ncbi_ref_seq_id(vi.variant_name)

            molecular_consequence_list = vi.get_molecular_consequence_list(
                ncbi_ref_seq_id, missense_only=True
            )
            # If there are no missense variants, skip
            if len(molecular_consequence_list) == 0:
                continue

            mutation_descs = VariantInfo.get_mutation_descs(vi.variant_name)
            ncbi_g_name = VariantInfo.get_ncbi_g_name(vi.variant_name)

            # there is one incorrect entry in DSC2 variants
            if g_name != ncbi_g_name:
                warnings.warn(
                    f"Gene mismatch for {variant_id}: {g_name} vs {ncbi_g_name}"
                )
                continue

            p_mutation = VariantInfo.get_p_mutation(vi.variant_name)
            #NOTE: the following will miss missense mutations where
            # more than one residue is mutated
            # e.g. p.Gly1094_His1095delinsValAsn in DSG2
            if (
                p_mutation is None
                or not is_missense_mutation(
                    p_mutation, allow_truncation=False
                )
            ):
                warnings.warn(
                    f"""
                    Protein mutation not found or not missense for {variant_id}
                    {p_mutation=}
                    Molecular consequence might be "missense", but more than
                    one residue might be mutated.
                    Skipping...
                    """
                )
                continue

            traits = vi.get_variant_associated_traits()

            if len(traits) == 0:
                traits = ["not provided"]

            # Look through all the assertions that contributed to aggregate
            # clinical significance
            clinical_assertions = vi.get_all_clinical_assertions(
                sort_by_date=True, date_format=DATE_FORMAT
            )

            last_significance = [""] # these are actually all assertions
            last_assertion_comment = [""]

            if len(clinical_assertions) > 0:

                last_significance = [
                    clinical_assertions[i]["asserted_germline_class"]
                    for i in range(len(clinical_assertions))
                ]
                last_assertion_comment = [
                    clinical_assertions[i]["assertion_comment"]
                    for i in range(len(clinical_assertions))
                ]

            last_assertion_comment_zip = [
                f"{s} :- {c}" for s, c in zip(last_significance, last_assertion_comment)
            ]

            # if agg. significance is conflicting, last assertion should be in allowed
            if (
                agg_significance == "Conflicting classifications of pathogenicity"
                and last_significance[0] not in CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE
            ):
                continue


        for variant_id, variant_info in clinvar_variants.items():

            variant_archive = variant_info.get("VariationArchive", {})
            variant_name = variant_archive.get("@VariationName", "")
            variant_type = variant_archive.get("@VariationType", "")
            classified_record = variant_archive.get("ClassifiedRecord", {})

            ###################################################################
            # get aggregate significance
            # - pathogenic or likely pathogenic only
            # - depending on user input, include VUS
            ###################################################################

            germline_classification = classified_record.get(
                "Classifications", {}
            ).get("GermlineClassification", {})

            agg_significance = germline_classification.get(
                "Description", ""
            )
            if agg_significance not in CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE:
                continue

            ###################################################################
            # get molecular consequence
            # - missense only
            ###################################################################
            hgvs_list = classified_record.get(
                "SimpleAllele", {}
            ).get("HGVSlist", {})

            # essentially points to which isoform the variant is described for
            # assumes format "NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)"
            # i.e. RefSeqID(GeneName):c.DNAChange (p.ProteinChange)
            ncbi_ref_seq_id = variant_name.split(":")[0].split("(")[0]

            molecular_consequence_list = get_molecular_consequence_list(
                hgvs_list, ncbi_ref_seq_id, missense_only=True
            )

            # If there are no missense variants, skip
            if len(molecular_consequence_list) == 0:
                continue

            ###################################################################
            # get protein mutation
            ###################################################################
            # assumes format "NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)"
            # i.e. RefSeqID(GeneName):c.DNAChange (p.ProteinChange)
            mutation_descs = variant_name.split(":")[1].split(" ")
            ncbi_g_name = variant_name.split(":")[0].split("(")[1].replace(")", "")

            # there is one incorrect entry in DSC2 variants
            if g_name != ncbi_g_name:
                warnings.warn(
                    f"Gene mismatch for {variant_id}: {g_name} vs {ncbi_g_name}"
                )
                continue

            # get mutation
            p_mutation = None
            for desc in mutation_descs:

                desc = desc.replace("(", "").replace(")", "")

                if desc.startswith("p."):
                    p_mutation = desc.replace("p.", "")

            #NOTE: the following will miss missense mutations where
            # more than one residue is mutated
            # e.g. p.Gly1094_His1095delinsValAsn in DSG2
            if (
                p_mutation is None
                or not is_missense_mutation(
                    p_mutation, allow_truncation=False
                )
            ):
                warnings.warn(
                    f"""
                    Protein mutation not found or not missense for {variant_id}
                    {p_mutation=}
                    Molecular consequence might be "missense", but more than
                    one residue might be mutated.
                    Skipping...
                    """
                )
                continue

            afm_patho_score = ""
            afm_pathogenicity = ""

            if (
                args.include_AF_missense
                and p_mutation in af_missense_dict.get(p_name, {})
            ):

                afm_patho_score = af_missense_dict[p_name][p_mutation][
                    "patho_score"
                ]

                afm_pathogenicity = af_missense_dict[p_name][p_mutation][
                    "v_pathogenicity"
                ]

            ###################################################################
            # get disease or trait information
            ###################################################################
            traits = get_variant_associated_traits(germline_classification)

            if len(traits) == 0:
                traits = ["not provided"]

            # Look through all the assertions that contributed to aggregate
            # clinical significance
            clinical_assertions = get_all_clinical_assertions(
                classified_record, sort_by_date=True, date_format=DATE_FORMAT
            )

            # last_significance = ""
            # last_assertion_comment = ""

            # if len(clinical_assertions) > 0:

            #     last_significance = clinical_assertions[0][
            #         "asserted_germline_class"
            #     ]
            #     last_assertion_comment = clinical_assertions[0][
            #         "assertion_comment"
            #     ]

            # # if agg. significance is conflicting, last assertion should be in allowed
            # if (
            #     agg_significance == "Conflicting classifications of pathogenicity"
            #     and last_significance not in CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE
            # ):
            #     continue

            last_significance = [""]
            last_assertion_comment = [""]

            if len(clinical_assertions) > 0:

                last_significance = [clinical_assertions[i]["asserted_germline_class"] for i in range(len(clinical_assertions))]
                last_assertion_comment = [clinical_assertions[i]["assertion_comment"] for i in range(len(clinical_assertions))]

            import re
            for comment in last_assertion_comment:
                # search all text in brackets
                pmids = re.findall(r"\(.*?\)", comment)
                if pmids is None:
                    continue
                # pritn all group matches
                for pmid_str in pmids:
                    if "21177847" in pmid_str:
                        print(pmids)

                if all(["PMID" not in pmid_str for pmid_str in pmids]):
                    continue
                if pmids:
                    # pmid_list = pmids.group(0).replace("(PMID: ","").replace(")","").split(",")
                    pmid_list = [
                        pmid_str.replace("PMID:", "").replace(")", "").replace("(", "").replace("; Invitae", "").split(",")
                        for pmid_str in pmids if "PMID" in pmid_str
                    ]
                    pmid_list = [pmid for sublist in pmid_list for pmid in sublist]
                    pubmed_ids[g_name].update([pmid.strip() for pmid in pmid_list])

            last_assertion_comment_zip = [
                f"{s} :- {c}" for s, c in zip(last_significance, last_assertion_comment)
            ]

            # if agg. significance is conflicting, last assertion should be in allowed
            if (
                agg_significance == "Conflicting classifications of pathogenicity"
                and last_significance[0] not in CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE
            ):
                continue

            ###################################################################
            # Adjust residue number based on alignment with modeled sequence
            ###################################################################
            ref_seq_file = os.path.join(
                args.clinvar_output_dir, f"{g_name}_{ncbi_ref_seq_id}_nuccore.json"
            )

            clinvar_pairwise_alignment_file = os.path.join(
                args.pairwise_alignments_dir,
                f"{p_name}_clinvar_vs_modeled.fasta"
            )

            p_sequence = get_ncbi_ref_seq(ncbi_ref_seq_id, ref_seq_file)

            clinvar_psa_map = handle_pairwise_alignment(
                p_name=p_name,
                sseq=p_sequence,
                qseq=modeled_seq,
                pairwise_alignment_file=clinvar_pairwise_alignment_file,
                ignore_warnings=True
            )

            wt_aa, res_num, mut_aa = split_missense_mutation(
                p_mutation, return_type="all"
            )

            res_num = int(res_num)
            res_num = clinvar_psa_map.get(res_num, res_num)
            p_mutation = f"{wt_aa}{res_num}{mut_aa}"

            row_dict = {
                "variant_id": variant_id,
                "protein": p_name,
                "gene": g_name,
                "uniprot_id": uniprot_id,
                "ncbi_ref_seq_id": ncbi_ref_seq_id,
                "p_mutation": p_mutation,
                "molecular_consequence": "\n".join(molecular_consequence_list),
                "trait_or_effect": "\n".join(traits),
                "clinical_significance": agg_significance,
                "last_significance": "\n\n".join(last_significance),
                "last_assertion_comment": "\n\n".join(last_assertion_comment_zip),
                "variant_type": variant_type,
                "afm_patho_score": afm_patho_score,
                "afm_pathogenicity": afm_pathogenicity,
            }

            df_rows.append(row_dict)
            ncbi_ref_seq_ids.add(ncbi_ref_seq_id)

        if len(ncbi_ref_seq_ids) == 0:
            print(f"No missense variants found for {g_name}")
            continue

        if len(ncbi_ref_seq_ids) != 1:
            raise NotImplementedError(
                f"Multiple RefSeq IDs found for {g_name}: {ncbi_ref_seq_ids}."
            )

    # ###########################################################################
    # # Monkey patching
    # ###########################################################################

    # df = pd.DataFrame(df_rows)

    # del df["molecular_consequence"]
    # del df["protein"]
    # # del df["last_significance"]
    # del df["variant_type"]

    # for idx, row in df.iterrows():
    #     trait_or_effects = row["trait_or_effect"].split("\n")
    #     trait_or_effects = [
    #         t for t in trait_or_effects if t not in [
    #             "not provided",
    #             "",
    #             "Cardiovascular phenotype",
    #             "not specified",
    #         ]
    #     ]

    #     if any(
    #         "Arrhythmogenic right ventricular dysplasia" in t
    #         for t in trait_or_effects
    #     ):

    #         if "Cardiomyopathy" in trait_or_effects:
    #             trait_or_effects.remove("Cardiomyopathy")

    #         if "Arrhythmogenic right ventricular cardiomyopathy" in trait_or_effects:
    #             trait_or_effects.remove(
    #                 "Arrhythmogenic right ventricular cardiomyopathy"
    #             )

    #         trait_or_effects = [
    #             t.replace(
    #                 "Arrhythmogenic right ventricular dysplasia",
    #                 "ARVD"
    #             ) for t in trait_or_effects
    #         ]

    #     if len(trait_or_effects) == 0:
    #         # look in comments for ARVC or similar terms
    #         comment = row["last_assertion_comment"].lower()
    #         terms_to_look_for = [
    #             "arvc",
    #             "arvd",
    #             "arrhythmia",
    #             "cardiomyopathy",
    #         ]

    #         if any(term in comment for term in terms_to_look_for):
    #             trait_or_effects = [
    #                 term.upper() for term in terms_to_look_for
    #                 if term in comment
    #             ]

    #         else:
    #             trait_or_effects = [""]

    #     df.at[idx, "trait_or_effect"] = "\n".join(trait_or_effects)

    #     if row["clinical_significance"] == "Conflicting classifications of pathogenicity":
    #         df.at[idx, "clinical_significance"] = (
    #             row["last_significance"]
    #         )

    # # sort by mutated residue number for each protein
    # df["residue_number"] = df["p_mutation"].apply(
    #     split_missense_mutation, return_type="res_num"
    # )
    # df = df.sort_values(by=["gene", "residue_number", "p_mutation"])
    # df = df.drop_duplicates(
    #     subset=["gene", "p_mutation"]
    # ).reset_index(drop=True)

    # df = df.drop(columns=["residue_number"])

    # del df["last_significance"]

    # # rename columns
    # df = df.rename(columns={
    #     "gene": "Gene",
    #     "uniprot_id": "Uniprot ID",
    #     "ncbi_ref_seq_id": "NCBI RefSeq ID",
    #     "p_mutation": "Mutation",
    #     "trait_or_effect": "Disease association",
    #     "clinical_significance": "ClinVar clinical significance",
    #     "last_assertion_comment": "All submission comments",
    #     # "last_significance": "Most recent clinical significance",
    #     "variant_id": "ClinVar Variant ID",
    #     # "molecular_consequence": "Molecular consequence",
    #     "afm_patho_score": "AlphaMissense score",
    #     "afm_pathogenicity": "AlphaMissense pathogenicity",
    # })

    # # if not args.include_AF_missense:
    # #     df = df.drop(columns=[
    # #         "AlphaMissense score", "AlphaMissense pathogenicity"
    # #     ])

    # print(df.head())

    # out_name = "clinvar_missense_variants"
    # if args.include_VUS:
    #     out_name += "_with_VUS"

    # df_file = os.path.join(
    #     args.clinvar_output_dir, f"{out_name}.xlsx"
    # )
    # df.to_excel(df_file, index=False)

    # print(f"ClinVar missense variants saved to {df_file}")

    # pubmed_ids = {
    #     g: [int(pmid) for pmid in pubmed_ids[g] if pmid.isdigit()]
    #     for g in pubmed_ids.keys()
    # }
    # for g_name, pmids in pubmed_ids.items():
    #     print(g_name)
    #     for pmid in pmids:
    #         print(pmid)