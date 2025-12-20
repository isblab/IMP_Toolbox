import os
from string import Template
from typing import Any
import xmltodict
import yaml
import pprint
import warnings
import argparse
import pandas as pd
from tqdm import tqdm
from datetime import datetime
import xml.etree.ElementTree as ET
from IMP_Toolbox.utils_imp_toolbox.api_helpers import (
    request_session,
    request_result,
)
from IMP_Toolbox.utils_imp_toolbox.file_helpers import (
    read_fasta,
    read_json,
    write_json,
)
from IMP_Toolbox.utils_imp_toolbox.special_helpers import (
    get_mapped_residue,
    handle_pairwise_alignment,
)
from IMP_Toolbox.pre_processing.mutations.utils_mutation import (
    split_missense_mutation,
    is_missense_mutation,
    get_ncbi_ref_seq,
)
from IMP_Toolbox.pre_processing.mutations.af_missense import (
    af_missense_df_to_dict,
    fetch_fasta_dict_for_af_missense,
    fetch_af_missense_data,
)
from IMP_Toolbox.pre_processing.mutations.mutation_constants import (
    AF_MISSENSE_CSV_SUFFIX,
    AF_MISSENSE_PAIR_ALN_SUFFIX,
    AF_MISSENSE_AA_SUBSTITUTIONS_TSV,
    CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE,
    CLINVAR_DF_COLUMNS,
    CLINVAR_TEMPLATE_QUERY_DETAIL,
    CLINVAR_TEMPLATE_QUERY_ID,
    DATE_FORMAT,
    API_URLS,
)

def get_variant_ids_from_clinvar(
    gene_name:str,
    api_parameters: dict,
    ignore_error: bool = False,
    max_retries: int = 3,
) -> list:
    """ Fetch variant IDs from ClinVar for a given gene name.

    TODO: Error handling is messed up here, fix it.

    Args:

        gene_name (str):
            The gene name to search for.

        api_parameters (dict):
            API parameters for the request.

        ignore_error (bool, optional):
            Whether to ignore errors. Defaults to False.

        max_retries (int, optional):
            Maximum number of retries for the request. Defaults to 3.

    Returns:

        list:
            List of ClinVar variant IDs.
    """

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(API_URLS["ncbi_esearch"], params=api_parameters)
    result = request_result(response, gene_name)

    if result is None:
        # try without returnmode
        api_parameters.pop("retmode")
        response = req_sess.get(API_URLS["ncbi_esearch"], params=api_parameters)

        if response.status_code == 200:
            result = response.text

        elif ignore_error:
            # print(response.url)
            return []

        else:
            # print(response.url)
            raise ValueError(
                f"""Error fetching variant ids for {gene_name}
                {response.status_code}: {response.text}
                """
            )

    # print(response.url)
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
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
) -> dict:
    """ Fetch variant details from ClinVar for a list of variant IDs.
    Note: If you parallelize it, make sure to not exceed API's rate limits.
    We get info in XML format from the API, which we convert and store as JSON.

    Args:

        variant_ids (list):
            List of ClinVar variant IDs.

        api_parameters (dict):
            API parameters for the request.

        ignore_error (bool, optional):
            Whether to ignore errors. Defaults to False.

        max_retries (int, optional):
            Maximum number of retries for the request. Defaults to 3.

    Returns:

        dict:
            Dictionary with variant IDs as keys and their details as values.
    """

    variant_id_batches = [
        variant_ids[i:i + 100] for i in range(0, len(variant_ids), 100)
    ]

    variant_info_dict = {}

    for idx, batch in enumerate(tqdm(variant_id_batches)):

        api_parameters["id"] = ",".join(batch)

        req_sess = request_session(max_retries=max_retries)
        response = req_sess.get(API_URLS["ncbi_efetch"], params=api_parameters)

        if response.status_code == 200:
            result = response.text
        else:
            if ignore_error:
                return {}
            else:
                raise ValueError(f"Error fetching variant details for {batch}")

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

def fetch_clinvar_variant_data(
    save_path: str,
    input_data: str | list,
    api_parameters: dict = {},
    ignore_error: bool = False,
    max_retries: int = 3,
    overwrite: bool = False,
) -> dict | list:
    """ Fetch ClinVar variant data for given gene name or variant IDs.

    Uses `get_variant_ids_from_clinvar` or `get_variant_details_from_clinvar`.

    Args:

        save_path (str):
            Path to save the fetched variant data.

        input_data (str | list):
            Gene name (str) or list of variant IDs (list).

        api_parameters (dict, optional):
            API parameters for the request. Defaults to {}.

        ignore_error (bool, optional):
            Whether to ignore errors. Defaults to False.

        max_retries (int, optional):
            Maximum number of retries for the request. Defaults to 3.

        overwrite (bool, optional):
            Whether to overwrite existing data at save_path. Defaults to False.

    Returns:
        dict | list:
            Fetched variant data.
    """

    data_type_dict = {str: "ids", list: "details"}

    try:
        data_type = data_type_dict[type(input_data)]

    except KeyError:
        raise ValueError(
            f"Invalid input data type: {type(input_data)}. "
            "Expected str (gene name) or list (variant ids)."
        )

    if os.path.exists(save_path) and overwrite is False:
        variant_data = read_json(save_path)
        return variant_data

    get_func = eval(f"get_variant_{data_type}_from_clinvar")

    variant_data = get_func(
        input_data,
        api_parameters=api_parameters,
        ignore_error=ignore_error,
        max_retries=max_retries,
    )

    write_json(save_path, variant_data)

    return variant_data

class VariantInfo:

    setter_dict: dict = {dict: lambda x: [x], list: lambda x: x}

    def __init__(
        self,
        p_name: str,
        g_name: str,
        variant_id: str,
        variant_info: dict,
    ):
        self.p_name = p_name
        self.g_name = g_name
        self.variant_id = variant_id
        self.variant_archive = variant_info.get("VariationArchive", {})
        self.set_variant_name()
        self.set_variant_type()
        self.set_classified_record()
        self.set_hgvs_list()
        self.set_germline_classification()
        self.set_agg_significance()
        self.set_trait_set()
        self.set_clinical_assertion_list()

        self.ncbi_ref_seq_id = self.get_ncbi_ref_seq_id()
        self.molecular_consequence_list = self.get_molecular_consequence_list()
        self.mutation_descs = self.get_mutation_descs()
        self.p_mutation = self.get_p_mutation()
        self.ncbi_g_name = self.get_ncbi_g_name()
        self.traits = self.get_variant_associated_traits()
        self.clinical_assertions = self.get_all_clinical_assertions(
            sort_by_date=True,
            date_format=DATE_FORMAT,
        )
        self.all_significances = self.get_all_significances()
        self.all_assertion_comments = self.get_all_assertion_comments()

    def set_variant_name(self):
        """ Set the variant name from the variant archive.
        The varient name of our relevance is in the format:
        NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)
        RefSeqID(GeneName):c.DNAChange (p.ProteinChange)
        """
        self.variant_name = self.variant_archive.get("@VariationName", "")

    def set_variant_type(self):
        """ Set the variant type from the variant archive.
        The variant type of our relevance is "single nucleotide variant"
        """
        self.variant_type = self.variant_archive.get("@VariationType", "")

    def set_classified_record(self) -> dict:
        """ Set the classified record from the variant archive.
        This has the germline classification information.
        """
        self.classified_record = self.variant_archive.get("ClassifiedRecord", {})

    def set_hgvs_list(self) -> dict:
        """ Set the HGVS list from the classified record.
        The HGVS list contains the molecular consequence information.
        """
        self.hgvs_list = self.classified_record.get("SimpleAllele", {}).get("HGVSlist", {})

    def set_germline_classification(self) -> dict:
        """ Set the germline classification from the classified record.
        This has the aggregate clinical significance information.
        """
        self.germline_classification = self.classified_record.get(
            "Classifications", {}
        ).get("GermlineClassification", {})

    def set_agg_significance(self) -> str:
        """ Set the aggregate clinical significance.
        This is the overall clinical significance for the variant provided by
        ClinVar based on all contributing submissions.
        Of our relevance are "Pathogenic" and "Likely pathogenic" and
        potentially "VUS" depending on user input.
        """
        self.agg_significance = self.germline_classification.get(
            "Description", ""
        )

    def get_ncbi_ref_seq_id(self, ignore_warnings: bool = True) -> str | None:
        """ Get the NCBI RefSeq ID from the variant name.
        Essentially points to which isoform the variant is described for.
        assumes format "NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)"
        i.e. RefSeqID(GeneName):c.DNAChange (p.ProteinChange)

        Args:

            ignore_warnings (bool, optional):
                Whether to ignore warnings. Defaults to True.

        Returns:

            str | None:
                The NCBI RefSeq ID.
        """

        ncbi_ref_seq_id = None
        ncbi_ref_seq = self.variant_name.split(":")

        if self.g_name in ncbi_ref_seq[0]:
            ncbi_ref_seq_id = ncbi_ref_seq[0].replace(f"({self.g_name})", "")

        elif ignore_warnings is False:
            warnings.warn(
                f"{self.g_name} not found in variant name {self.variant_name}"
            )

        return ncbi_ref_seq_id

    def get_molecular_consequence_list(
        self,
        missense_only: bool = True,
    ) -> list:
        """ Get the molecular consequence list from the HGVS list for a given
        RefSeq ID.

        At least one molecular consequence should be "missense variant".

        Args:

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

            # BUG: potential bug
            # can cause issues if the aa sequence changed in the newer
            # versions of same base RefSeq ID
            # "@sequenceAccessionVersion" provides the full ID, but it can be
            # different from the one in the title
            if (
                self.ncbi_ref_seq_id is None
                or seq_accession != self.ncbi_ref_seq_id.split(".")[0]
            ):
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

    def get_mutation_descs(self, ignore_warnings: bool = True) -> list:
        """ Get the mutation descriptors from the variant name.

        Assuming the format "NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)"
        the mutation descriptors are in the part after the colon,
        i.e. RefSeqID(GeneName):c.DNAChange (p.ProteinChange)
        then, DNAChange and ProteinChange are the mutation descriptors.

        Args:

            ignore_warnings (bool, optional):
                Whether to ignore warnings. Defaults to True.

        Returns:

            list:
                List of mutation descriptors.
        """

        mutation_descs = []
        ncbi_ref_seq = self.variant_name.split(":")

        if len(ncbi_ref_seq) > 2:
            print(f"Unexpected variant name format: {self.variant_name}")

        if self.g_name in ncbi_ref_seq[0] and len(ncbi_ref_seq) == 2:
            mutation_descs = ncbi_ref_seq[1].split(" ")

        elif ignore_warnings is False:
            warnings.warn(
                f"{self.g_name} not found in variant name {self.variant_name}"
            )

        return mutation_descs

    def get_ncbi_g_name(self, ignore_warnings: bool = True) -> str:
        """ Get the NCBI gene name from the variant name.

        Assuming the format "NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)"
        the gene name is in the part before the colon,
        i.e. RefSeqID(GeneName):c.DNAChange (p.ProteinChange)
        then, the gene name is GeneName.

        Args:

            ignore_warnings (bool, optional):
                Whether to ignore warnings. Defaults to True.

        Returns:

            str:
                The NCBI gene name.
        """

        ncbi_g_name = None
        ncbi_ref_seq = self.variant_name.split(":")

        if g_name in ncbi_ref_seq[0]:
            ncbi_g_name = ncbi_ref_seq[0].split("(")[1].replace(")", "")

        elif ignore_warnings is False:
            warnings.warn(
                f"{g_name} not found in variant name {self.variant_name}"
            )

        return ncbi_g_name

    def get_p_mutation(self) -> str:
        """ Get the protein mutation from the mutation descriptors.

        Assumes format "NM_001943.5(DSG2):c.137G>T (p.Arg46Leu)"
        i.e. RefSeqID(GeneName):c.DNAChange (p.ProteinChange)
        The protein mutation is Arg46Leu in this case.

        Returns:

            str:
                The protein mutation.
        """

        p_mutation = None

        for desc in self.mutation_descs:
            desc = desc.replace("(", "").replace(")", "")
            if desc.startswith("p."):
                p_mutation = desc.replace("p.", "")

        return p_mutation

    def set_trait_set(self):
        """ Set the trait set from the germline classification."""

        trait_set = self.germline_classification.get("ConditionList", {}).get(
            "TraitSet", []
        )

        setter_type = type(trait_set)
        assert setter_type in VariantInfo.setter_dict, (
            f"Unexpected trait set format: {setter_type}"
            f"for trait_set: {trait_set}"
            f"in variant ID: {self.variant_id}"
        )

        self.trait_set = VariantInfo.setter_dict[setter_type](trait_set)

    @staticmethod
    def is_contributing_disease_trait(trait_dict: dict) -> bool:
        """ Check if a trait is a disease and contributes to agg. classification.

        Args:

            trait_dict (dict):
                Trait dictionary from ClinVar variant information.

        Returns:

            bool:
                True if the trait is a disease and contributes to agg.
        """

        _t_type = trait_dict.get("@Type", "")
        _t_contributes = trait_dict.get(
            "@ContributesToAggregateClassification", "false"
        )

        return (_t_type == "Disease" and _t_contributes == "true")

    @staticmethod
    def is_preferred_trait_name(trait_name: dict) -> bool:
        """ Check if a trait name is preferred.

        Args:

            trait_name (dict):
                Trait name dictionary from ClinVar variant information.

        Returns:

            bool:
                True if the trait name is preferred, False otherwise.
        """

        return (
            trait_name.get("ElementValue", {}).get("@Type", "")
            == "Preferred"
        )

    @staticmethod
    def extract_preferred_trait_names(trait_dict: dict) -> list:
        """ Extract preferred trait names from a trait dictionary.

        Args:

            trait_dict (dict):
                Trait dictionary from ClinVar variant information.

        Returns:

            list:
                List of preferred trait names.
        """

        preferred_trait_names = []
        trait_ = trait_dict.get("Trait", {})

        trait_type = type(trait_)
        assert trait_type in VariantInfo.setter_dict, (
            f"Unexpected trait format: {trait_type}"
            f"for trait: {trait_}"
            f"in variant ID: {trait_dict}"
        )

        trait_ = VariantInfo.setter_dict[trait_type](trait_)

        for _trait_ in trait_:

            trait_names = _trait_.get("Name", [])
            trait_name_type = type(trait_names)
            assert trait_name_type in VariantInfo.setter_dict, (
                f"Unexpected trait name format: {trait_name_type}"
                f"for trait_names: {trait_names}"
                f"in variant ID: {trait_dict}"
            )
            trait_names = VariantInfo.setter_dict[trait_name_type](trait_names)

            for trait_name in trait_names:

                if not isinstance(trait_name, dict):
                    continue

                if VariantInfo.is_preferred_trait_name(trait_name):
                    preferred_trait_names.append(
                        trait_name.get("ElementValue", {}).get("#text", "")
                    )

        return preferred_trait_names

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

        for trait_dict in self.trait_set:
            # trait should be a disease and contribute to agg. classification
            if VariantInfo.is_contributing_disease_trait(trait_dict) is False:
                continue

            preferred_trait_names = VariantInfo.extract_preferred_trait_names(
                trait_dict=trait_dict
            )

            traits.extend(preferred_trait_names)

        if len(set(traits)) == 0:
            traits = ["not provided"]

        return list(set(traits))

    def set_clinical_assertion_list(self):
        """ Set the clinical assertion list from the classified record."""

        assertions = self.classified_record.get(
            "ClinicalAssertionList", {}
        ).get("ClinicalAssertion", [])


        assertions_type = type(assertions)

        assert assertions_type in VariantInfo.setter_dict, (
            f"Unexpected clinical assertion format: {assertions_type}"
            f"for clinical_assertion_list: {assertions}"
            f"in variant ID: {self.variant_id}"
        )

        self.clinical_assertion_list = (
            VariantInfo.setter_dict[assertions_type](assertions)
        )

    def get_all_clinical_assertions(
        self,
        sort_by_date: bool=True,
        date_format:str=DATE_FORMAT,
    ) -> list:
        """ Get all clinical assertions from a classified record.

        Clinical assertion represents a single submission to ClinVar corresponding
        to the variant. Multiple submissions can be present for a single variant.
        We only look at those that contribute to the agg. clinical significance.

        if `sort_by_date` is True, the assertions are sorted by the date they
        were last evaluated. If the evaluation date is not available, the date
        last updated is used instead.

        Look through all the assertions that contributed to aggregate clinical
        significance

        Args:

            sort_by_date (bool, optional):
                Whether to sort the assertions by date. Defaults to True.

            date_format (str, optional):
                Date format for parsing dates. Defaults to DATE_FORMAT.

        Returns:

            list:
                List of clinical assertions.
        """

        clinical_assertions = []

        for clinical_assertion in self.clinical_assertion_list:
            # we look at only those that contribute
            contributes_to_agg_class = clinical_assertion.get(
                "@ContributesToAggregateClassification", "false"
            )
            if contributes_to_agg_class == "false":
                continue

            last_updated = clinical_assertion.get("@DateLastUpdated", "")
            assertion_class = clinical_assertion.get("Classification", {})

            germline_class = assertion_class.get("GermlineClassification", "")
            last_evaluated = assertion_class.get("@DateLastEvaluated", "")
            assertion_comment = assertion_class.get("Comment", "")

            if isinstance(assertion_comment, dict):
                assertion_comment = assertion_comment.get("#text", "")

            last_updated = datetime.strptime(last_updated, date_format)
            try:
                last_evaluated = datetime.strptime(last_evaluated, date_format)
            except ValueError:
                last_evaluated = last_updated

            clinical_assertions.append({
                "last_updated_date": last_updated,
                "asserted_germline_class": germline_class,
                "assertion_comment": assertion_comment,
                "last_evaluated_date": last_evaluated,
            })

        if sort_by_date:
            clinical_assertions = sorted(
                clinical_assertions,
                key=lambda x: x["last_evaluated_date"],
                reverse=True
            )

        return clinical_assertions

    def get_all_significances(self) -> list:
        """ Get all clinical significances from the clinical assertions.

        Returns:

            list:
                List of clinical significances from all studies.
        """

        all_significances = [""]

        if len(self.clinical_assertions) > 0:

            all_significances = [
                self.clinical_assertions[i]["asserted_germline_class"]
                for i in range(len(self.clinical_assertions))
            ]

        return all_significances

    def get_all_assertion_comments(self) -> list:
        """ Get all assertion comments from the clinical assertions.

        Returns:

            list:
                List of assertion comments provided for all studies.
        """

        all_comments = [""]

        if len(self.clinical_assertions) > 0:

            all_comments = [
                self.clinical_assertions[i]["assertion_comment"]
                for i in range(len(self.clinical_assertions))
            ]
            assert len(all_comments) == len(self.all_significances), (
                "Length of assertion comments and significances do not match."
            )
            all_comments = [
                f"{sign} :- {comment}"
                for sign, comment in zip(self.all_significances, all_comments)
            ]

        return all_comments

    def is_invalid_variant(self) -> bool:
        """ Check if the variant is invalid based on certain conditions.

        The conditions checked are:
        - Aggregate clinical significance is not in allowed list.
        - NCBI RefSeq ID is None.
        - Molecular consequence is not valid (i.e. not a missense variant).
        - Provided gene name doesn't match NCBI gene name (weird but true, DSC2)
        - Mutation is not a single missense mutation.
        - If agg. significance is "Conflicting classifications of pathogenicity",
            the most recent significance is not in allowed list.

        Returns:

            bool:
                True if the variant is invalid, False otherwise.
        """

        condition_val = (
            self.agg_significance not in CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE
            or self.ncbi_ref_seq_id is None
            or len(self.molecular_consequence_list) == 0
            or self.g_name != self.ncbi_g_name
            or is_missense_mutation(self.p_mutation) is False
            or (self.agg_significance
                == "Conflicting classifications of pathogenicity"
                and self.all_significances[0]
                    not in CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE)
        )

        return condition_val

    def update_p_mutation(
        self,
        modeled_seq: str,
        ref_seq_file: str,
        pairwise_alignment_file: str,
        ignore_warnings: bool = True,
    ):
        """ Update the protein mutation based on the modeled sequence.

        The ClinVar protein mutation is based on the NCBI reference sequence.
        The modeled sequence might differ due to different isoforms, this
        function updates the residue number in the protein mutation accordingly.

        Args:

            modeled_seq (str):
                Sequence of the modeled protein.

            ref_seq_file (str):
                Path to the NCBI reference sequence JSON file.

            pairwise_alignment_file (str):
                Path to the pairwise alignment file.

            ignore_warnings (bool, optional):
                Whether to ignore warnings. Defaults to True.

        Returns:

            str | None:
                Updated protein mutation or None if p_mutation is None.
        """

        if self.p_mutation is None:
            return None

        p_sequence = get_ncbi_ref_seq(vi.ncbi_ref_seq_id, ref_seq_file)

        clinvar_psa_map = handle_pairwise_alignment(
            p_name=self.p_name,
            sseq=p_sequence,
            qseq=modeled_seq,
            pairwise_alignment_file=pairwise_alignment_file,
            ignore_warnings=True
        )

        split_mut = split_missense_mutation(
            p_mutation=self.p_mutation,
            return_type="all"
        )

        if split_mut is None:
            self.p_mutation = None
            return None

        wt_aa, res_num, mut_aa = split_mut
        res_num = int(res_num)

        res_num_mapped, warn_msg = get_mapped_residue(
            psa_map=clinvar_psa_map,
            codon_number=res_num,
        )

        if ignore_warnings is False and len(warn_msg) > 0:
            warnings.warn(f"""{warn_msg} (Protein: {p_name})""")

        if len(clinvar_psa_map) == 0:
            res_num_mapped = res_num
        else:
            try:
                res_num_mapped = clinvar_psa_map[res_num]
            except KeyError:
                warnings.warn(
                    f"""
                    Residue number {res_num} not found in
                    pairwise alignment map for protein {p_name}.
                    Skipping...
                    """
                ) if ignore_warnings is False else None
                res_num_mapped = None

        self.p_mutation = (
            None if res_num_mapped is None
            else f"{wt_aa}{res_num_mapped}{mut_aa}"
        )

        return self.p_mutation

    def make_variant_dict(self):
        """ Create a dictionary representation of the variant information."""

        self.variant_dict = {
            "variant_id": self.variant_id,
            "protein": self.p_name,
            "gene": self.g_name,
            "ncbi_ref_seq_id": self.ncbi_ref_seq_id,
            "p_mutation": self.p_mutation,
            "molecular_consequence": "\n".join(self.molecular_consequence_list),
            "trait_or_effect": "\n".join(self.traits),
            "clinical_significance": self.agg_significance,
            "all_significances": "\n\n".join(self.all_significances),
            "all_assertion_comments": "\n\n".join(self.all_assertion_comments),
            "variant_type": self.variant_type,
            # "afm_patho_score": afm_patho_score,
            # "afm_pathogenicity": afm_pathogenicity,
        }

    def add_to_variant_dict(
        self,
        key: str,
        value: Any,
        overwrite: bool = True
    ):
        """ Add a key-value pair to the variant dictionary.

        Args:

            key (str):
                The key to add to the variant dictionary.

            value (Any):
                The value to add to the variant dictionary.

            overwrite (bool, optional):
                Whether to overwrite the existing key in the variant dictionary.
                Defaults to True.
        """

        warn_msg = ""

        if not hasattr(self, "variant_dict"):
            warn_msg += "variant_dict not found, creating new one. "
            self.make_variant_dict()

        warn_msg_dict = {
            True: "Overwriting existing key in variant_dict. ",
            False: "Not overwriting existing key in variant_dict. ",
        }

        if key in self.variant_dict:
            warn_msg += f"Key {key} already exists in variant_dict. "
            warn_msg += warn_msg_dict[overwrite]
            if overwrite is False:
                warnings.warn(warn_msg)
                return

        if warn_msg != "":
            warnings.warn(warn_msg)

        self.variant_dict[key] = value

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
        "--af_missense_tsv",
        type=str,
        required=False,
        default=AF_MISSENSE_AA_SUBSTITUTIONS_TSV,
        help="Path to AlphaMissense aa substitutions TSV file for offline mode.",
    )
    parser.add_argument(
        "--odp_sequences_fasta",
        type=str,
        default="/home/omkar/Projects/cardiac_desmosome/data/sequences/odp_protein_sequences.fasta",
        # default=odp_sequences_fasta,
        help="Fasta file containing modeled protein sequences.",
    )
    args = parser.parse_args()

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

    #######################################################################
    # Fetch AlphaMissense data
    #######################################################################
    if args.include_AF_missense:

        af_missense_df_gen = fetch_af_missense_data(
            args.alpha_missense_dir,
            uniprot_bases,
            mode=args.af_missense_mode,
            overwrite=False,
            af_missense_tsv=args.af_missense_tsv,
        )

        for af_missense_df, uniprot_base in af_missense_df_gen:

            if af_missense_df.empty:
                warnings.warn(f"No AlphaMissense data found for {uniprot_base}.")
                continue

            AF_MISSENSE_CSV = os.path.join(
                args.alpha_missense_dir, f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
            )

            if os.path.exists(AF_MISSENSE_CSV) and not args.overwrite:
                warnings.warn(f"File {AF_MISSENSE_CSV} already exists.")
                continue

            af_missense_df.to_csv(AF_MISSENSE_CSV, index=False)

            print(f"Saved AlphaMissense variants to {AF_MISSENSE_CSV}")

    for p_name, uniprot_id in protein_uniprot_map.items():

        print(f"Processing {p_name}...")

        uniprot_base = uniprot_id.split("-")[0]
        g_name = protein_gene_map[p_name]

        AF_MISSENSE_CSV = os.path.join(
            args.alpha_missense_dir, f"{uniprot_base}{AF_MISSENSE_CSV_SUFFIX}.csv"
        )
        AF_MISSENSE_PAIR_ALN_FASTA = os.path.join(
            args.pairwise_alignments_dir, f"{p_name}{AF_MISSENSE_PAIR_ALN_SUFFIX}.fasta"
        )
        CLINVAR_VARIANTS_JSON = os.path.join(
            args.clinvar_output_dir, f"{g_name}_clinvar_variants1.json"
        )
        CLINVAR_VARIANT_IDS_JSON = os.path.join(
            args.clinvar_output_dir, f"{g_name}_clinvar_variant_ids.json"
        )
        CLINVAR_PAIR_ALN_FASTA = os.path.join(
            args.pairwise_alignments_dir, f"{p_name}_clinvar_vs_modeled.fasta"
        )

        af_missense_dict = {}
        modeled_seq = odp_sequences.get(protein_uniprot_map[p_name], None)
        if modeled_seq is None:
            warnings.warn(f"Sequence not found for {p_name} ({uniprot_id}).")
            continue

        #######################################################################
        # Get AlphaMissense scores for variants in the protein
        #######################################################################
        if args.include_AF_missense:

            af_missense_ref_seq = fasta_dict.get(uniprot_base, None)
            if af_missense_ref_seq is None:
                warnings.warn(f"Sequence not found for {p_name} ({uniprot_base}).")
                continue

            # won't align if sequences are identical
            afm_psa_map = handle_pairwise_alignment(
                p_name=p_name,
                sseq=af_missense_ref_seq,
                qseq=modeled_seq,
                pairwise_alignment_file=AF_MISSENSE_PAIR_ALN_FASTA,
                ignore_warnings=True
            )

            if not os.path.exists(AF_MISSENSE_CSV):
                warnings.warn(
                    f"""AlphaMissense data file not found for {p_name}:
                    {AF_MISSENSE_CSV}"""
                )
                continue

            af_missense_df = pd.read_csv(AF_MISSENSE_CSV)

            af_missense_dict = af_missense_df_to_dict(
                p_name,
                af_missense_df,
                af_missense_dict=af_missense_dict,
                afm_psa_map=afm_psa_map, # may not be appropriate to use in some cases
                ignore_warnings=True,
            )

        print(f"Fetching ClinVar variants for {g_name}...")
        #######################################################################
        # Fetch variant ids from ClinVar
        #######################################################################
        api_parameters = CLINVAR_TEMPLATE_QUERY_ID.copy()
        api_parameters["term"] = Template(api_parameters["term"]).substitute(
            gene_name=g_name
        )

        variant_ids = fetch_clinvar_variant_data(
            save_path=CLINVAR_VARIANT_IDS_JSON,
            input_data=g_name,
            api_parameters=api_parameters,
            ignore_error=False,
            overwrite=False,
        )

        if len(variant_ids) == 0:
            print(f"No variants found for {g_name}")
            continue

        #######################################################################
        # Fetch variant details from ClinVar
        #######################################################################
        clinvar_variants = fetch_clinvar_variant_data(
            save_path=CLINVAR_VARIANTS_JSON,
            input_data=variant_ids,
            api_parameters=CLINVAR_TEMPLATE_QUERY_DETAIL,
            ignore_error=False,
            overwrite=False,
        )


        #######################################################################
        # Collect relevant information about missense variants
        #######################################################################
        ncbi_ref_seq_ids = set()
        for variant_id, variant_info in clinvar_variants.items():

            vi = VariantInfo(
                p_name=p_name,
                g_name=g_name,
                variant_id=variant_id,
                variant_info=variant_info,
            )

            REF_SEQ_JSON = os.path.join(
                args.clinvar_output_dir, f"{g_name}_{vi.ncbi_ref_seq_id}_nuccore.json"
            )

            vi.update_p_mutation(
                modeled_seq=modeled_seq,
                ref_seq_file=REF_SEQ_JSON,
                pairwise_alignment_file=CLINVAR_PAIR_ALN_FASTA,
                ignore_warnings=True,
            )

            if vi.is_invalid_variant():
                continue

            vi.make_variant_dict()
            vi.add_to_variant_dict("uniprot_id", uniprot_id)

            df_rows.append(vi.variant_dict)
            ncbi_ref_seq_ids.add(vi.ncbi_ref_seq_id)

        if len(ncbi_ref_seq_ids) == 0:
            print(f"No missense variants found for {g_name}")
            continue

        if len(ncbi_ref_seq_ids) != 1:
            raise NotImplementedError(
                f"Multiple RefSeq IDs found for {g_name}: {ncbi_ref_seq_ids}."
            )

    df = pd.DataFrame(df_rows)
    df = df.rename(columns=CLINVAR_DF_COLUMNS)
    # print(df.head())

    out_name = "clinvar_missense_variants"
    if args.include_VUS:
        out_name += "_with_VUS"

    df_file = os.path.join(
        args.clinvar_output_dir, f"{out_name}.xlsx"
    )
    df.to_excel(df_file, index=False)

    print(f"ClinVar missense variants saved to {df_file}")