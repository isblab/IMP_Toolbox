import re


DATE_FORMAT = r"%Y-%m-%d"

AA_MAP_THREE_TO_ONE = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V"
}

AA_MAP_ONE_TO_THREE = {v: k for k, v in AA_MAP_THREE_TO_ONE.items()}

AMINO_ACID_MAP = {
    **AA_MAP_THREE_TO_ONE,
    **AA_MAP_ONE_TO_THREE
}

AFM_PATHOGENICITY = {
    "LBen": "Likely benign",
    "Amb": "Uncertain",
    "LPath": "Likely pathogenic",
    "benign": "Likely benign",
    "ambiguous": "Uncertain",
    "pathogenic": "Likely pathogenic",
}
# https://www.ebi.ac.uk/training/online/courses/alphafold/classifying-the-effects-of-missense-variants-using-alphamissense/understanding-pathogenicity-scores-from-alphamissense/
#
# 0 to 0.34: likely benign
# 0.34 to 0.564: uncertain
# 0.564 to 1: likely pathogenic

CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE = [
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
    "Conflicting classifications of pathogenicity",
]

POLYPHEN_ORDER = [
    "probably damaging",
    "possibly damaging",
    "benign",
    "unknown",
]

SIFT_ORDER = [
    "deleterious",
    "deleterious - low confidence",
    "tolerated",
    "tolerated - low confidence",
]

VARIANT_PREDICTORS = [
    "SIFT",
    "PolyPhen",
]

UNIPROT_ALLOWED_CLINICAL_SIGNIFICANCE = CLINVAR_ALLOWED_CLINICAL_SIGNIFICANCE

# add new ids as and when required
UNIPROT_PUBMED_TO_IGNORE = [
    "14607462",
    "16839424",
    "20301308",
    "20301310",
    "20301486",
    "20301717",
    "20301725",
    "21810866",
    "23788249",
    "23994779",
    "25173338",
    "25356965",
    "27854360",
    "34012068",
    "35802134",
    "20301373",
    "20301466",
    "20301690",
    "21267010",
    "35378257"
]

API_URLS = {
    "uniprot_variant": "https://www.ebi.ac.uk/proteins/api/variation/uniprot_id",
    "af_missense_csv": "https://alphafold.ebi.ac.uk/files/AF-uniprot_id-F1-aa-substitutions.csv",
    # "af_missense_json": "https://alphafold.ebi.ac.uk/api/annotations/uniprot_id.json?type=MUTAGEN",
    "af_missense_res": "https://alphamissense.hegelab.org/hotspotapi?uid=uniprot_id&resi=res_num",
    "ncbi_esearch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
    "ncbi_efetch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
    "lovd_variant_json": "https://databases.lovd.nl/shared/api/rest.php/variants/gene_name?format=application/json",
    "lovd_variant_txt": "https://databases.lovd.nl/shared/download/all/gene/gene_name",
}

LOVD_PREFERRED_TRANSCRIPT_IDS = {
    "PKP2": 25790,
    "DSP": 6721,
    "DSG2": 6717,
    "DSC2": 6702,
    "JUP": 10288,
    "FZD4": '00008236',
}

LOVD_ALLOWED_CLINICAL_SIGNIFICANCE = [
    "pathogenic",
    "likeply pathogenic",
    "pathogenic (dominant)",
    "likely pathogenic (dominant)",
]

TRUNCATION_NOTATIONS = ["*", "Ter", "X", "Term"]

ALLOWED_AA = set(AMINO_ACID_MAP.keys()).union(set(TRUNCATION_NOTATIONS))

MISSENSE_REGEX = {
    "three_letter": {
        "regex": r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Term|Ter|X)$",
        "condition": lambda x: (
            "p.(" not in x and "p." not in x and "(" not in x
            and sum(char in ALLOWED_AA for char in re.split(r'\d+', x)) == 2
            and len(re.split(r'\d+', x)[0]) == 3
            and len(re.split(r'\d+', x)[-1]) in [1, 3, 4]
        ),
    },
    "one_letter": {
        "regex": r"^([A-Z])(\d+)([A-Z]|\*|Term|Ter|X)$",
        "condition": lambda x: (
            "p.(" not in x and "p." not in x and "(" not in x
            and sum(char in ALLOWED_AA for char in re.split(r'\d+', x)) == 2
            and len(re.split(r'\d+', x)[0]) == 1
            and len(re.split(r'\d+', x)[-1]) in [1, 3, 4]
        ),
    },
    "p_three_letter": {
        "regex": r"^p\.\(([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Term|Ter|X)\)$",
        "condition": lambda x: "p.(" in x,
    },
    "p_one_letter": {
        "regex": r"^p\.\(([A-Z])(\d+)([A-Z]|\*|Term|Ter|X)\)$",
        "condition": lambda x: "p.(" in x,
    },
    "p_three_letter_noparen": {
        "regex": r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*|Term|Ter|X)$",
        "condition": lambda x: "p." in x and "(" not in x,
    },
    "p_one_letter_noparen": {
        "regex": r"^p\.([A-Z])(\d+)([A-Z]|\*|Term|Ter|X)$",
        "condition": lambda x: "p." in x and "(" not in x,
    }
}

ARVC_ALLOWED_CLINICAL_SIGNIFICANCE = [
    "Pathogenic",
    "Likely pathogenic",
]

ARVC_REF_SEQUENCE_IDS = {
    "PKP2": "NM_004572.4",
    "DSP": "NM_004415.4",
    "DSG2": "NM_001943.5",
    "DSC2": "NM_024422.6",
    "JUP": "NM_021991.4",
}

AF_MISSENSE_CSV_SUFFIX = "_alpha_missense_variants"
AF_MISSENSE_PAIR_ALN_SUFFIX = "_afm_vs_modeled"

# from https://console.cloud.google.com/storage/browser/dm_alphamissense
import os
AF_MISSENSE_AA_SUBSTITUTIONS_TSV = f"/data/{os.getlogin()}/Projects/IMP_Toolbox/AlphaMissense_aa_substitutions.tsv.gz"
# AF_MISSENSE_ISOFORMS_AA_SUBSTITUTIONS_TSV = "/data/omkar/Projects/IMP_Toolbox/AlphaMissense_isoforms_aa_substitutions.tsv.gz"

CLINVAR_TEMPLATE_QUERY_ID = {
    "db": "clinvar",
    "term": "$gene_name[gene]",
    "retmax": 10000,
    "retmode": "json",
}

CLINVAR_TEMPLATE_QUERY_DETAIL = {
    "db": "clinvar",
    "rettype": "vcv",
    "is_variationid": "true",
    "from_esearch": "true"
}

CLINVAR_DF_COLUMNS = {
    "gene": "Gene",
    "uniprot_id": "Uniprot ID",
    "ncbi_ref_seq_id": "NCBI RefSeq ID",
    "p_mutation": "Mutation",
    "trait_or_effect": "Disease association",
    "clinical_significance": "ClinVar clinical significance",
    "all_assertion_comments": "All submission comments",
    "all_significances": "All clinical significances",
    "variant_id": "ClinVar Variant ID",
    "molecular_consequence": "Molecular consequence",
}

AF_MISSENSE_COLUMNS = {
    "afm_patho_score": "AlphaMissense score",
    "afm_pathogenicity": "AlphaMissense pathogenicity",
}