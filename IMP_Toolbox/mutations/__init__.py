"""
mutations
===========

- Assisting with the mutation related tasks such as-
  - fetching and processing mutation data from alpha-missense, ClinVar, etc.

### Use cases:

#### 1. Obtain AlphaMissense data for given UniProt IDs

```bash
python af_missense.py \\
  --uniprot_ids P60709,P68133 \\
  --alpha_missense_dir /path/to/alpha_missense_dir \\
  --mode online
```

- This command fetches AlphaMissense data for the specified UniProt IDs
  (`P60709` and `P68133`) and saves the data in the specified directory
  (`/path/to/alpha_missense_dir`).

- The `--mode online` flag indicates that the data should be fetched from the
  online AlphaMissense API. Alternatively, you can use `--mode offline` to fetch
  data from a local tsv file. To obtain this file you only need to run the command
  with `--mode offline --af_missense_tsv /path/to/AlphaMissense_aa_substitutions.tsv.gz`
  once, and it will download the file to the specified path.

### Classes

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
    class VariantInfo {
        + dict setter_dict
        - \_\_init__(self, p_name, g_name, variant_id, variant_info) None
        + set_variant_name(self)
        + set_variant_type(self)
        + set_classified_record(self) dict
        + set_hgvs_list(self) dict
        + set_germline_classification(self) dict
        + set_agg_significance(self) str
        + get_ncbi_ref_seq_id(self, ignore_warnings) str | None
        + get_molecular_consequence_list(self, missense_only) list
        + get_mutation_descs(self, ignore_warnings) list
        + get_ncbi_g_name(self, ignore_warnings) str
        + get_p_mutation(self) str
        + set_trait_set(self)
        + @staticmethod is_contributing_disease_trait(trait_dict) bool$
        + @staticmethod is_preferred_trait_name(trait_name) bool$
        + @staticmethod extract_preferred_trait_names(trait_dict) list$
        + get_variant_associated_traits(self) list
        + set_clinical_assertion_list(self)
        + get_all_clinical_assertions(self, sort_by_date, date_format) list
        + get_all_significances(self) list
        + get_all_assertion_comments(self) list
        + is_invalid_variant(self) bool
        + update_p_mutation(self, modeled_seq, ref_seq_file, pairwise_alignment_file, ignore_warnings)
        + make_variant_dict(self)
        + add_to_variant_dict(self, key, value, overwrite)
    }
```

"""