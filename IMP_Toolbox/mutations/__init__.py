"""
mutations
===========

- Assisting with the mutation related tasks such as-
  - fetching and processing mutation data from alpha-missense, ClinVar, etc.

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