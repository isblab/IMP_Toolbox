import os
import argparse

import pandas as pd
from IMP_Toolbox.utils.file_helpers import read_json, read_fasta
from IMP_Toolbox.mutations.clinvar_mutations import process_clinvar_variant_data

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Fetch missense variants from ClinVar."
    )
    parser.add_argument(
        "--protein_uniprot_map",
        type=str,
        default="./input/proteins.json",
        help="Path to JSON file mapping protein names to UniProt IDs.",
    )
    parser.add_argument(
        "--protein_gene_map",
        type=str,
        default="./input/genes.json",
        help="Path to JSON file mapping protein names to gene names.",
    )
    parser.add_argument(
        "--odp_sequences_fasta",
        type=str,
        default="./output/protein_sequences.fasta",
        help="Fasta file containing modeled protein sequences.",
    )
    parser.add_argument(
        "--include_VUS",
        action="store_true",
        default=False,
        help="Include variants of uncertain significance (VUS).",
    )
    parser.add_argument(
        "--include_AF_missense",
        action="store_true",
        default=False,
        help="Include AlphaMissense pathogenicity scores.",
    )
    parser.add_argument(
        "--clinvar_output_dir",
        type=str,
        default="./output/mutations/clinvar",
        help="Directory to save ClinVar variant information.",
    )
    parser.add_argument(
        "--alpha_missense_dir",
        type=str,
        default="./output/mutations/alpha_missense",
        help="Directory to save AlphaMissense variant information.",
    )
    parser.add_argument(
        "--pairwise_alignments_dir",
        type=str,
        default="./output/mutations/pairwise_alignments",
        help="Directory to save pairwise alignments.",
    )
    parser.add_argument(
        "--af_missense_mode",
        type=str,
        choices=["online", "offline"],
        default="online",
        help=f"Mode to fetch AlphaMissense data. Choices are: %(choices)s",
    )
    parser.add_argument(
        "--af_missense_tsv",
        type=str,
        required=False,
        help="Path to AlphaMissense aa substitutions TSV file for offline mode.",
    )
    args = parser.parse_args()

    protein_uniprot_map = read_json(args.protein_uniprot_map)
    protein_uniprot_map = {k: v for k, v in protein_uniprot_map.items() if v is not None}
    protein_gene_map = read_json(args.protein_gene_map)
    protein_gene_map = {k: v for k, v in protein_gene_map.items() if v is not None}
    protein_sequences = read_fasta(args.odp_sequences_fasta)

    os.makedirs(args.clinvar_output_dir, exist_ok=True)

    df = process_clinvar_variant_data(
        protein_uniprot_map=protein_uniprot_map,
        protein_gene_map=protein_gene_map,
        protein_sequences=protein_sequences,
        clinvar_output_dir=args.clinvar_output_dir,
        pairwise_alignments_dir=args.pairwise_alignments_dir,
        include_VUS=args.include_VUS,
        include_AF_missense=args.include_AF_missense,
        alpha_missense_dir=args.alpha_missense_dir,
        af_missense_mode=args.af_missense_mode,
        af_missense_tsv=args.af_missense_tsv,
    )

    if not isinstance(df, pd.DataFrame) or df.empty:
        print("No missense variants found for any protein. Exiting.")
        exit(0)

    out_name = "clinvar_missense_variants"
    if args.include_VUS:
        out_name += "_with_VUS"

    df_file = os.path.join(args.clinvar_output_dir, f"{out_name}.xlsx")
    df.to_excel(df_file, index=False)

    print(f"ClinVar missense variants saved to {df_file}")