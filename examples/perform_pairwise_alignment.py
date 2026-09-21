import argparse
from IMP_Toolbox.sequence.sequence_alignment import PairwiseSequenceAlignment
from IMP_Toolbox.sequence.sequence import (
    query_uniprot_api_for_sequences,
    only_uniprot_id_as_header,
    fasta_str_to_dict,
)
# for P60709 and Q92747

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-p1",
        "--protein1",
        type=str,
        required=True,
        help="Uniprot id of the first protein",
    )

    parser.add_argument(
        "-p2",
        "--protein2",
        type=str,
        required=True,
        help="Uniprot id of the second protein",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/pairwise_alignment.fasta",
        help="Path to output fasta file containing pairwise alignment of the two proteins",
    )

    args = parser.parse_args()

    fasta = query_uniprot_api_for_sequences(uniprot_ids=[args.protein1, args.protein2])
    fasta = only_uniprot_id_as_header(fasta_str=fasta)
    sequences = fasta_str_to_dict(fasta)
    seq1, seq2 = sequences[args.protein1], sequences[args.protein2]

    pairwise_alignment = PairwiseSequenceAlignment(
        seq1=seq1,
        seq2=seq2,
        moltype="prot",
        program="stretcher",
        header1=args.protein1, # optional, default is "query 1-${end}"
        header2=args.protein2, # optional, default is "subject 1-${end}"
    )

    psa_map = pairwise_alignment.fetch_pairwise_alingment_map(
        pairwise_alignment_file=args.output,
        overwrite=True
    )
    print(pairwise_alignment.pairwise_alignment)

    seq1_aln = pairwise_alignment.get_alignment_attribute(attribute="qaln")
    seq2_aln = pairwise_alignment.get_alignment_attribute(attribute="saln")
    print(f"Pairwise alignment of {args.protein1} and {args.protein2}:")
    print(seq1_aln)
    print(seq2_aln)
    print("*"*50)

    percent_identity = pairwise_alignment.get_alignment_attribute(attribute="pidentity")
    print(f"Percentage identity: {percent_identity}%")

    percent_similarity = pairwise_alignment.get_alignment_attribute(attribute="psimilarity")
    print(f"Percentage similarity: {percent_similarity}%")

    # See all alignment attributes in IMP_Toolbox.constants.sequence_constants.PSAAttribute
    gap = PairwiseSequenceAlignment.get_gap(
        qseq=seq1_aln,
        sseq=seq2_aln,
        start=1,
        end=len(seq1),
        reference="qseq",
        as_percentage=True,
    )
    print(f"Percentage gap: {gap}%")