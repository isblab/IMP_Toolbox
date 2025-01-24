import sys
from argparse import ArgumentParser
import yaml
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.AFInput import AFInput
from utils import read_json, read_fasta

if __name__ == "__main__":
    args = ArgumentParser()
    args.add_argument(
        "-p",
        "--prtotein_sequences",
        type=str,
        required=False,
        default="./output/protein_sequences.fasta",
        help="fasta file containing all sequences"
    )
    args.add_argument(
        "-u",
        "--uniprot",
        type=str,
        required=False,
        default="./input/proteins.json",
        help="json file containing protein names and uniprot ids"
    )
    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/af_input",
        help="output directory for alphafold input"
    )
    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/af_server_targets.yaml",
        help="input yaml file containing the target proteins"
    )
    args.add_argument(
        "-n",
        "--nucleotide_sequences",
        type=str,
        required=False,
        default=None,
        help="fasta file containing dna or rna sequences"
    )
    args = args.parse_args()
    proteins = read_json(args.uniprot)
    protein_sequences = read_fasta(args.prtotein_sequences)
    dna_sequences = read_fasta(args.nucleotide_sequences) if args.nucleotide_sequences else None
    input_yml = yaml.load(open("./input/af_server_targets.yaml"), Loader=yaml.FullLoader)

    af_input = AFInput(proteins, protein_sequences, input_yml, dna_sequences)
    af_input.output_dir = args.output
    af_input.sanity_check_uniprot(proteins)
    af_input.create_job_files()