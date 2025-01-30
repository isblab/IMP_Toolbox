import sys
from argparse import ArgumentParser
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.RigidBodies import RigidBodies
import os
import yaml

if __name__ == "__main__":

    args = ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/af_predictions.yaml",
        help="Path to input yaml file",
    )

    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/af_output",
        help="Path to output directory",
    )

    args = args.parse_args()

    os.makedirs(args.output, exist_ok=True)

    input_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)

    for pred_to_analyse in input_yaml:
        structure_path = pred_to_analyse.get("structure_path")
        data_path = pred_to_analyse.get("data_path")
        af_offset = pred_to_analyse.get("af_offset")

        file_name = os.path.basename(structure_path).split(".")[0]

        af_rigid = RigidBodies(
            structure_path=structure_path,
            data_path=data_path,
            af_offset=af_offset,
        )

        af_rigid.plddt_cutoff = 70
        af_rigid.pae_cutoff = 5
        af_rigid.pae_power = 1
        af_rigid.resolution = 0.5
        af_rigid.library = "igraph" # "networkx" is slower

        domains = af_rigid.predict_domains(
            num_res=5, # minimum number of residues in a domain
            num_proteins=1, # minimum number of proteins in a domain
            plddt_filter=True, # filter domains based on pLDDT score
        )

        # save the rigid bodies in PDB format
        af_rigid.save_rb_structures(
            domains=domains,
            output_dir=args.output
        )

        # save the rigid bodies in txt format
        af_rigid.save_rb(
            domains=domains,
            output_dir=args.output,
            output_format="txt"
        )
        print(f"saved rigid bodies for {file_name} to: {args.output}")