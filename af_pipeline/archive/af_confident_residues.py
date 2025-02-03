import sys
from argparse import ArgumentParser
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.archive.ConfidentPredictions import ConfidentPredictions
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
        default="./output/af_output/confident_residues",
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

        out_file = os.path.join(args.output, file_name + "_confident_regions.pdb")

        af_conf_struct = ConfidentPredictions(
            struct_file_path=structure_path,
            data_file_path=data_path,
            out_file=out_file,
        )

        af_conf_struct.plddt_cutoff = 70
        af_conf_struct.pae_cutoff = 5
        af_conf_struct.apply_plddt = True
        af_conf_struct.apply_pae = True

        af_conf_struct.save_confident_regions()
        print("saved confident residues to:", out_file)