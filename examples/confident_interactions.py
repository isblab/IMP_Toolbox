import sys
from argparse import ArgumentParser
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.Interaction import Interaction
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
        pred_selection = pred_to_analyse.get("af_offset")
        structure_path = pred_to_analyse.get("structure_path")

        data_path = pred_to_analyse.get("data_path")
        file_name = os.path.basename(structure_path).split(".")[0]

        af_interaction = Interaction(
            struct_file_path=structure_path,
            data_file_path=data_path,
        )

        af_interaction.plddt_cutoff = 70
        af_interaction.pae_cutoff = 5
        af_interaction.interaction_map_type = "contact"
        af_interaction.contact_threshold = 8
        interacting_region = {
            "A": (50, 200),
            "D": (1, 9),
        }

        contact_map = af_interaction.get_confident_interactions(
            interacting_region=interacting_region
        )

        seg_map, interacting_patches = af_interaction.get_interacting_patches(
            contact_map=contact_map,
            interacting_region=interacting_region,
            bandwidth=2,
        )

        print(seg_map.shape, interacting_patches)