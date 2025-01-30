from pprint import pprint
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

    input_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)

    interacting_regions = [
        [
            {
                "A": (50, 200),
                "D": (1, 9),
            },
            {
                "A": (600, 700),
                "E": (1, 20),
            },
            {
                "A": (1, 1200),
                "C": (1, 25)
            }
        ],
        [
            {
                "A": (1, 375),
                "B": (1, 127),
            }
        ]
    ]

    for _, pred_to_analyse in enumerate(input_yaml):
        pred_selection = pred_to_analyse.get("af_offset")
        structure_path = pred_to_analyse.get("structure_path")
        data_path = pred_to_analyse.get("data_path")

        dir_name = os.path.basename(structure_path).split(".")[0]
        output_dir = os.path.join(args.output, f"{dir_name}_patches")
        os.makedirs(output_dir, exist_ok=True)

        af_interaction = Interaction(
            struct_file_path=structure_path,
            data_file_path=data_path,
        )

        af_interaction.plddt_cutoff = 70
        af_interaction.pae_cutoff = 5
        af_interaction.interaction_map_type = "contact"
        af_interaction.contact_threshold = 8

        for interacting_region in interacting_regions[_]:

            contact_map = af_interaction.get_confident_interactions(
                interacting_region=interacting_region
            )

            seg_map, interacting_patches = af_interaction.get_interacting_patches(
                contact_map=contact_map,
                interacting_region=interacting_region,
                bandwidth=None,
            )

            ir_str = "_".join([f"{k}:{v[0]}-{v[1]}" for k, v in interacting_region.items()])

            pprint(interacting_patches)

            af_interaction.save_map(
                contact_map=seg_map,
                patches=interacting_patches,
                interacting_region=interacting_region,
                out_file=os.path.join(output_dir, f"patch_{ir_str}.html"),
                save_plot=True,
            )