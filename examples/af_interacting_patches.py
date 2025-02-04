import sys
from argparse import ArgumentParser
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.Interaction import Interaction
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
        default="./output/af_output/interacting_patches",
        help="Path to output directory",
    )

    args = args.parse_args()

    input_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)

    regions_of_interest = [
        [
            # {
            #     "A": (50, 200),
            #     "D": (1, 9),
            # },
            # {
            #     "A": (600, 700),
            #     "E": (1, 20),
            # },
            # {
            #     "A": (1, 1200),
            #     "C": (1, 25)
            # }
        ],
        [
            # {
            #     "A": (1, 375),
            #     "B": (1, 127),
            # }
        ],
        [
            {
                "A": (1600, 1800),
                "B": (1600, 1800),
            }
        ],
    ]

    for pred_idx, pred_to_analyse in enumerate(input_yaml):
        af_offset = pred_to_analyse.get("af_offset")
        structure_path = pred_to_analyse.get("structure_path")
        data_path = pred_to_analyse.get("data_path")

        regions_of_interest_ = regions_of_interest[pred_idx]

        af_interaction = Interaction(
            struct_file_path=structure_path,
            data_file_path=data_path,
            af_offset=af_offset,
            output_dir=args.output,
        )

        af_interaction.plddt_cutoff = 70
        af_interaction.pae_cutoff = 10
        af_interaction.interaction_map_type = "contact"
        af_interaction.contact_threshold = 10

        if not regions_of_interest_:
            regions_of_interest_ = af_interaction.create_regions_of_interest()

        for region_of_interest in regions_of_interest_:
            af_interaction.save_ppair_interaction(
                region_of_interest=region_of_interest,
                save_plot=True,
                plot_type="both",
            )