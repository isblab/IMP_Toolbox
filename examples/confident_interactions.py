import sys
from argparse import ArgumentParser
from matplotlib import pyplot as plt
import pandas as pd
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.AFoutput import ContactMap
from utils import write_json, get_key_from_res_range
import os
import yaml

if __name__ == "__main__":
    args = ArgumentParser()
    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="./input/af_output_analysis.yaml",
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
    # print(input_yaml)
    
    for pred_to_analyse in input_yaml:
        pred_selection = pred_to_analyse.get("selection")
        structure_path = pred_to_analyse.get("structure_path")
        data_path = pred_to_analyse.get("data_path")
        file_name = os.path.basename(structure_path).split(".")[0]
        
        af_contacts =  ContactMap(
            structure_path=structure_path,
            data_path=data_path,
            pred_selection=pred_selection,
        )
        af_contacts.contact_threshold = 8
        af_contacts.interaction_map_type = "contact"
        interacting_regions = af_contacts.get_interacting_regions()
        for interacting_region in interacting_regions:
            contact_map = af_contacts.get_contact_map(interacting_region)
            # plt.imshow(contact_map, cmap="viridis", interpolation="nearest")
            # plt.show()
            # interacting_residues = af_contacts.get_contacting_residues(interacting_region)
            segmented_map = af_contacts.get_interacting_patches(interacting_region, bandwidth=15)
    # print(interacting_regions)
    # plt.imshow(contact_map, cmap="viridis", interpolation="nearest")
    # print(interacting_residues)
    plt.imshow(segmented_map, cmap="tab20", interpolation="nearest")
    plt.show()