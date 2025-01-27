import sys
from argparse import ArgumentParser
import pandas as pd
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.AFoutput import RigidBodies
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
        structure_path = pred_to_analyse.get("structure_path")
        data_path = pred_to_analyse.get("data_path")
        selection = pred_to_analyse.get("selection")

        file_name = os.path.basename(structure_path).split(".")[0]

        af = RigidBodies(
            structure_path=structure_path,
            data_path=data_path,
            output_dir=args.output,
            pred_selection=selection,
        )

        af.plddt_cutoff = 70
        af.pae_cutoff = 5
        af.pae_power = 1
        af.resolution = 0.1
        af.library = "igraph"

        domains = af.predict_domains()
        domains = af.plddt_filtered_domains(domains, selected=True)

        rbs = af.get_rigid_bodies(domains, num_proteins=2, selected=True)
        write_json(f"{args.output}/{file_name}.json", rbs)

        # txt output
        with open(f"{args.output}/{file_name}.txt", "w") as f:
            for idx, rb in enumerate(rbs):
                f.write(f"Rigid body {idx+1}\n")
                for ch_id, res_range in rb.items():
                    res_range = get_key_from_res_range(res_range)
                    f.write(f"{ch_id}: {res_range}\n")
                f.write("\n")

        domains = pd.DataFrame(domains)
        domains.to_csv(f"{args.output}/{file_name}.csv", index=False)