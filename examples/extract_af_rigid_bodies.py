import sys
from argparse import ArgumentParser
import pandas as pd
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.AFoutput import RigidBodies
from utils import read_json, write_json
import os

if __name__ == "__main__":
    args = ArgumentParser()
    args.add_argument(
        "-s",
        "--structure",
        type=str,
        required=False,
        default="./input/af_structure.cif",
        help="Path to AF2/3 structure file",
    )
    args.add_argument(
        "-d",
        "--data",
        type=str,
        required=False,
        default="./input/af_data.json",
        help="Path to AF2/3 data file",
    )
    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="./output/af_output",
        help="Path to output directory",
    )
    args.add_argument(
        "--modeled_region",
        type=str,
        required=False,
        help="Modeled region in the form of a dictionary",
    )
    args = args.parse_args()

    af = RigidBodies(
        structure_path=args.structure,
        data_path=args.data,
        output_dir=args.output,
        # modeled_region=read_json(args.modeled_region),
    )

    af.plddt_cutoff = 70
    af.pae_cutoff = 5
    af.pae_power = 1
    af.resolution = 0.1
    af.library = "igraph"

    os.makedirs(args.output, exist_ok=True)

    domains = af.predict_domains()
    # domains = af.plddt_filtered_domains(domains)
    # domains = pd.DataFrame(domains)
    # domains.to_csv(f"{args.output}/rigid_bodies.csv", index=False)
    rbs = af.get_rigid_bodies(domains, num_proteins=2)
    write_json(f"{args.output}/rigid_bodies.json", rbs)