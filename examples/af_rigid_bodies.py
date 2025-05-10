import sys
from argparse import ArgumentParser
from set_up import IMP_TOOLBOX
sys.path.append(IMP_TOOLBOX)
from af_pipeline.RigidBodies import RigidBodies
import os
import yaml
import time

if __name__ == "__main__":
    start = time.time()

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
        default="./output/af_output/rigid_bodies",
        help="Path to output directory",
    )

    args = args.parse_args()

    os.makedirs(args.output, exist_ok=True)

    input_yaml = yaml.load(open(args.input), Loader=yaml.FullLoader)

    for pred_to_analyse in input_yaml:
        structure_path = pred_to_analyse.get("structure_path")
        data_path = pred_to_analyse.get("data_path")
        af_offset = pred_to_analyse.get("af_offset")
        average_atom_pae = pred_to_analyse.get("average_atom_pae", False)
        idr_chains = pred_to_analyse.get("idr_chains", [])

        file_name = os.path.basename(structure_path).split(".")[0]

        af_rigid = RigidBodies(
            structure_path=structure_path,
            data_path=data_path,
            af_offset=af_offset,
            idr_chains=idr_chains, # list of chains that are disordered for e.g. ["A", "B"]
            average_atom_pae=average_atom_pae, # average PTM for all chains
        )

        af_rigid.plddt_cutoff = 70
        af_rigid.plddt_cutoff_idr = 50 # you can set different cutoff for IDR
        af_rigid.pae_cutoff = 12
        af_rigid.pae_power = 1
        af_rigid.resolution = 0.5 # higher will result in strict partitioning hence smaller and more rigid domains
        af_rigid.library = "igraph" # "networkx" is slower
        af_rigid.patch_threshold = 0

        domains = af_rigid.predict_domains(
            num_res=5, # minimum number of residues in a domain
            num_proteins=1, # minimum number of proteins in a domain
            plddt_filter=True, # filter domains based on pLDDT score
        )

        # save the rigid bodies in txt format
        af_rigid.save_rigid_bodies(
            domains=domains,
            output_dir=args.output,
            output_format="txt",
            save_structure=True,
            structure_file_type="pdb", # save structure in pdb format
            no_plddt_filter_for_structure=False, # save structure with pLDDT filter
            pae_plot=True, # save PAE plot
        )
        print(f"saved rigid bodies for {file_name} to: {args.output}")

        # metrics for rigid bodies
        all_interface_residues = af_rigid.get_interface_residues(
            domains=domains,
            contact_threshold=8,
            as_matrix=False
        )

        af_rigid.get_average_pLDDT(domains=domains, chain_type="any")
        # af_rigid.get_average_pLDDT(domains=domains, chain_type="idr")
        # af_rigid.get_average_pLDDT(domains=domains, chain_type="r")

        af_rigid.get_ipLDDT(all_interface_residues=all_interface_residues, interface_type="any-any")
        # af_rigid.get_ipLDDT(all_interface_residues=all_interface_residues, interface_type="idr-idr")
        # af_rigid.get_ipLDDT(all_interface_residues=all_interface_residues, interface_type="r-r")
        # af_rigid.get_ipLDDT(all_interface_residues=all_interface_residues, interface_type="idr-any")
        # af_rigid.get_ipLDDT(all_interface_residues=all_interface_residues, interface_type="r-any")
        # af_rigid.get_ipLDDT(all_interface_residues=all_interface_residues, interface_type="idr-r")

        all_interface_residues = af_rigid.get_interface_residues(
            domains=domains,
            contact_threshold=8,
            as_matrix=True
        )

        af_rigid.get_ipae(
            all_interface_residues=all_interface_residues
        )
    time_taken = time.time() - start
    print(f"Time taken: {time_taken:.2f} seconds")