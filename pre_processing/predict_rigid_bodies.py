# TODO: merge small regions within a domain based on the coarse-grain level in IMP
import argparse
import random
import sys
from pae_to_domains.pae_to_domains import (
    parse_pae_file,
    domains_from_pae_matrix_igraph,
    domains_from_pae_matrix_networkx,
)
import os
from utils import read_json, write_json, get_key_from_res_range

# cd to the directory containing the script

#! THINGS TO NOTE:
# you need parser.py from the af_pipeline repository
# pass your own paths for the input_file, predictions_dir and result_dir as arguments
# this script will only work for monomer predictions
# job_name must be "{protein_name}_{everything_else}" and
# protein_name should be the same in modeled_residues.json (case insensitive)

class DomainPredictor:

    def __init__(self):
        self.model = 0
        self.pae_cutoff = 5.0
        self.plddt_cutoff = 70
        self.resolution = 1.0
        self.pae_power = 1.0
        self.library = "igraph"
        self.modeled_regions = None
        self.predictions_dir = None
        self.result_dir = None
        self.model_path = None
        self.pae_path = None

    def set_modeled_regions(self, input_file):
        """Set the modeled regions for the protein

        Args:
            modeled_regions (dict): Dictionary of modeled regions
        """

        modeled_regions = read_json(input_file)
        self.modeled_regions = modeled_regions

    def set_model_path(self):
        """Set the path to the model file"""

        self.model_path = os.path.join(
            self.predictions_dir,
            self.job_name,
            f"fold_{self.job_name}_model_{self.model}.cif",
        )

        assert os.path.exists(
            self.model_path
        ), f"Model path {self.model_path} does not exist"

    def set_pae_path(self):
        """Set the path to the PAE file"""

        self.pae_path = os.path.join(
            self.predictions_dir,
            self.job_name,
            f"fold_{self.job_name}_full_data_{self.model}.json",
        )

        assert os.path.exists(self.pae_path), f"PAE path {self.pae_path} does not exist"

    def set_predictions_dir(self, path):
        """Set the path to the directory containing the predictions and PAE files
        Args:
            path (str): Path to the predictions directory
        """

        self.predictions_dir = path

    def set_result_dir(self, path):
        """Set the path to the directory to save the results

        Args:
            path (str): Path to the result directory
        """

        self.result_dir = path
        os.makedirs(self.result_dir, exist_ok=True)

    def set_output_type(self, output_type):
        """Set the output type for the domains
        set to "csv", "json", "txt" or "cxc"

        Args:
            output_type (str): Type of output file
        """

        assert output_type in ["csv", "json", "txt", "cxc"], "Invalid output type. Expected 'csv', 'json', 'txt' or 'cxc'"

        self.output_type = output_type

    def predict_domains(self):
        """Predict domains from a PAE file.
        Please Check pae_to_domains.py for more details

        Args:
            pae_file (str): Path to the PAE file

        Raises:
            ValueError: Invalid library specified. Use 'igraph' or 'network'

        Returns:
            domains (list): List of residues in each domain
        """

        pae_path = self.pae_path
        pae_matrix = parse_pae_file(pae_path)

        if self.library == "igraph":
            f = domains_from_pae_matrix_igraph

        elif self.library == "networkx":
            f = domains_from_pae_matrix_networkx

        else:
            raise ValueError("Invalid library specified. Use 'igraph' or 'network")

        domains = f(
            pae_matrix,
            pae_power=self.pae_power,
            pae_cutoff=self.pae_cutoff,
            graph_resolution=self.resolution,
        )

        return domains

    def high_plddt_residues(self):
        """Get the residues with pLDDT >= self.plddt_cutoff

        Args:
            mdoeled_regions (str): Path to the modeled residues file (json)

        Returns:
            confident_residues (list): List of residues with pLDDT >= plddt_cutoff
        """

        protein_name = self.job_name.split("_")[0].capitalize()
        self.protein_name = protein_name

        # the region of the protein to be modeled
        modeled_regions = self.modeled_regions
        region_of_interest = modeled_regions[protein_name]["region_to_be_modeled"]

        # get the pLDDT values for the modeled region
        model_path = self.model_path
        pae_path = self.pae_path
        p_chain = "A"

        prediction = AfParser(model_path, pae_path)
        plddt = prediction.get_ca_plddt()

        plddt[p_chain] = plddt[p_chain][
            region_of_interest[0] - 1 : region_of_interest[1]
        ]

        confident_residues = [
            i + 1 for i, p in enumerate(plddt[p_chain]) if p >= self.plddt_cutoff
        ]

        return confident_residues

    def filter_domains_by_plddt(self):
        """Filter the predicted domains by pLDDT values

        Args:
            domains (list): List of predicted domains
            confident_residues (list): List of residues with pLDDT >= plddt_cutoff

        Returns:
            filtered_domains (list): List of domains with only confident_residues
        """

        filtered_domains = []
        domains = self.domains
        confident_residues = self.confident_residues

        for domain in domains:
            f_domain = list(set(domain).intersection(set(confident_residues)))

            if len(f_domain) > 0:
                filtered_domains.append(f_domain)

        return filtered_domains

    def get_rigid_bodies(self):
        """Predict domains for all monomer models in the predictions directory
        the output will be saved in the result directory
        Depending on the output type, the domains will be saved in a csv, json, txt or cxc file
        - cxc files can be opened in ChimeraX as:
            "runscript /path/to/file.cxc"
        """

        output_type = self.output_type
        job_names = [
            x
            for x in os.listdir(self.predictions_dir)
            if os.path.isdir(os.path.join(self.predictions_dir, x))
        ]

        for job_name in job_names:

            self.job_name = job_name
            self.set_pae_path()
            self.set_model_path()

            self.confident_residues = self.high_plddt_residues()
            self.domains = self.predict_domains()
            self.filtered_domains = self.filter_domains_by_plddt()

            save_path = os.path.join(
                self.result_dir, f"{job_name}_model_{self.model}.{output_type}"
            )
            self.save_domains(save_path=save_path)

    def save_domains(self, save_path):
        """Save the predicted pLDDT filtered domains to a file

        Args:
            save_path (str): Path to save the domains

        if output_type is "txt" or "csv", the domains will be saved in the following format:
            1,2,3,5,6,7
            10,11,12,13,14
        where each line or row is a domain

        if output_type is "json", the domains will be saved as a list of lists

        if output_type is "cxc", the domains will be saved as a cxc file (ChimeraX commands)
        """

        confident_residues = self.confident_residues
        domains = self.domains
        filtered_domains = self.filtered_domains
        output_type = self.output_type

        file_ext = os.path.splitext(save_path)[1][1:]
        assert file_ext == output_type, f"File extension does not match file type. Expected {output_type}, got {file_ext}"

        if file_ext == "csv":
            with open(save_path, "wt") as outfile:
                for d in filtered_domains:
                    outfile.write(",".join([str(res) for res in sorted(d)]) + "\n")

        elif file_ext == "json":
            write_json(save_path, filtered_domains)

        elif file_ext == "txt":
            with open(save_path, "w") as f:
                for d in filtered_domains:
                    # d_range = get_key_from_res_range(d)
                    f.write(",".join([str(res) for res in sorted(d)]) + "\n")

        elif file_ext == "cxc":
            plddt_islands = get_key_from_res_range(confident_residues)
            with open(save_path, "w") as f:

                model_name = f"fold_{self.job_name}_model_{self.model}.cif"
                model_path = os.path.join(
                    self.predictions_dir, self.job_name, model_name
                )
                model_path = os.path.abspath(model_path)
                f.write(f"open {model_path}\n")

                for d in domains:
                    while True:
                        color = "%06x" % random.randint(0, 0xFFFFFF)
                        if not (color.startswith("ff") or color.startswith("FF")):
                            break
                    d_range = get_key_from_res_range(d)
                    f.write(f"col #1/A:{d_range} #{color}\n")

                if len(plddt_islands) > 0:
                    f.write(
                        f'sel #1/A:{self.modeled_regions[self.protein_name]["region_to_be_modeled"][0]}'
                        + f'-{self.modeled_regions[self.protein_name]["region_to_be_modeled"][1]}\n'
                    )
                    f.write("sel ~sel; hide sel r\n")
                    f.write(f"sel #1/A:{plddt_islands} ; sel ~sel ; col sel #FF0000\n")
                else:
                    print(
                        f"No confident regions (pLDDT >= {self.plddt_cutoff}) found for {self.job_name} in the given region"
                    )

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument(
        "--af_pipeline_path",
        type=str,
        required=True,
        default="../../af_pipeline",
    )
    args.add_argument(
        "--input_file",
        type=str,
        required=False,
        default="../../../../data/modeled_residues.json",
    )
    args.add_argument(
        "--predictions_dir",
        type=str,
        required=False,
        default="../../../../data/AF_predictions/AF_monomer",
    )
    args.add_argument(
        "--result_dir",
        type=str,
        required=False,
        default="../../../misc_data/domains"
    )
    args.add_argument(
        "--output_type",
        type=str,
        required=False,
        default="cxc"
    )

    args = args.parse_args()

    af_pipeline_path = os.path.abspath(args.af_pipeline_path)
    if os.path.exists(af_pipeline_path):
        sys.path.insert(0, af_pipeline_path)
        from parser import AfParser
    else:
        raise FileNotFoundError(f"Directory {af_pipeline_path} not found")

    DP = DomainPredictor()
    DP.set_modeled_regions(args.input_file)
    DP.set_predictions_dir(args.predictions_dir)
    DP.set_result_dir(args.result_dir)
    DP.set_output_type(args.output_type)
    DP.get_rigid_bodies()
