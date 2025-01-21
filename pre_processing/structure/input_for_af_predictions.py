# Description: write alphafold json file for prediction from target protein pairs
# Input: yaml file containing the target proteins
# Output: json files for alphafold prediction
# Use case:
# you have a large number of proteins for which you want to predict the structure (monomer or complex)
# so, the number of jobs for alphafold prediction is large
# This script will create the input json files for AF server (20 jobs per file)

#TODO: better way to define the input to the script or the input dictionary can be generated programmatically

#NOTE: currently, it is helpful to have copilot do the job of defining the input
#      (once you start typing the dictionary, copilot will suggest the rest of the dictionary)

import yaml
import argparse
import sys
sys.path.append("../")
import os
import json
from utils import read_json, read_fasta

class AFInput:
    """Create input for AlphaFold prediction
    """

    def __init__(self, to_be_predicted):
        self.to_be_predicted = to_be_predicted
        self.all_sequences = read_fasta("../output/all_sequences.fasta")
        self.proteins = read_json("../inputs/cardiac_desmosome_proteins.json")
        self.output_dir = "../output/af_input"


    def get_uniprot_ids(self, proteins):
        """Get the uniprot ids for the proteins in the list

        Args:
            proteins (list): list of protein names

        Returns:
            list: list of uniprot ids
        """

        assert all(
            p_name in self.proteins
            for p_name
            in proteins
        ), "Missing UniProt assignments for the proteins"

        uniprot_list = [
            self.proteins[p_name]
            for p_name
            in proteins
        ]

        assert all(
            uniprot_id in self.all_sequences
            for uniprot_id
            in uniprot_list
        ), "Uniprot ids not found in the fasta file"

        return uniprot_list


    def generate_job_name(self, proteins, kwargs):
        """Generate a job name for the given protein list and kwargs describing the job

        Args:
            proteins (list): list of protein names
            kwargs (dict): dictionary of job description

        Returns:
            str: job name
        """
        # add other arguments in the valid list as needed (for example, "seeds")
        assert all(
            arg in ["copies", "fragments"]
            for arg
            in kwargs.keys()
        ), "Invalid argument in kwargs"

        copy_list = kwargs.get("copies", [1] * len(proteins))
        uniprot_list = self.get_uniprot_ids(proteins)
        fragment_list = kwargs.get("fragments", [(1, len(self.all_sequences[uniprot_id])) for uniprot_id in uniprot_list])
        copy_list = [str(copy) for copy in copy_list]
        fragment_list = [f"{str(fragment[0])}to{str(fragment[1])}" for fragment in fragment_list]

        job_name = (
            " ".join(
                [
                    f"{p_name}-{copy}-{fragment}"
                    for p_name, copy, fragment
                    in zip(proteins, copy_list, fragment_list)
                ]
            )
        )
        return job_name


    def write_to_json(self, sets_of_20, file_name ,output_dir="../output/af_input"):
        """Write the sets of 20 jobs to json files

        Args:
            sets_of_20 (list): list of sets of 20 jobs
        """

        for i, a_set in enumerate(sets_of_20):
            with open(os.path.join(output_dir, f"{file_name}_{i}.json"), "w") as f:
                json.dump(a_set, f, indent=4)


    def make_job(self, job_name, combined_info):
        """Create a dictionary for a job

        Args:
            job_name (str): name of the job
            combined_info (list): combined list containing protein, uniprot_id, copy, fragment

        Returns:
            dict: dictionary for the job
        """

        proteins, uniprot_list, copy_list, fragment_list = zip(*combined_info)
        af_dict = {
                "name": job_name,
                "modelSeeds": [],
                "sequences": []
            }

        for idx, uniprot_id in enumerate(uniprot_list):
            start = fragment_list[idx][0] - 1
            end = fragment_list[idx][1]
            seq = self.all_sequences[uniprot_id]
            af_dict["sequences"].append(
                    {
                        "proteinChain": {
                            "sequence": seq[start:end],
                            "count": copy_list[idx]
                        }
                    }
                )

        return af_dict


    def sanity_check(self, proteins, uniprot_list, copy_list, fragment_list):
        """Check if the input lists are of same length

        Args:
            proteins (list): list of protein names
            uniprot_list (list): list of uniprot ids
            copy_list (list): list of copy numbers
            fragment_list (list): list of sequence fragments (start, end)
        """

        assert (
                len(proteins) == len(uniprot_list) == len(copy_list) == len(fragment_list)
            ), "proteins, uniprot_list, copy_list and fragment_list should be of same length"


    def create_af_input(self, file_name="af_input", job_name=None):
        """Create input for AlphaFold prediction
        """

        af_dict_list = []  # all jobs for AF server

        for jobs in self.to_be_predicted:
            proteins = jobs["proteins"]
            uniprot_list = self.get_uniprot_ids(proteins)

            copy_list = jobs.get(
                "copies",
                [1] * len(proteins)
            )

            # job_name = self.generate_job_name(
            #     proteins,
            #     {"copies": copy_list}
            # )

            fragment_list = jobs.get(
                "fragments",
                [
                    (1, len(self.all_sequences[uniprot_id]))
                    for uniprot_id
                    in uniprot_list
                ]
            )

            # you can add more arguments in the kwargs as needed
            if job_name is None:
                job_name = self.generate_job_name(
                    proteins,
                    {
                        "copies": copy_list,
                        "fragments": fragment_list
                    }
                )
            else:
                job_name = job_name

            self.sanity_check(proteins, uniprot_list, copy_list, fragment_list)

            combined_info = zip(proteins, uniprot_list, copy_list, fragment_list)

            af_dict = self.make_job(job_name, combined_info) # a single job for AF server

            af_dict_list.append(af_dict)

        # AF server at a time accepts 20 jobs; so divide the list into sets of 20
        sets_of_20 = [
            af_dict_list[i:i+20]
            for i
            in range(0, len(af_dict_list), 20)
        ]

        os.makedirs(self.output_dir, exist_ok=True)
        self.write_to_json(sets_of_20, file_name=file_name)


if __name__ == "__main__":

    args = argparse.ArgumentParser()
    args.add_argument(
        "-s",
        "--sequences",
        type=str,
        required=False,
        default="../output/all_sequences.fasta",
        help="fasta file containing all sequences"
    )
    args.add_argument(
        "-p",
        "--proteins",
        type=str,
        required=False,
        default="../inputs/cardiac_desmosome_proteins.json",
        help="json file containing protein names and uniprot ids"
    )
    args.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="../output/af_input",
        help="output directory for alphafold input"
    )
    args.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        default="../inputs/af_server_targets.yaml",
        help="input yaml file containing the target proteins"
    )
    args = args.parse_args()

    monomer_targets = yaml.load(open(args.input), Loader=yaml.FullLoader)["monomer_targets"]
    pairwise_targets = yaml.load(open(args.input), Loader=yaml.FullLoader)["pairwise_targets"]
    full_odp_targets = yaml.load(open(args.input), Loader=yaml.FullLoader)["full_odp_targets"]

    monomer_input = AFInput(monomer_targets)
    monomer_input.proteins = read_json(args.proteins)
    monomer_input.all_sequences = read_fasta(args.sequences)
    monomer_input.output_dir = args.output

    pairwise_input = AFInput(pairwise_targets)
    pairwise_input.proteins = read_json(args.proteins)
    pairwise_input.all_sequences = read_fasta(args.sequences)
    pairwise_input.output_dir = args.output

    full_odp_input = AFInput(full_odp_targets)
    full_odp_input.proteins = read_json(args.proteins)
    full_odp_input.all_sequences = read_fasta(args.sequences)
    full_odp_input.output_dir = args.output

    # write the input json for alphafold prediction
    monomer_input.create_af_input(file_name="af_input_monomer")
    pairwise_input.create_af_input(file_name="af_input_pairwise")
    full_odp_input.create_af_input(file_name="af_input_full_odp", job_name="full_odp")