# Description: write alphafold json file for prediction from target gene pairs
# Input: to be defined in the main function
# Output: json files for alphafold prediction

#TODO: better way to define the input to the script

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


    def get_uniprot_ids(self, protein_tuple):
        """Get the uniprot ids for the proteins in the tuple

        Args:
            protein_tuple (tuple): tuple of protein names

        Returns:
            list: list of uniprot ids
        """

        assert all(
            p_name in self.proteins
            for p_name
            in protein_tuple
        ), "Missing UniProt assignments for the proteins"

        uniprot_list = [
            self.proteins[p_name]
            for p_name
            in protein_tuple
        ]

        assert all(
            uniprot_id in self.all_sequences
            for uniprot_id
            in uniprot_list
        ), "Uniprot ids not found in the fasta file"

        return uniprot_list


    def generate_job_name(self, protein_tuple, kwargs):
        """Generate a job name for the given protein tuple and kwargs describing the job

        Args:
            protein_tuple (tuple): tuple of protein names
            kwargs (dict): dictionary of job description

        Returns:
            str: job name
        """

        assert all(
            arg in ["copies", "fragments"]
            for arg
            in kwargs.keys()
        ), "Invalid argument in kwargs"

        copy_list = kwargs.get("copies", [1] * len(protein_tuple))
        uniprot_list = self.get_uniprot_ids(protein_tuple)
        fragment_list = kwargs.get("fragments", [(1, len(self.all_sequences[uniprot_id])) for uniprot_id in uniprot_list])
        copy_list = [str(copy) for copy in copy_list]
        fragment_list = [f"{str(fragment[0])}to{str(fragment[1])}" for fragment in fragment_list]

        job_name = (
            " ".join(
                [
                    f"{p_name}-{copy}-{fragment}"
                    for p_name, copy, fragment
                    in zip(protein_tuple, copy_list, fragment_list)
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

        protein_tuple, uniprot_list, copy_list, fragment_list = zip(*combined_info)
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


    def sanity_check(self, protein_tuple, uniprot_list, copy_list, fragment_list):
        """Check if the input lists are of same length

        Args:
            protein_tuple (tuple): tuple of protein names
            uniprot_list (list): list of uniprot ids
            copy_list (list): list of copy numbers
            fragment_list (list): list of sequence fragments (start, end)
        """

        assert (
                len(protein_tuple) == len(uniprot_list) == len(copy_list) == len(fragment_list)
            ), "protein_tuple, uniprot_list, copy_list and fragment_list should be of same length"


    def create_af_input(self, file_name="af_input"):
        """Create input for AlphaFold prediction
        """

        af_dict_list = []  # all jobs for AF server

        for protein_tuple, input_info in self.to_be_predicted.items():

            uniprot_list = self.get_uniprot_ids(protein_tuple)

            copy_list = input_info.get(
                "copies",
                [1] * len(protein_tuple)
            )

            job_name = self.generate_job_name(
                protein_tuple,
                {"copies": copy_list}
            )

            fragment_list = input_info.get(
                "fragments",
                [
                    (1, len(self.all_sequences[uniprot_id]))
                    for uniprot_id
                    in uniprot_list
                ]
            )

            job_name = self.generate_job_name(
                protein_tuple,
                {
                    "copies": copy_list,
                    "fragments": fragment_list
                }
            )

            self.sanity_check(protein_tuple, uniprot_list, copy_list, fragment_list)

            combined_info = zip(protein_tuple, uniprot_list, copy_list, fragment_list)

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

    # pairwise prediction
    pairwise_targets = {
        ("Dsc2a", "Dsg2"): {
            "copies": [1, 1]
        },
        ("Dsc2a", "Pkp2a"): {
            "copies": [1, 1]
        },
        ("Dsc2a", "Dp1"): {
            "copies": [1, 1]
        },
        ("Dsc2a", "Pg"): {
            "copies": [1, 1]
        },
        ("Pg", "Pkp2a"): {
            "copies": [1, 1]
        },
        ("Pg", "Dsg2"): {
            "copies": [1, 1]
        },
        ("Pg", "Dp1"): {
            "copies": [1, 1]
        },
        ("Dp1", "Dsg2"): {
            "copies": [1, 1]
        },
        ("Dp1", "Pkp2a"): {
            "copies": [1, 1]
        },
        ("Dsg2", "Pkp2a"): {
            "copies": [1, 1]
        },
    }

    # monomer prediction
    monomer_targets = {
        ("Dsc2a",): {
            "copies": [1]
        },
        ("Dsg2",): {
            "copies": [1]
        },
        ("Pg",): {
            "copies": [1]
        },
        ("Pkp2a",): {
            "copies": [1]
        },
        ("Dp1",): {
            "copies": [1]
        }
    }

    # full ODP prediction
    full_odp_targets = {
        ("Dsc2a", "Dsg2", "Pg", "Pkp2a", "Dp1"): {
            "copies": [1, 1, 2, 2, 2],
            "fragments": [(718, 901), (634, 941), (1, 745), (1, 837), (1, 584)]
        }
    }

    # write jobs to json
    AFInput(monomer_targets).create_af_input(file_name="af_input_monomer")
    AFInput(pairwise_targets).create_af_input(file_name="af_input_pairwise")
    AFInput(full_odp_targets).create_af_input(file_name="af_input_full_odp")