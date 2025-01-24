# NOTE:
# In the input yaml file, the following keys are required:
# job_cycle, entities, name, type
# job_cycle = [
#     {
#         "entities": [
#             {
#                 "name": "protein_1",
#                 "type": "proteinChain",
#             }
#         ]
#     }
# ]

# The remaining keys are optional and set to default values if not provided:
# job_name
#   format: a python string
#   default: see generate_job_name method
# modelSeeds
#   format: a list of integers
#   empty list
# range (what part of the sequence to use)
#   format: a tuple or a list of integers
#   default: (1, len(sequence)+1)
# count (how many copies of the entity)
#   format: an integer
#   default: 1
# glycans
#   format: a list of tuples where tuple = (residue, position)
#   default: empty list
# modifications
#   format: a list of tuples where tuple = (modificationType, basePosition) or (ptmType, ptmPosition)
#   default: empty list

from collections import defaultdict
import os
import json

# source: https://github.com/google-deepmind/alphafold/blob/main/server/README.md
allowed_modifications = {
    "proteinChain": [
        "CCD_SEP", 
        "CCD_TPO", 
        "CCD_PTR", 
        "CCD_NEP", 
        "CCD_HIP", 
        "CCD_ALY", 
        "CCD_MLY", 
        "CCD_M3L", 
        "CCD_MLZ", 
        "CCD_2MR", 
        "CCD_AGM", 
        "CCD_MCS", 
        "CCD_HYP", 
        "CCD_HY3", 
        "CCD_LYZ", 
        "CCD_AHB", 
        "CCD_P1L", 
        "CCD_SNN", 
        "CCD_SNC", 
        "CCD_TRF", 
        "CCD_KCR", 
        "CCD_CIR", 
        "CCD_YHA"
    ],
    "dnaSequence": [
        "CCD_5CM", 
        "CCD_C34", 
        "CCD_5HC", 
        "CCD_6OG", 
        "CCD_6MA", 
        "CCD_1CC", 
        "CCD_8OG", 
        "CCD_5FC", 
        "CCD_3DR"
    ],
    "rnaSequence": [
        "CCD_PSU", 
        "CCD_5MC", 
        "CCD_OMC", 
        "CCD_4OC", 
        "CCD_5MU", 
        "CCD_OMU", 
        "CCD_UR3", 
        "CCD_A2M", 
        "CCD_MA6", 
        "CCD_6MZ", 
        "CCD_2MG", 
        "CCD_OMG", 
        "CCD_7MG", 
        "CCD_RSQ"
    ]
}

allowed_small_molecules = {
    "ligand": [
        "CCD_ADP", 
        "CCD_ATP", 
        "CCD_AMP", 
        "CCD_GTP", 
        "CCD_GDP", 
        "CCD_FAD", 
        "CCD_NAD", 
        "CCD_NAP", 
        "CCD_NDP", 
        "CCD_HEM", 
        "CCD_HEC", 
        "CCD_PLM", 
        "CCD_OLA", 
        "CCD_MYR", 
        "CCD_CIT", 
        "CCD_CLA", 
        "CCD_CHL", 
        "CCD_BCL", 
        "CCD_BCB"
    ],
    "ion": [
        "MG", 
        "ZN", 
        "CL", 
        "CA", 
        "NA", 
        "MN", 
        "K", 
        "FE", 
        "CU", 
        "CO"
    ]
}

allowed_entity_types = [
    "proteinChain",
    "dnaSequence",
    "rnaSequence",
    "ligand",
    "ion",
]

class AFInput:
    """Class to create input files for AlphaFold prediction
    """

    def __init__(self, proteins: dict, protein_sequences: dict, input_yml: dict, nucleic_acid_sequences=None):
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.proteins = proteins
        self.input_yml = input_yml
        self.all_jobs = defaultdict(list)
        self.job_name = None
        self.output_dir = None

    def sanity_check_uniprot(self, proteins: dict):
        """Check if the proteins have uniprot assignments

        Args:
            proteins (dict): dictionary containing protein names and uniprot ids
        """

        uniprot_list = []
        for p_name in proteins:
            uniprot_id = self.proteins[p_name]
            if uniprot_id is not None:
                uniprot_list.append(uniprot_id)
            else:
                uniprot_list.append(p_name)

        # does protein_sequences.fasta contain all the uniprot ids?
        assert all(
            uniprot_id in self.protein_sequences
            for uniprot_id
            in uniprot_list
        ); "Uniprot ids not found in the fasta file"

    def write_to_json(self, sets_of_20: list, file_name: str):
        """Write the sets of 20 jobs to json files

        Args:
            sets_of_20 (list): list of sets of 20 jobs
            file_name (str): name of the file
        """

        for i, job_set in enumerate(sets_of_20):

            save_path = os.path.join(
                self.output_dir,
                f"{file_name}_set_{i}.json"
            )

            with open(save_path, "w") as f:
                json.dump(job_set, f, indent=4)

    def generate_job_name(self, entities: list):
        """Generate a job name from the entities

        If not provided in input yaml, the job name is generated
        from the entity names

        If the job name is too long, an exception is raised
        you can set the job name in the input yaml file
        or provide a shorter job name as an argument 'af_input.job_name'

        Args:
            entities (list): list of entities

        Raises:
            Exception: Job name is too long

        Returns:
            str: job name
        """

        if self.job_name is None:
            entity_names = [
                entity["name"]
                for entity
                in entities
            ]

            entity_ranges = []
            for entity in entities:
                if "range" in entity:
                    entity_ranges.append(
                        f"{entity['range'][0]}to{entity['range'][1]}"
                    )
                else:
                    entity_seq = self.get_entity_sequence(entity)
                    entity_ranges.append(
                        f"1to{len(entity_seq)}"
                    )

            job_name = "_".join(
                [
                    f"{name}_{range_}"
                    for name, range_
                    in zip(entity_names, entity_ranges)
                ]
            )

            if len(job_name) >= 100: # AF server limits job name to 100 characters
                raise Exception(
                    "Job name is too long, please provide a shorter job name"
                )

        else:
            job_name = self.job_name + "_" + str(len(self.all_jobs[self.job_cycle]))

        return job_name

    def add_small_molecule(self, sequences: list, entity: dict):
        """ Add a small molecule entity to the sequences (ion or ligand)

        whether an entity is a valid small molecule?
        if yes, the entity is added to the sequences

        Args:
            sequences (list): list of sequences (This is AlphaFold 'sequences')
            entity (dict): entity dictionary
        """

        small_molecule = entity["name"]
        entity_type = entity["type"]

        assert entity_type in ["ion", "ligand"]; "Invalid small molecule, must be ion or ligand"
        assert small_molecule in allowed_small_molecules[entity_type]; "Invalid small molecule"

        count = entity["count"]
        sequences.append(
                    {
                        entity_type: {
                            entity_type: small_molecule,
                            "count": count
                        }
                    }
                )

    def add_nucleic_acid_chain(self, sequences: list, entity: dict):
        """Add a nucleic_acid entity to the sequences (dnaSequence or rnaSequence)

        Whether an entity is a valid nucleic_acid?
        If yes, the entity is added to the sequences

        The modifications are added to the sequences if provided

        Args:
            sequences (list): list of sequences (this is AlphaFold 'sequences')
            entity (dict): entity dictionary

        Raises:
            Exception: Invalid modification type
        """

        nucleic_acid = entity["name"]
        entity_type = entity["type"]

        assert nucleic_acid in self.nucleic_acid_sequences; "Invalid nucleic_acid sequence or name"
        assert entity_type in ["dnaSequence", "rnaSequence"]; "Invalid entity type, must be dnaSequence or rnaSequence"

        count = entity["count"] if "count" in entity else 1
        # sequence = self.nucleic_acid_sequences[nucleic_acid]
        sequence = self.get_entity_sequence(entity)
        start, end = entity["range"] if "range" in entity else (1, len(sequence)+1)

        if start < 1 or end > len(sequence)+1:
            raise Exception(f"Invalid range for {nucleic_acid}: {start} to {end}; sequence length is {len(sequence)}")

        modifications = entity["modifications"] if "modifications" in entity else []

        if not all(
            [
            mod[0] in allowed_modifications[entity_type]
            for mod in modifications
            ]
        ):
            raise Exception("Invalid modification type")

        if len(modifications) > 0:
            modifications = [
                {
                    "modificationType": modification[0],
                    "basePosition": modification[1]-start+1,
                }
                for modification
                in modifications
                ]
        sequences.append(
                    {
                        entity_type: {
                            "sequence": sequence[start-1:end-1],
                            "count": count,
                            "modifications": modifications
                        }
                    }
                )

    def add_protein_chain(self, sequences, entity):
        """Add a protein chain entity to the sequences (proteinChain)

        Args:
            sequences (list): list of sequences (This is AlphaFold 'sequences')
            entity (dict): entity dictionary

        Raises:
            Exception: Invalid modification type
        """

        protein = entity["name"]
        entity_type = entity["type"]
        assert entity_type == "proteinChain"; "This is not a protein chain"

        count = entity["count"] if "count" in entity else 1
        # sequence = self.get_protein_sequence(entity)
        sequence = self.get_entity_sequence(entity)
        start, end = entity["range"] if "range" in entity else (1, len(sequence)+1)

        if start < 1 or end > len(sequence)+1:
            raise Exception(f"Invalid range for {protein}: {start} to {end}; sequence length is {len(sequence)}")

        glycans = entity["glycans"] if "glycans" in entity else []

        if len(glycans) > 0:
            glycans = [
                        {
                            "residue": glycan[0],
                            "position": glycan[1]-start+1,
                        }
                        for glycan
                        in glycans
                    ]
        modifications = entity["modifications"] if "modifications" in entity else []

        if not all(
            [
            mod[0]
            in allowed_modifications["proteinChain"]
            for mod in modifications
            ]
        ):
            raise Exception("Invalid modification type")

        if len(modifications) > 0:
            modifications = [
            {
                "ptmType": modification[0],
                "ptmPosition": modification[1]-start+1,
            }
            for modification
            in modifications
            ]

        sequences.append(
                    {
                        "proteinChain": {
                            "sequence": sequence[start-1:end-1],
                            "count": count,
                            "glycans": glycans,
                            "modifications": modifications
                        }
                    }
                )

    def get_entity_sequence(self, entity):
        """Get the sequence of the entity

        Args:
            entity (dict): entity dictionary

        Raises:
            Exception: Sequence not found

        Returns:
            str: sequence of the entity
        """

        if entity["type"] == "proteinChain":
            return self.get_protein_sequence(entity)

        elif entity["type"] == "dnaSequence" or entity["type"] == "rnaSequence":
            return self.nucleic_acid_sequences[entity["name"]]

        else:
            raise Exception("Sequence not found")

    def get_protein_sequence(self, entity):
        """Get the protein sequence from the entity dictionary

        Args:
            entity (dict): entity dictionary

        Raises:
            Exception: Protein sequence not found

        Returns:
            str: protein sequence
        """

        assert entity["type"] == "proteinChain"; "This is not a protein chain"
        uniprot_id = self.proteins[entity["name"]]

        # for unknown proteins, the uniprot id is set to None in the proteins.json file
        if uniprot_id is not None:
            sequence = self.protein_sequences[uniprot_id]

        else:
            if entity["name"] not in self.protein_sequences:
                raise Exception("Protein sequence not found")

            # if the protein is novel, add the sequence directly to
            # protein_sequences.fasta with the protein name as the header
            sequence = self.protein_sequences[entity["name"]]

        return sequence

    def make_sequences(self, entities):
        """Make sequences from the entities

        Args:
            entities (list): list of entities

        Raises:
            Exception: Invalid entity type

        Returns:
            list: list of sequences (This is AlphaFold 'sequences')
        """

        sequences = []
        for entity in entities:

            entity_type = entity["type"]
            if entity_type not in allowed_entity_types:
                raise Exception("Invalid entity type")

            if entity_type == "proteinChain":
                self.add_protein_chain(sequences, entity)

            elif entity_type == "dnaSequence" or entity_type == "rnaSequence":
                self.add_nucleic_acid_chain(sequences, entity)

            elif entity_type == "ligand" or entity_type == "ion":
                self.add_small_molecule(sequences, entity)

        return sequences

    def make_job(self, job_cycle):
        """Make a single AF-server job

        Args:
            job_cycle (str): job type

        Raises:
            Exception: No entities were found in the job
        """

        jobs = self.input_yml[job_cycle]
        for job in jobs:
            if "entities" not in job:
                raise Exception("No entities were found in the job")

            job_name = job["name"] if "name" in job else self.generate_job_name(job["entities"])
            model_seeds = job["modelSeeds"] if "modelSeeds" in job else []
            sequences = self.make_sequences(job["entities"])

            self.all_jobs[job_cycle].append(
                {
                    "name": job_name,
                    "modelSeeds": model_seeds,
                    "sequences": sequences
                }
            )

    def create_job_files(self):
        """Create job files for AlphaFold prediction
        for each job cycle, make_job is called to create the jobs that are
        stored in 'all_jobs'

        all_jobs = {
            "job_cycle_1": [
                {
                    "name": "job_name",
                    "modelSeeds": [list of model seeds],
                    "sequences": [list of sequences]
                },
                ...
            ],
            ...
        }

        for each job cycle, the jobs are divided into sets of 20 and written
        to json files

        Raises:
            Exception: Output directory not set
        """

        if self.output_dir is None:
            raise Exception("Output directory not set")

        for job_cycle, jobs in self.input_yml.items():
            self.job_cycle = job_cycle
            self.make_job(job_cycle)
            print(f"Created jobs for {job_cycle}: {len(self.all_jobs[job_cycle])}")

        for job_cycle, jobs in self.all_jobs.items():
            sets_of_20 = [
                jobs[i:i+20]
                for i
                in range(0, len(jobs), 20)
            ]

            os.makedirs(self.output_dir, exist_ok=True)
            self.write_to_json(sets_of_20, file_name=job_cycle)