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
#   format: a list of integers
#   default: [1, len(sequence)+1]
# count (how many copies of the entity)
#   format: an integer
#   default: 1
# glycans
#   format: a list of lists where inner_list = (residue, position)
#   default: empty list
# modifications
#   format: a list of lists where inner_list = (modificationType, basePosition) or (ptmType, ptmPosition)
#   default: empty list

import os
import json
import random
from utils import read_json

allowed_things = read_json(
    os.path.join(os.path.dirname(__file__), "allowed_af_things.json")
)
allowed_modifications = allowed_things["allowed_modifications"]
allowed_small_molecules = allowed_things["allowed_small_molecules"]
allowed_entity_types = allowed_things["allowed_entity_types"]


class AFInput:

    def __init__(
        self,
        input_yml: dict,
        protein_sequences: dict,
        nucleic_acid_sequences: dict | None = None,
        proteins: dict = {},
    ):
        self.input_yml = input_yml
        self.proteins = proteins
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.output_dir = None
        self.job_cycles = {}

    def create_job_cycles(self):
        """Create job cycles from the input yaml file"""

        for job_cycle, jobs_info in self.input_yml.items():
            print("Creating job cycle", job_cycle, "\n")
            af_cycle = AFCycle(
                jobs_info=jobs_info,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                proteins=self.proteins,
            )

            af_cycle.update_cycle()
            self.job_cycles[job_cycle] = af_cycle.cycle

    def write_to_json(self, sets_of_20: list, file_name: str):
        """Write the sets of 20 jobs to json files

        Args:
            sets_of_20 (list): list of sets of 20 jobs
            file_name (str): name of the file
        """

        for i, job_set in enumerate(sets_of_20):

            save_path = os.path.join(self.output_dir, f"{file_name}_set_{i}.json")

            with open(save_path, "w") as f:
                json.dump(job_set, f, indent=4)
            print(f"{len(job_set)} jobs written for {file_name}_set_{i}")

    def write_job_files(self):
        """Write job files to the output directory"""

        for job_cycle, jobs in self.job_cycles.items():
            sets_of_20 = [jobs[i : i + 20] for i in range(0, len(jobs), 20)]

            os.makedirs(self.output_dir, exist_ok=True)
            self.write_to_json(sets_of_20, file_name=job_cycle)
        print("\nAll job files written to", self.output_dir)


class AFCycle:

    def __init__(
        self,
        jobs_info: list,
        protein_sequences: dict,
        nucleic_acid_sequences: dict | None = None,
        proteins: dict = {},
    ):
        self.jobs_info = jobs_info
        self.proteins = proteins
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.cycle = []

    def update_cycle(self):
        """Create jobs for the cycle"""

        for job in self.jobs_info:
            af_job = AFJob(
                job_info=job,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                proteins=self.proteins,
            )
            af_job.add_job_name()
            af_job.add_model_seeds()
            af_job.add_sequences()

            if af_job.job_name is None:
                af_job.generate_job_name()

            af_job.update_job()
            self.add_jobs(af_job.job)

    def add_jobs(self, job: dict):
        """Seed the jobs"""
        if len(job["modelSeeds"]) == 0:
            self.cycle.append(job)

        else:
            for seed in job["modelSeeds"]:
                job_copy = job.copy()
                job_copy["modelSeeds"] = [seed]
                self.cycle.append(job_copy)


class AFJob:

    def __init__(
        self,
        job_info: dict,
        protein_sequences: dict,
        nucleic_acid_sequences: dict | None = None,
        proteins: dict = {},
    ):
        self.job_info = job_info
        self.proteins = proteins
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.job = {}
        self.job_name = None
        self.model_seeds = []
        self.sequences = []
        self.name_fragments = []

    def add_job_name(self):
        """Create a job from the job info"""

        if "name" in self.job_info:
            self.job_name = self.job_info["name"]

    def add_model_seeds(self):
        """Update the model seeds"""

        if "modelSeeds" in self.job_info:

            if isinstance(self.job_info["modelSeeds"], int):
                self.generate_seeds()

            elif isinstance(self.job_info["modelSeeds"], list):
                self.model_seeds = self.job_info["modelSeeds"]

            else:
                raise Exception("modelSeeds must be an integer or a list")

    def generate_seeds(self):
        """Generate model seeds"""

        assert isinstance(
            self.job_info["modelSeeds"], int
        ), "modelSeeds must be an integer"

        num_seeds = self.job_info["modelSeeds"]
        self.model_seeds = random.sample(range(1, 10 * num_seeds), num_seeds)

    def add_sequences(self):
        """Update the sequences"""

        for entity_info in self.job_info["entities"]:
            entity = AFSequence(
                entity_info=entity_info,
                protein_sequences=self.protein_sequences,
                nucleic_acid_sequences=self.nucleic_acid_sequences,
                proteins=self.proteins,
            )
            entity.update_real_sequence()
            entity.update_sequence()
            self.name_fragments.append(entity.get_name_fragment())
            self.sequences.append(entity.sequence)

    def generate_job_name(self):
        """Generate a job name"""

        job_name = "_".join(self.name_fragments)
        self.job_name = job_name

    def update_job(self):
        """Update the job"""

        self.job = {
            "name": self.job_name,
            "modelSeeds": self.model_seeds,
            "sequences": self.sequences,
        }


class Entity:

    def __init__(
        self,
        entity_info: dict,
        protein_sequences: dict,
        nucleic_acid_sequences: dict,
        proteins: dict = {},
    ):
        self.entity_info = entity_info
        self.proteins = proteins
        self.protein_sequences = protein_sequences
        self.nucleic_acid_sequences = nucleic_acid_sequences
        self.entity_name = entity_info["name"]
        self.entity_type = entity_info["type"]
        self.real_sequence = None
        self.start = 1
        self.end = None
        self.glycans = None
        self.modifications = None

    def get_real_sequence(self):
        """Get the protein or dna or rna sequence of the entity"""

        if self.entity_type == "proteinChain":

            try:
                uniprot_id = self.proteins[self.entity_name]
                real_sequence = self.protein_sequences[uniprot_id]

            except KeyError:
                try:
                    real_sequence = self.protein_sequences[self.entity_name]
                except:
                    raise Exception(
                        f"Could not find the entity sequence for {self.entity_name}"
                    )

        elif self.entity_type == "dnaSequence":
            real_sequence = self.nucleic_acid_sequences[self.entity_name]

        elif self.entity_type == "rnaSequence":
            real_sequence = self.nucleic_acid_sequences[self.entity_name]

        else:
            real_sequence = None

        self.real_sequence = real_sequence
        return self.real_sequence

    def get_entity_range(self):
        """Update the start of the entity"""

        if "range" in self.entity_info:

            assert (
                len(self.entity_info["range"]) == 2
            ), "Invalid range; must be a list of two integers (start and end)"

            start, end = self.entity_info["range"]

        elif self.real_sequence is not None:
            start, end = 1, len(self.real_sequence)

        else:
            start, end = 1, 1

        self.start, self.end = start, end
        return self.start, self.end

    def get_glycans(self):
        """Get the glycans of the entity"""

        if "glycans" in self.entity_info and self.entity_type == "proteinChain":
            glycans = self.entity_info["glycans"]
            glycans = [
                {
                    "residues": glycan[0],
                    "position": glycan[1] - self.start + 1,
                }
                for glycan in glycans
            ]

        else:
            glycans = []

        self.glycans = glycans
        self.sanity_check_glycans()

        return self.glycans

    def sanity_check_glycans(self):
        """Sanity check the glycans"""

        assert (
            self.entity_type == "proteinChain"
        ), "Glycans are only supported for protein chains"

        for glycan in self.glycans:
            glyc_pos = glycan["position"]

            if glyc_pos < 1 or glyc_pos > len(self.real_sequence):
                raise Exception(
                    f"Invalid glycan position at {glyc_pos} in {self.entity_name}"
                )

    def get_modifications(self):
        """Get the modifications of the entity"""

        if "modifications" in self.entity_info:
            modifications = self.entity_info["modifications"]

            if self.entity_type == "proteinChain":
                modifications = [
                    {
                        "ptmType": mod[0],
                        "ptmPosition": mod[1] - self.start + 1,
                    }
                    for mod in modifications
                ]

            elif self.entity_type == "dnaSequence" or self.entity_type == "rnaSequence":
                modifications = [
                    {
                        "modificationType": mod[0],
                        "basePosition": mod[1] - self.start + 1,
                    }
                    for mod in modifications
                ]

            else:
                raise Exception("modifications are not supported for this entity type")

        else:
            modifications = []

        self.modifications = modifications
        self.sanity_check_modifications()

        return modifications

    def get_small_molecule(self):
        """Get the small molecule of the entity"""

        if self.entity_type in ["ligand", "ion"]:
            small_molecule = self.entity_name

        else:
            small_molecule = None

        self.sanity_check_small_molecule()
        return small_molecule

    def sanity_check_modifications(self):
        """Sanity check the modifications"""

        assert self.entity_type in [
            "proteinChain",
            "dnaSequence",
            "rnaSequence",
        ], "Modifications are only supported for protein chains, dna and rna sequences"

        if self.entity_type == "proteinChain":
            if not all(
                [
                    mod["ptmType"] in allowed_modifications["proteinChain"]
                    for mod in self.modifications
                ]
            ):
                raise Exception("Invalid modification type")

        elif self.entity_type == "dnaSequence" or self.entity_type == "rnaSequence":
            if not all(
                [
                    mod["modificationType"] in allowed_modifications[self.entity_type]
                    for mod in self.modifications
                ]
            ):
                raise Exception("Invalid modification type")

        for mod in self.modifications:
            mod_pos = (
                mod["ptmPosition"]
                if self.entity_type == "proteinChain"
                else mod["basePosition"]
            )

            if mod_pos < 1 or mod_pos > len(self.real_sequence):
                raise Exception(
                    f"Invalid modification at {mod_pos} in {self.entity_name}"
                )

    def sanity_check_small_molecule(self):
        """Sanity check the small molecule"""

        assert self.entity_type in [
            "ligand",
            "ion",
        ], "Invalid entity type for small molecule"

        if self.entity_name not in allowed_small_molecules[self.entity_type]:
            raise Exception("Invalid small molecule")

        else:
            pass


class AFSequence(Entity):

    def __init__(
        self,
        entity_info: dict,
        protein_sequences: dict,
        nucleic_acid_sequences: dict,
        proteins: dict = {},
    ):
        super().__init__(
            entity_info=entity_info,
            protein_sequences=protein_sequences,
            nucleic_acid_sequences=nucleic_acid_sequences,
            protiens=proteins,
        )
        self.name = self.entity_name
        self.type = self.entity_type
        self.real_sequence = self.update_real_sequence()
        self.count = self.add_entity_count()
        self.sequence = {}

    def update_real_sequence(self):
        """Update the real sequence of the entity"""
        self.get_real_sequence()
        self.get_entity_range()
        if self.type in ["proteinChain", "dnaSequence", "rnaSequence"]:
            self.real_sequence = self.real_sequence[self.start - 1 : self.end]

    def add_entity_count(self):
        """Update the count of the entity"""
        if "count" in self.entity_info:
            count = self.entity_info["count"]
        else:
            count = 1
        return count

    def get_name_fragment(self):
        """Get the name fragments of the entity"""
        return f"{self.name}_{self.count}_{self.start}to{self.end}"

    def update_sequence(self):
        """Update the sequence of the entity"""
        if self.type == "proteinChain":
            self.sequence = {
                self.type: {
                    "sequence": self.real_sequence,
                    "glycans": self.get_glycans(),
                    "modifications": self.get_modifications(),
                    "count": self.count,
                }
            }
        elif self.type in ["dnaSequence", "rnaSequence"]:
            self.sequence = {
                self.type: {
                    "sequence": self.real_sequence,
                    "modifications": self.get_modifications(),
                    "count": self.count,
                }
            }
        elif self.type in ["ligand", "ion"]:
            self.sequence = {
                self.type: {self.type: self.get_small_molecule(), "count": self.count}
            }
