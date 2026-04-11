import os
import psa
import warnings
from typing import Any
from IMP_Toolbox.utils.file_helpers import read_fasta
from IMP_Toolbox.constants.sequence_constants import (
    PSAProgram,
    PSAAttribute,
    PSATerm,
    MolType,
)
from IMP_Toolbox.constants.imp_toolbox_constants import (
    FileFormat,
    MiscStrEnum,
)
from IMP_Toolbox.sequence.sequence import (
    fasta_str_to_dict,
)

class PairwiseSequenceAlignment:
    """ Handle pairwise sequence alignment from two sequences. """

    seq1: str
    """ First sequence. """

    seq2: str
    """ Second sequence. """

    program: PSAProgram
    """ Program to use for pairwise sequence alignment. """

    moltype: str = MolType.PROT.value
    """ Molecule type for alignment. Default is "protein". """

    pairwise_alignment: psa.PairwiseAlignment | None
    """ Pairwise alignment object. """

    def __init__(self, seq1: str, seq2: str, moltype: str = MolType.PROT.value):
        self.seq1 = seq1
        self.seq2 = seq2
        self.program = PSAProgram.STRETCHER
        self.moltype = moltype
        self.pairwise_alignment = None

    def perform_alignment(self):
        """ Perform pairwise sequence alignment using psa module.

        ## Returns:

        - **psa.PairwiseAlignment**:<br />
            Pairwise alignment object
        """

        self.pairwise_alignment = psa.align(
            program=self.program,
            moltype=self.moltype,
            qseq=self.seq1,
            sseq=self.seq2,
        )

        return self.pairwise_alignment

    def fetch_pairwise_alingment_map(
        self,
        pairwise_alignment_file: str | None = None,
        overwrite: bool = False,
    ) -> dict:
        """ Fetch the pairwise alignment map.

        > [!NOTE]
        > If `perform_alingment` has been called already, the pairwise alignment
        > map will be fetched from the `pairwise_alignment` attribute and
        > pairwise_alignment_file will be ignored if `overwrite` is False.

        ## Arguments:

        - **pairwise_alignment_file (str, optional):**:<br />
            Path to the pairwise alignment file. If not provided, the pairwise
            alignment must have been performed already.

        - **overwrite (bool, optional):**:<br />
            Whether to overwrite the pairwise alignment file if it already exists.
            Defaults to False.

        ## Returns:
        - **dict**:<br />
            Mapping of codon numbers to residue numbers
        """

        alignment_file_exists = (
            pairwise_alignment_file is not None and
            os.path.exists(pairwise_alignment_file)
        )

        if alignment_file_exists is False or overwrite is True:

            if self.alignment_performed is False:
                self.perform_alignment()

            psa_map = self.get_pairwise_alignment_map(
                pairwise_alignment=self.pairwise_alignment,
            )

            self.save_pairwise_alignment(
                save_path=pairwise_alignment_file,
                overwrite=overwrite,
            )

        elif alignment_file_exists is True and overwrite is False:

            if os.path.exists(pairwise_alignment_file) and overwrite is False:
                psa_map = self.get_pairwise_alignment_map(
                    pairwise_alignment_file=pairwise_alignment_file,
                )

        else:
            raise ValueError(
                """Invalid combination of pairwise_alignment_file and overwrite arguments.
                Either provide a valid pairwise_alignment_file or
                set overwrite to True to perform alignment and save it to the specified file or
                call `perform_alignment` first and then `fetch_pairwise_alignment_map`.
                """
            )

        return psa_map

    @property
    def alignment_performed(self) -> bool:
        """ Check if the pairwise alignment is performed.

        ## Returns:

        - **bool**:<br />
            True if the pairwise alignment is performed, False otherwise.
        """

        return self.pairwise_alignment is not None

    def verify_alignment_performed(self):
        """ Verify if the the pariwise alignment is performed else raise Exception."""

        if self.alignment_performed is False:
            raise Exception("Please perform alignment before accessing it")

    def verify_alignment_attribute(self, attribute: str):
        """ Verify if the attribute is available in PairwiseAlignment object

        ## Arguments:

        - **attribute (str)**:<br />
            Attribute to check in the PairwiseAlignment object
        """

        self.verify_alignment_performed()

        if hasattr(self.pairwise_alignment, attribute) is False:
            raise Exception(
                f"Attribute {attribute} not found in pairwise alignment object"\
                f"Available attributes are: "\
                f"{list(PSAAttribute)}"
            )

    def save_pairwise_alignment(self, save_path: str, overwrite: bool = False):
        """ Save pairwise alignment to the specified path.

        ## Arguments:

        - **save_path (str)**:<br />
            Path to save the pairwise alignment in fasta format. Should end with
            .fasta or .fa.

        - **overwrite (bool, optional):**:<br />
            Whether to overwrite the file if it already exists. Defaults to False.
        """

        ext = os.path.splitext(save_path)[1].lstrip('.')
        assert ext in [FileFormat.FASTA, FileFormat.FA], (
            f"Expected file extension to be .fasta or .fa, but got {ext}"
        )

        if not os.path.exists(save_path) or overwrite:

            self.verify_alignment_performed()
            os.makedirs(os.path.dirname(save_path),exist_ok=True)

            with open(save_path, "w") as f:
                f.write(self.pairwise_alignment.fasta())

    def get_alignment_attribute(self, attribute: str) -> Any:
        """ Get value of an attribute from PairwiseAlingnment object.

        ## Arguments:

        - **attribute (str)**:<br />
            Attribute to get from the PairwiseAlignment object. Should be one of
            the attributes defined in PSAAttribute enum.

        ## Returns:

        - **Any**:<br />
            Value of the requested attribute from the PairwiseAlignment object.
        """

        self.verify_alignment_performed()
        self.verify_alignment_attribute(attribute)

        return getattr(self.pairwise_alignment, attribute)

    @staticmethod
    def get_pairwise_alignment_map(
        pairwise_alignment_file: str | None = None,
        pairwise_alignment: psa.PairwiseAlignment | None = None,
    ) -> dict:
        """ Given a pairwise alignment, return one-one map of aligned residues.

        ## Arguments:

        - **pairwise_alignment_file (str)**:<br />
            Path to the pairwise alignment file

        ## Returns:

        - **dict**:<br />
            Mapping of codon numbers to residue numbers
        """

        if (
            (pairwise_alignment is None and pairwise_alignment_file is None) or
            (pairwise_alignment is not None and pairwise_alignment_file is not None)
        ):
            raise ValueError(
                """
                Only one of the following should be provided:
                - pairwise_alignment_file (path to saved pairwise alignment fasta)
                - pairwise_alignment (psa.PairwiseAlignment object)
                """
            )

        if isinstance(pairwise_alignment, psa.PairwiseAlignment):
            pairwise_alignment_dict = fasta_str_to_dict(
                fasta_str=pairwise_alignment.fasta()
            )

        if isinstance(pairwise_alignment_file, str):
            pairwise_alignment_dict = read_fasta(pairwise_alignment_file)

        if len(pairwise_alignment_dict) != 2:
            raise ValueError(
                f"""
                Expected 2 sequences in the pairwise alignment file,
                found {len(pairwise_alignment_dict)}.
                """
            )

        psa_map = {}
        q_count, s_count = 0, 0
        qseq, sseq = list(pairwise_alignment_dict.values())

        increment_dict = {
            (True, True): lambda q, s: (q + 1, s + 1),
            (True, False): lambda q, s: (q + 1, s),
            (False, True): lambda q, s: (q, s + 1),
            (False, False): lambda q, s: (q, s),
        }

        init_case = {
            True: lambda d: d + 1,
            False: lambda d: d,
        }

        for q_res, s_res in zip(qseq, sseq):

            inc_case = (s_res != "-", q_res != "-")
            q_count, s_count = increment_dict[inc_case](q_count, s_count)

            if s_res != "-" and q_res != "-":
                psa_map[init_case[q_count==0](q_count)] = (
                    init_case[s_count==0](s_count)
                )

        return psa_map

    @staticmethod
    def get_mapped_residue(
        psa_map: dict,
        codon_number: int,
        p_name: str = "",
    ) -> tuple[int | None, str]:
        """ Given a pairwise alignment map and a codon number, return the
        mapped residue number.

        ## Arguments:

        - **psa_map (dict)**:<br />
            Pairwise alignment map of codon number to residue number

        - **codon_number (int)**:<br />
            Codon number to map

        - **p_name (str, optional):**:<br />
            Protein name for warning messages. Defaults to "".

        ## Returns:

        - **tuple[int | None, str]**:<br />
            Mapped residue number and warning message (if any)
        """

        warn_msg = ""

        if len(psa_map) == 0:
            res_num_mapped = codon_number
        else:
            try:
                res_num_mapped = psa_map[codon_number]
            except KeyError:
                warn_msg += (
                    f"""
                    Residue number {codon_number} not found in
                    pairwise alignment map for protein {p_name}.
                    Skipping...
                    """
                )
                res_num_mapped = None

        return res_num_mapped, warn_msg

    @staticmethod
    def get_closest_mapped_residue(
        psa_map: dict,
        codon_number: int,
        which: str = MiscStrEnum.LOWER,
    ) -> int:
        """ Given a pairwise alignment map and a codon number, return the
        closest mapped residue number.

        ## Arguments:

        - **psa_map (dict)**:<br />
            Pairwise alignment map of codon number to residue number

        - **codon_number (int)**:<br />
            Codon number to map

        - **which (str, optional):**:<br />
            Whether to return the closest lower or higher residue number.
            Defaults to "lower".

        ## Returns:

        - **int | None**:<br />
            Closest mapped residue number.
        """

        mapped_residues = sorted(psa_map.keys())

        if which == MiscStrEnum.LOWER:
            lower_residues = [
                res for res in mapped_residues if res <= codon_number
            ]
            if lower_residues:
                closest_residue = psa_map[lower_residues[-1]]
            else:
                warnings.warn(
                    f"""
                    No lower residue found for codon number {codon_number}
                    in the pairwise alignment map.
                    Returning higher residue instead.
                    """
                )
                closest_residue = PairwiseSequenceAlignment.get_closest_mapped_residue(
                    psa_map=psa_map,
                    codon_number=codon_number,
                    which=MiscStrEnum.UPPER
                )

        elif which == MiscStrEnum.UPPER:
            higher_residues = [
                res for res in mapped_residues if res >= codon_number
            ]
            if higher_residues:
                closest_residue = psa_map[higher_residues[0]]
            else:
                warnings.warn(
                    f"""
                    No higher residue found for codon number {codon_number}
                    in the pairwise alignment map.
                    Returning lower residue instead.
                    """
                )
                closest_residue = PairwiseSequenceAlignment.get_closest_mapped_residue(
                    psa_map=psa_map,
                    codon_number=codon_number,
                    which=MiscStrEnum.LOWER
                )

        else:
            raise ValueError(
                f"which must be one of {MiscStrEnum.LOWER.value, MiscStrEnum.UPPER.value}, but got {which}"
            )

        return closest_residue

    @staticmethod
    def get_sequence_identity(
        qseq: str,
        sseq: str,
        start: int,
        end: int,
        reference: str = PSATerm.QSEQ,
        as_percentage: bool=True
    ):
        """ Calculate sequence identity between two aligned sequences in a given range.

        ## Arguments:

        - **qseq (str)**:<br />
            Aligned query sequence (with gaps)

        - **sseq (str)**:<br />
            Aligned subject sequence (with gaps)

        - **start (int)**:<br />
            Start position (1-based) of query or subject sequence (depending on reference)

        - **end (int)**:<br />
            End position (1-based) of query or subject sequence (depending on reference)

        - **reference (str, optional):**:<br />
            Reference sequence for start and end positions.
            Should be either 'qseq' or 'sseq'. Defaults to 'qseq'.

        - **as_percentage (bool, optional):**:<br />
            Whether to return sequence identity as percentage. Defaults to True.

        ## Returns:

        - **float**:<br />
            Sequence identity percentage (if as_percentage=True) or number of
            matches (if as_percentage=False)
        """

        seq_count = 0
        start_idx = end_idx = None
        ref_seq = {
            PSATerm.QSEQ: qseq,
            PSATerm.SSEQ: sseq,
        }
        for i, a in enumerate(ref_seq[reference]):
            if a != '-':
                seq_count += 1
            if seq_count == start and start_idx is None:
                start_idx = i
            if seq_count == end:
                end_idx = i
                break

        if start_idx is None or end_idx is None:
            return 0.0

        aligned_qseq = qseq[start_idx:end_idx+1]
        aligned_sseq = sseq[start_idx:end_idx+1]

        aligned_qseq = qseq[start-1:end]
        aligned_sseq = sseq[start-1:end]

        matches = sum(
            1 for a, b in zip(aligned_qseq, aligned_sseq)
            if (a == b and a != '-' and b != '-')
        )
        length = len(aligned_qseq)

        if length == 0:
            return 0.0

        if as_percentage:
            identity = (matches / length) * 100
        else:
            identity = matches

        return identity

    @staticmethod
    def get_gap(
        qseq: str,
        sseq: str,
        start: int,
        end: int,
        reference: str = PSATerm.QSEQ,
        as_percentage: bool=True
    ) -> float | int:
        """ Calculate gap percentage between two aligned sequences in a given range.

        ## Arguments:

        - **qseq (str)**:<br />
            Query sequence (aligned)

        - **sseq (str)**:<br />
            Subject sequence (aligned)

        - **start (int)**:<br />
            Start position (1-based) of query or subject sequence (depending on reference)

        - **end (int)**:<br />
            End position (1-based) of query or subject sequence (depending on reference)

        - **reference (str, optional):**:<br />
            Reference sequence for start and end positions.
            Should be either 'qseq' or 'sseq'. Defaults to 'qseq'.

        - **as_percentage (bool, optional):**:<br />
            Whether to return gap percentage. Defaults to True.

        ## Returns:

        - **float | int**:<br />
            Gap percentage (if as_percentage=True) or number of gaps (if as_percentage=False)
        """

        seq_count = 0
        start_idx = end_idx = None
        ref_seq = {
            PSATerm.QSEQ: qseq,
            PSATerm.SSEQ: sseq,
        }
        for i, a in enumerate(ref_seq[reference]):
            if a != '-':
                seq_count += 1
            if seq_count == start and start_idx is None:
                start_idx = i
            if seq_count == end:
                end_idx = i
                break

        if start_idx is None or end_idx is None:
            return 100.0

        aligned_qseq = qseq[start_idx:end_idx+1]
        aligned_sseq = sseq[start_idx:end_idx+1]

        aligned_qseq = qseq[start-1:end]
        aligned_sseq = sseq[start-1:end]

        gaps = sum(
            1 for a, b in zip(aligned_qseq, aligned_sseq)
            if (a == '-' or b == '-')
        )
        length = len(aligned_qseq)

        if length == 0:
            return 100.0

        if as_percentage:
            gap_percentage = (gaps / length) * 100
        else:
            gap_percentage = gaps

        return gap_percentage