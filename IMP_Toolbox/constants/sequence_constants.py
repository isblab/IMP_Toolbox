
from enum import StrEnum, auto

UNIPROT_ISOFORM_SEPARATOR = "-"

class MolType(StrEnum):
    PROT = auto()
    NUCL = auto()

class PSATerm(StrEnum):
    QSEQ = auto()
    SSEQ = auto()

class PSAProgram(StrEnum):
    """
    Enum for supported pairwise sequence alignment programs.
    """
    NEEDLE = auto()
    STRETCHER = auto()
    WATER = auto()

class PSAAttribute(StrEnum):
    """
    Enum for supported attributes of pairwise sequence alignment objects.
    """
    QID = auto()
    SID = auto()
    QSEQ = auto()
    SSEQ = auto()
    QALN = auto()
    SALN = auto()
    QSTART = auto()
    QEND = auto()
    SSTART = auto()
    SEND = auto()
    LENGTH = auto()
    SCORE = auto()
    NIDENTITY = auto()
    PIDENTITY = auto()
    NSIMILARITY = auto()
    PSIMILARITY = auto()
    NGAPS = auto()
    PGAPS = auto()
    MOLTYPE = auto()
    PROGRAM = auto()
    GAPOPEN = auto()
    GAPEXTEND = auto()
    MATRIX = auto()
    RAW = auto()