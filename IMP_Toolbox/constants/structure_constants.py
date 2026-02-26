from enum import StrEnum, auto

class BestStructureCol(StrEnum):
    """
    Enum for supported columns in best structures dataframe.
    """
    UNIPROT_ID = auto()
    PROTEIN = auto()
    CHAIN_ID = auto()
    PDB_ID = auto()
    COVERAGE = auto()
    UNP_START = auto()
    UNP_END = auto()
    START = auto()
    END = auto()
    RESOLUTION = auto()