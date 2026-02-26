from enum import StrEnum, auto
from dataclasses import dataclass
from string import Template

MAX_API_RETRIES = 3
CHIMERAX_RUN_CMD = "flatpak run edu.ucsf.rbvi.ChimeraX"

class MiscStrEnum(StrEnum):
    """
    Enum for miscellaneous string constants.
    """
    LOWER = auto()
    UPPER = auto()
    TRUE = auto()
    FALSE = auto()

class FileFormat(StrEnum):
    """
    Enum for supported file formats.
    """
    FASTA = auto()
    FA = auto()
    PDB = auto()
    CIF = auto()
    MRC = auto()
    TXT = auto()
    HTML = auto()

@dataclass
class APIurl:

    uniprot_rest_api = "https://rest.uniprot.org/"
    uniprot_rest_api_sequence = Template(
        f"{uniprot_rest_api}" +
        "/uniprotkb/accessions?accessions=${ids}&format=fasta"
    )

    pdbe_graph_api = "https://www.ebi.ac.uk/pdbe/graph-api"
    pdbe_api_best_structures = Template(
        f"{pdbe_graph_api}" +
        "/mappings/best_structures/${uniprot_id}"
    )

    emdb_ftp = "https://ftp.ebi.ac.uk/pub/databases/emdb"
    emdb_ftp_map = Template(
        f"{emdb_ftp}/structures" +
        "/${emdb_id_hyphen}/map/${emdb_id_underscore}.map.gz"
    )
    emdb_ftp_mask = Template(
        f"{emdb_ftp}/structures" +
        "/${emdb_id_hyphen}/masks/${mask_name}.map"
    )

class CorrelationMetric(StrEnum):
    """
    Enum for supported correlation metrics for fitting GMMs to density maps.
    """
    CAM = auto()
    CORRELATION = auto()
    OVERLAP = auto()

@dataclass
class ChimeraXDefaults:

    correlation_metric: str | CorrelationMetric = CorrelationMetric.OVERLAP
    envelop_val: str | MiscStrEnum = MiscStrEnum.TRUE
    shift_val: str | MiscStrEnum = MiscStrEnum.TRUE
    rotate_val: str | MiscStrEnum = MiscStrEnum.TRUE
    zeros_val: str | MiscStrEnum = MiscStrEnum.FALSE
    fitmap_max_steps: int = 2000

class ChimeraXVolumeLevelType(StrEnum):
    """
    Enum for supported volume thresholding level types in ChimeraX.
    """
    LEVEL = "level"
    SD_LEVEL = "sdLevel"
    RMS_LEVEL = "rmsLevel"

@dataclass
class ChimeraXCommands:
    fitmap_cmd: str = Template(
        "fitmap #${model_count} inmap #1 " +
        "metric ${metric} " +
        "shift ${shift_val} " +
        "rotate ${rotate_val} " +
        "envelope ${envelop_val} " +
        "maxSteps ${fitmap_max_steps} " +
        "zeros ${zeros_val}"
    )
    measure_correlation_cmd: str = Template(
        "measure correlation #${model_count} in_map #1 " +
        "envelope ${envelop_val}"
    )
    open_cmd = Template("open ${file_path}")
    volume_threshold_cmd = Template(
        "volume #${model_count} ${level_type} ${level_val}"
    )

class ChimeraXCommand(StrEnum):
    """
    Enum for supported ChimeraX command names.
    """
    FITMAP = "fitmap"
    MEASURE_CORRELATION = "measure correlation"

@dataclass
class ChimeraXLogPatterns:

    correlation = {
        ChimeraXCommand.FITMAP: {
            "look_for": "Fit map",
            "regex1": r"Fit map (.+) in map (.+) using (\d+) points",
            "regex2": r"  correlation\s*=\s*([0-9]+\.[0-9]+),\s*correlation\s+about\s+mean\s*=\s*([0-9]+\.[0-9]+),\s*overlap\s*=\s*([0-9]+\.[0-9]+)"
        },
        ChimeraXCommand.MEASURE_CORRELATION: {
            "look_for": "Correlation of",
            "regex1": r"Correlation\s+of\s+(\S+)\s+#(\d+)\s+above\s+level\s+([0-9.]+e[+-]?[0-9]+)\s+in\s+(\S+)\s+#(\d+)",
            "regex2": r"correlation\s*=\s*([0-9]+\.[0-9]+),\s*correlation\s+about\s+mean\s*=\s*([0-9]+\.[0-9]+)"
        }
    }

@dataclass
class GMMParams:
    """
    Dataclass for GMM parameters.
    """
    threshold: float = 0.95
    correlation_metric : str | CorrelationMetric = CorrelationMetric.CAM
    gmm_mrc_name = Template("${f_name}_${n_gaussian}")
    log_file_name = Template("${f_name}_gmm_selection_log_chimerax")