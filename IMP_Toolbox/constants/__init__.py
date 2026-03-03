"""
constants
===========

- Constants used across the IMP Toolbox.

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
    class APIurl {
        <<dataclass>>
        + str uniprot_rest_api
        + Any uniprot_rest_api_sequence
        + str pdbe_graph_api
        + Any pdbe_api_best_structures
        + str emdb_ftp
        + Any emdb_ftp_map
        + Any emdb_ftp_mask
    }

    class ChimeraXCommand {
        <<enumeration>>
        + Any FITMAP
        + str MEASURE_CORRELATION
        + Any OPEN
        + Any CLOSE
        + Any VOLUME
        + Any VOLUME_THRESHOLD
        + str VOLUME_ADD
        + Any RENAME
    }

    class ChimeraXVolumeLevelType {
        <<enumeration>>
        + str LEVEL
        + str SD_LEVEL
        + str RMS_LEVEL
    }

    class FileFormat {
        <<enumeration>>
        + Any FASTA
        + Any FA
        + Any PDB
        + Any CIF
        + Any MRC
        + Any TXT
        + Any HTML
    }

```


```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram

    class ChimeraXDefaults {
        <<dataclass>>
        + str | CorrelationMetric correlation_metric
        + str | MiscStrEnum envelop_val
        + str | MiscStrEnum shift_val
        + str | MiscStrEnum rotate_val
        + str | MiscStrEnum zeros_val
        + int fitmap_max_steps
    }

    class CorrelationMetric {
        <<enumeration>>
        + Any CAM
        + Any CORRELATION
        + Any OVERLAP
    }

    class MiscStrEnum {
        <<enumeration>>
        + Any LOWER
        + Any UPPER
        + Any TRUE
        + Any FALSE
        + Any NUM_PTS
        + Any MODEL_MRC
        + Any REF_MRC
        + str SAMPLE_A
        + str SAMPLE_B
        + Any NUM_GAUSSIANS
    }

    class GMMParams {
        <<dataclass>>
        + float threshold
        + str | CorrelationMetric correlation_metric
        + Any gmm_mrc_name
        + Any log_file_name
    }

    class MRCComparisonMode {
        <<enumeration>>
        + Any PASANI
        + Any CHIMERAX
    }

    ChimeraXDefaults *-- CorrelationMetric

    ChimeraXDefaults *-- MiscStrEnum

    GMMParams *-- CorrelationMetric
```

```mermaid
---
config:
    class:
        hideEmptyMembersBox: true
---
classDiagram
    class PSAAttribute {
        <<enumeration>>
        + Any QID
        + Any SID
        + Any QSEQ
        + Any SSEQ
        + Any QALN
        + Any SALN
        + Any QSTART
        + Any QEND
        + Any SSTART
        + Any SEND
        + Any LENGTH
        + Any SCORE
        + Any NIDENTITY
        + Any PIDENTITY
        + Any NSIMILARITY
        + Any PSIMILARITY
        + Any NGAPS
        + Any PGAPS
        + Any MOLTYPE
        + Any PROGRAM
        + Any GAPOPEN
        + Any GAPEXTEND
        + Any MATRIX
        + Any RAW
    }

    class PSAProgram {
        <<enumeration>>
        + Any NEEDLE
        + Any STRETCHER
        + Any WATER
    }

    class PSATerm {
        <<enumeration>>
        + Any QSEQ
        + Any SSEQ
    }

    class BestStructureCol {
        <<enumeration>>
        + Any UNIPROT_ID
        + Any PROTEIN
        + Any CHAIN_ID
        + Any PDB_ID
        + Any COVERAGE
        + Any UNP_START
        + Any UNP_END
        + Any START
        + Any END
        + Any RESOLUTION
    }
```


"""