
DOCKER_BASE_COMMAND = (
    """sudo docker run --rm -v $path_to_mount:$path_to_mount -it andrpet/cocomaps-backend:0.0.19 python /app/coco2/begin.py $config_path"""
)

TEMPLATE_CONFIG = {
    "HBOND_DIST": 3.9,
    "HBOND_ANGLE": 90,
    "SBRIDGE_DIST": 4.5,
    "WBRIDGE_DIST": 3.9,
    "CH_ON_DIST": 3.6,
    "CH_ON_ANGLE": 110,
    "CUT_OFF": 5,
    "APOLAR_TOLERANCE": 0.5,
    "POLAR_TOLERANCE": 0.5,
    "PI_PI_DIST": 5.5,
    "PI_PI_THETA": 80,
    "PI_PI_GAMMA": 90,
    "ANION_PI_DIST": 5,
    "LONEPAIR_PI_DIST": 5,
    "AMINO_PI_DIST": 5,
    "CATION_PI_DIST": 5,
    "METAL_DIST": 3.2,
    "HALOGEN_THETA1": 165,
    "HALOGEN_THETA2": 120,
    "C_H_PI_DIST": 5.0,
    "C_H_PI_THETA1": 120,
    "C_H_PI_THETA2": 30,
    "NSOH_PI_DIST": 4.5,
    "NSOH_PI_THETA1": 120,
    "NSOH_PI_THETA2": 30
}

ATOM_COL_NAMES = {
    "H-bond": ["Atom 1", "Atom 2"],
    "Salt_bridge": ["Atom 1", "Atom 2"],
    "SS_bond": ["Atom 1", "Atom 2"],
    "C-H_ON": ["Atom 1", "Atom 2"],
    "Polar_vdw": ["Atom 1", "Atom 2"],
    "Apolar_vdw": ["Atom 1", "Atom 2"],
    "Halogen_bond": ["Atom 1", "Atom 2"],
    "Proximal": ["Atom 1", "Atom 2"],
    # "Metal_Mediated": ["Atom 1", "Atom 2", "Metal Identity"], # Not implemented
    # "Water_Mediated": ["Atom 1", "Atom 2", "Water Identity"], # Not implemented
    "Lone_pair_pi": ["Lone_pair Atom","Ring From", "Lone_pair From"],
    "Cation_pi": ["Cation Atom", "Cation From", "Ring From"],
    "Anion_pi": ["Anion Atom", "Anion From", "Ring From"],
    "Amino_pi": ["Polar Atom", "Polar From", "Ring From"],
    "N-S-O-H_pi": ["C Atom", "C Atom From", "Ring From"],
    "C-H_pi": ["C Atom", "C Atom From", "Ring From"],
}

VALID_INTERACTIONS = {
    # "Clash": "Clash", # there is a bug here
    "H-bond": "H-bond",
    "Salt_bridge": "Salt-bridge",
    "SS_bond": "S-S Bond",
    "C-H_ON": "CH-O/N bond",
    "Polar_vdw": "Polar vdW contact",
    "Apolar_vdw": "Apolar vdW contact",
    "Halogen_bond": "Halogen bond",
    "Proximal": "Proximal contact",
    # "Metal_Mediated": "Metal mediated contact", # not implemented
    # "Water_Mediated": "Water mediated contact", # not implemented
    "Lone_pair_pi": "lp-π interaction",
    "Cation_pi": "Cation-π interaction",
    "Anion_pi": "Anion-π interaction",
    "Amino_pi": "Amino-π interaction",
    "N-S-O-H_pi": "O/N/SH-π interaction",
    "C-H_pi": "CH-π interaction",
    "pi-pi": "π-π interaction",
}

REP_ATOMS = {
    "GLY": "CA",
    "ALA": "CB",
    "VAL": "CB",
    "LEU": "CB",
    "ILE": "CB",
    "THR": "CB",
    "SER": "CB",
    "MET": "CB",
    "CYS": "CB",
    "PRO": "CB",
    "PHE": "CB",
    "TYR": "CB",
    "TRP": "CB",
    "HIS": "CB",
    "LYS": "CB",
    "ARG": "CB",
    "ASP": "CB",
    "GLU": "CB",
    "ASN": "CB",
    "GLN": "CB",
}