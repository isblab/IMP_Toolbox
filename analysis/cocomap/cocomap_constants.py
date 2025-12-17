
DOCKER_BASE_COMMAND = (
    """docker run --rm -v $path_to_mount:$path_to_mount -it andrpet/cocomaps-backend:0.0.19 python /app/coco2/begin.py $config_path"""
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
