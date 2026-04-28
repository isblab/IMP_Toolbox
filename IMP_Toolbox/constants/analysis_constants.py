import numpy as np

PAIR_SEP = "|"
RES_RANGE_SEP = "-"
MOL_COPY_SEP = "_"
MOL_RANGE_SEP = ":"
REGEX_MOLNAME = r'^([A-Za-z0-9]+)(?:_(\d+))?(?::([\d-]+))?$'
EXPECTED_MOLNAME_FORMATS = [
    "MOL_COPYIDX:RESSTART-RESEND",
    "MOL_COPYIDX",
    "MOL",
    "MOL:RESSTART-RESEND",
]
F_DTYPES = {16: np.float16, 32: np.float32, 64: np.float64}
I_DTYPES = {8: np.int8, 16: np.int16, 32: np.int32, 64: np.int64}