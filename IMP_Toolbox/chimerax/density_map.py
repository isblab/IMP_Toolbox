import re
from IMP_Toolbox.constants.imp_toolbox_constants import (
    ChimeraXCommand,
    ChimeraXLogPatterns,
)

def parse_chimerax_correlation_log(
    chimerax_log: str,
    command: str | ChimeraXCommand = ChimeraXCommand.FITMAP,
) -> list:
    """ Parse the log file from ChimeraX output.

    ## Arguments:

    - **chimerax_log (str)**:<br />
        Path to the ChimeraX log file.

    - **command (str, optional):**:<br />
        The command for which to parse the log.
        Options are "fitmap" and "measure correlation". Default is "fitmap".

    ## Returns:

    - **list**:<br />
        A list of dictionaries containing the parsed correlation results.
    """

    attrs = ChimeraXLogPatterns.correlation

    assert command in attrs, (
        f"Command {command} not recognized. Valid options are: {list(attrs.keys())}"
    )

    log_lines = []
    with open(chimerax_log, 'r') as f:
        log_lines = f.readlines()

    correlation_list = []

    for idx, line in enumerate(log_lines):
        if (
            idx + 1 >= len(log_lines) or
            line.startswith(attrs[command]["look_for"]) is False
        ):
            continue

        regex1 = attrs[command]["regex1"]
        regex2 = attrs[command]["regex2"]

        match1 = re.match(regex1, line)
        match2 = re.match(regex2, log_lines[idx + 1].strip())

        if command == ChimeraXCommand.FITMAP:
            model_map, ref_map, num_points = match1.groups()
            correlation, cam, overlap = match2.groups()

        elif command == ChimeraXCommand.MEASURE_CORRELATION:
            model_map, ref_map = match1.groups()
            correlation, cam = match2.groups()

        correlation_list.append({
            "model_mrc": model_map.replace("mrc", "").strip(),
            "ref_mrc": ref_map.replace("mrc", "").strip(),
            "correlation": correlation,
            "cam": cam,
        })

    return correlation_list