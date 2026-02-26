import mrcfile
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from IMP_Toolbox.utils.api_helpers import request_session
from IMP_Toolbox.constants.imp_toolbox_constants import (
    APIurl,
    MAX_API_RETRIES,
)

def fetch_emdb_map(emdb_id: str, max_retries: int = MAX_API_RETRIES) -> bytes:
    """ Fetch density map from EMDB

    ## Arguments:

    - **emdb_id (str)**:<br />
        EMDB ID for which the density map is to be fetched.
        This should be in the format "EMD-XXXX" where XXXX is a 4 digit number.

    - **max_retries (int, optional):**:<br />
        Maximum number of retries for the API request. Default is 3.

    ## Returns:

    - **bytes**:<br />
        Density map file content in bytes.
    """

    EMDB_MAP_URL = APIurl.emdb_ftp_map.substitute(
        emdb_id_hyphen=emdb_id,
        emdb_id_underscore=emdb_id.lower().replace("-", "_")
    )

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(EMDB_MAP_URL)

    if response.status_code == 200:
        print("Successfully fetched EMDB map for given EMDB id")
        return response.content

    else:
        raise Exception("Error while requesting EMDB map for given EMDB id")

def fetch_emdb_mask(
    emdb_id: str,
    mask_name: str,
    max_retries: int = MAX_API_RETRIES,
) -> bytes:
    """ Fetch mask from EMDB

    ## Arguments:

    - **emdb_id (str)**:<br />
        EMDB ID for which the mask is to be fetched.
        This should be in the format "EMD-XXXX" where XXXX is a 4 digit number.

    - **mask_name (str)**:<br />
        Name of the mask to be fetched. This should be in the format "mask_XX".

    - **max_retries (int, optional):**:<br />
        Maximum number of retries for the API request. Default is 3.

    ## Returns:

    - **bytes**:<br />
        Mask file content in bytes.
    """

    if mask_name is None:
        raise ValueError("Mask name must be provided to fetch EMDB mask")

    EMDB_MASK_URL = APIurl.emdb_ftp_mask.substitute(
        emdb_id_hyphen=emdb_id,
        mask_name=mask_name
    )

    req_sess = request_session(max_retries=max_retries)
    response = req_sess.get(EMDB_MASK_URL)

    if response.status_code == 200:
        print("Successfully fetched EMDB mask for given EMDB id and mask name")
        return response.content

    else:
        raise Exception(
            "Error while requesting EMDB mask for given EMDB id and mask name"
        )

def extract_voxel_data(mrc_files: list) -> tuple:
    """ Extract grid points and value points from mrc files.

    ## Arguments:

    - **mrc_files (list)**:<br />
        List of paths to mrc files.

    ## Returns:

    - **tuple**:<br />
        A tuple containing two lists: grid_points, grid_values and voxel sizes.
    """

    grid_points = []
    grid_values = []
    voxel_sizes = []

    for _idx, mrc_path in enumerate(mrc_files):
        with mrcfile.open(mrc_path, 'r', permissive=True) as mrc:
            origin = mrc.header.origin
            voxel_size = mrc.voxel_size
            xvals = voxel_size.x * np.arange(mrc.data.shape[2]) + origin.x
            yvals = voxel_size.y * np.arange(mrc.data.shape[1]) + origin.y
            zvals = voxel_size.z * np.arange(mrc.data.shape[0]) + origin.z
            grid_points.append((xvals, yvals, zvals))
            grid_values.append(mrc.data.transpose(2, 1, 0).copy())
            voxel_sizes.append((voxel_size.x, voxel_size.y, voxel_size.z))

    return grid_points, grid_values, voxel_sizes

def interpolate_from_voxel_data(voxel_data: tuple) -> dict:
    """ Get interpolators from voxel data.

    ## Arguments:

    - **voxel_data (tuple)**:<br />
        A tuple containing three lists: grid_points, grid_values and voxel_sizes.
        Each of which is a list of length equal to the number of mrc files.

    ## Returns:

    - **dict**:<br />
        A dictionary containing the interpolators for each mrc file in the voxel data.
    """

    interpolators = []
    grid_points, grid_values, _ = voxel_data

    for (xvals, yvals, zvals), values in zip(grid_points, grid_values):

        interpolator = RegularGridInterpolator(
            points=(xvals, yvals, zvals),
            values=values,
            bounds_error=False,
            fill_value=0.0
        )
        interpolators.append(interpolator)

    return interpolators

def get_meshgrid_from_voxel_data(
    voxel_data: tuple,
    return_type: str = "coordinates",
    voxel_size: float | None = None,
) -> tuple | np.ndarray:
    """ Get meshgrid from voxel data.

    ## Arguments:

    - **voxel_data (tuple)**:<br />
        Voxel data tuple containing grid points, grid values and voxel sizes.

    - **return_type (str, optional):**:<br />
        Type of the return value. Can be "arrays" or "coordinates".

    - **voxel_size (float, optional):**:<br />
        Voxel size to be used for generating the meshgrid if the voxel size from
        the voxel data is not to be used. If None, the voxel size from the voxel
        data will be used. Default is None.

    ## Returns:

    - **tuple | np.ndarray**:<br />
        - A tuple containing three 3D numpy arrays representing the x, y, and z
        coordinates of the meshgrid. OR
        - A 2D numpy array of shape (N, 3) containing the coordinates of the
        meshgrid points and a tuple representing the shape of the meshgrid.
    """

    grid_points, _, voxel_sizes = voxel_data

    if voxel_size is None:
        voxel_size = voxel_sizes[0][0]

    grid_func = lambda grid_points, axis_idx: np.hstack(
        [data[axis_idx] for data in grid_points]
    )

    grid_pts_stack_x = grid_func(grid_points, 0)
    grid_pts_stack_y = grid_func(grid_points, 1)
    grid_pts_stack_z = grid_func(grid_points, 2)

    x, y, z = np.meshgrid(
        np.arange(
            np.min(grid_pts_stack_x),
            np.max(grid_pts_stack_x) + voxel_size,
            voxel_size
        ),
        np.arange(
            np.min(grid_pts_stack_y),
            np.max(grid_pts_stack_y) + voxel_size,
            voxel_size
        ),
        np.arange(
            np.min(grid_pts_stack_z),
            np.max(grid_pts_stack_z) + voxel_size,
            voxel_size
        ),
        indexing='ij'
    )

    if return_type == "arrays":
        return x, y, z

    elif return_type == "coordinates":
        return (
            np.array([
                (i, j, k)
                for i, j, k in zip(x.flatten(), y.flatten(), z.flatten())
            ]),
            x.shape
        )

    else:
        raise ValueError("Invalid return type. Must be 'arrays' or 'coordinates'.")

def get_interpolated_values_from_voxel_data(
    interpolators: list,
    voxel_data: tuple,
    voxel_size: float | None = None,
) -> list:
    """ Get interpolated values from voxel data at desired points.

    ## Arguments:

    - **voxel_data (tuple)**:<br />
        Voxel data tuple containing grid points, grid values and voxel sizes.

    - **desired_points (np.ndarray)**:<br />
        A 2D numpy array of shape (N, 3) containing the desired points
        at which to get the interpolated values.

    - **voxel_size (float, optional):**:<br />
        Voxel size to be used for generating the desired points if the desired
        points are not provided. If None, the voxel size from the voxel data
        will be used. Default is None.

    ## Returns:
    - **list**:<br />
        A list of 3D numpy arrays containing the interpolated values at the desired points
        for each interpolator.
    """

    if voxel_size is None:
        _, _, voxel_sizes = voxel_data
        voxel_size = voxel_sizes[0][0]

    desired_points, _shape = get_meshgrid_from_voxel_data(
        voxel_data=voxel_data,
        return_type="coordinates",
        voxel_size=voxel_size,
    )

    return [i(desired_points).reshape(_shape) for i in interpolators]

def measure_correlation(
    u: np.ndarray,
    v: np.ndarray,
    zeros: bool = False,
):
    """ Compute correlation metric between two vectors.

    Similar to the implementation in ChimeraX's fitmap command.
    See `overlap_and_correlation` function in ChimeraX source code.

    ## Arguments:

    - **u (np.ndarray)**:<br />
        vector 1.

    - **v (np.ndarray)**:<br />
        vector 2.

    - **zeros (bool, optional):**:<br />
        Whether to include zero values in the correlation calculation.
        Default is False.

    ## Returns:

    - **tuple**:<br />
        - overlap: sum of element-wise product of u and v.
        - correlation: correlation between u and v.
    """

    u, v = u.flatten().astype(np.float64), v.flatten().astype(np.float64)

    pts = len(u)

    if zeros is False:
        mask = (u != 0) & (v != 0)
        u, v = u[mask], v[mask]

    overlap = np.sum(u * v, dtype=np.float64)
    m1 = np.sqrt(np.sum(u ** 2, dtype=np.float64))
    m2 = np.sqrt(np.sum(v ** 2, dtype=np.float64))
    corr = overlap / (m1 * m2) if m1 * m2 != 0 else 0

    return overlap, corr, pts

def get_correlation_metrics(
    voxel_data1: tuple,
    voxel_data2: tuple,
    voxel_size: float | None = None,
    zeros: bool = False,
) -> tuple:
    """ Get correlation metrics between two sets of voxel data.

    ## Arguments:

    - **voxel_data1 (tuple)**:<br />
        A tuple containing three lists: grid_points, grid_values and voxel_sizes
        for the first set of voxel data.

    - **voxel_data2 (tuple)**:<br />
        A tuple containing three lists: grid_points, grid_values and voxel_sizes
        for the second set of voxel data.

    - **zeros (bool, optional):**:<br />
        Whether to include zero values in the correlation calculation.
        Default is False.

    ## Returns:

    - **tuple**:<br />
        A tuple containing three values: overlap, correlation and correlation over mean.
    """

    voxel_data = (
        list(voxel_data1[0]) + list(voxel_data2[0]),
        list(voxel_data1[1]) + list(voxel_data2[1]),
        list(voxel_data1[2]) + list(voxel_data2[2]),
    )

    interpolators1 = interpolate_from_voxel_data(voxel_data=voxel_data1)
    interpolators2 = interpolate_from_voxel_data(voxel_data=voxel_data2)

    values1 = get_interpolated_values_from_voxel_data(
        interpolators=interpolators1,
        voxel_data=voxel_data,
        voxel_size=voxel_size,
    )
    values2 = get_interpolated_values_from_voxel_data(
        interpolators=interpolators2,
        voxel_data=voxel_data,
        voxel_size=voxel_size,
    )

    u = np.sum(values1, axis=0, dtype=np.float64)
    v = np.sum(values2, axis=0, dtype=np.float64)

    overlap, corr, pts = measure_correlation(
        u=u,
        v=v,
        zeros=zeros,
    )

    _, corr_over_mean, _ = measure_correlation(
        u=u - np.mean(u),
        v=v - np.mean(v),
        zeros=zeros,
    )

    return overlap, corr, corr_over_mean, pts