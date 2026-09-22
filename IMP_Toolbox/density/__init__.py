"""
density
===========

- Assisting with the density map related tasks such as-
  - fetching density maps from EMDB
  - comparing density maps

<hr>

### Examples:

<hr>

#### 1. Fetch EMDB data for a given EMDB id

```python
from IMP_Toolbox.density.emdb import fetch_emdb_map, fetch_emdb_mask

emdb_id = "EMD-1703"
savepath = os.path.join("./output", f"{emdb_id}.map.gz")
extracted_savepath = os.path.join("./output", f"{emdb_id}.map")
emdb_map = fetch_emdb_map(emdb_id=emdb_id, max_retries=3)
with open(savepath, "wb") as f:
    f.write(emdb_map)
with gzip.open(savepath, "rb") as f_in:
    with open(extracted_savepath, "wb") as f_out:
        f_out.write(f_in.read())

print(f"Extracted EMDB map saved in {os.path.abspath(extracted_savepath)}")

mask = "emd_1703_msk_1"
mask_savepath = os.path.join("./output", f"{mask}.map")
emdb_mask = fetch_emdb_mask(emdb_id=emdb_id, mask_name=mask_name)
with open(mask_savepath, "wb") as f:
    f.write(emdb_mask)
print(f"EMDB mask saved in {os.path.abspath(mask_savepath)}")
```

#### 2. Compare two density maps and get correlation metrics

```python
from IMP_Toolbox.density.compare import (
    extract_voxel_data,
    get_correlation_metrics,
)

mrc_file1 = "./path/to/map1.mrc"
mrc_file2 = "./path/to/map2.mrc"

voxel_data1 = extract_voxel_data(mrc_files=[mrc_file1])
voxel_data2 = extract_voxel_data(mrc_files=[mrc_file2])
overlap, corr, corr_over_mean, pts = get_correlation_metrics(
    voxel_data1=voxel_data1,
    voxel_data2=voxel_data2,
)
print(f'''
    Overlap: {overlap},
    Correlation: {corr},
    Correlation over mean: {corr_over_mean},
    Number of points: {pts}'''
)
```
"""