import os
import pandas as pd
from pathlib import Path
import zipfile
from IMP_Toolbox.utils.file_helpers import download
from IMP_Toolbox.constants.ptm_constants import API_URLS

def get_ptmd_data(
    savedir: str,
    key: str = "ptmd_total",
    overwrite: bool = False,
):

    api_url = API_URLS.get(key)
    if not api_url:
        raise ValueError(f"Invalid key: {key}")

    savepath = Path(savedir) / f"{key}.zip"
    if savepath.exists() and not overwrite:
        print(f"PTMD data already exists at {savepath}. Skipping download.")
        return
    os.makedirs(savedir, exist_ok=True)

    print("Downloading PTMD data ...")
    download(
        url=api_url,
        filename=savepath,
        chunksize=10240000,
    )

    assert savepath.exists(), "Failed to download PTMD data"

    with zipfile.ZipFile(savepath, "r") as zip_ref:
        zip_ref.extractall(Path(savedir) / key)

    print("PTMD data downloaded and extracted successfully")


if __name__ == "__main__":

    get_ptmd_data(
        savedir="/data/omg/Projects/IMP_Toolbox/ptms",
        key="ptmd_total",
        overwrite=False,
    )

    ptmd_txt = "/data/omg/Projects/IMP_Toolbox/ptms/ptmd_total/Total.txt"
    df = pd.read_csv(ptmd_txt, sep="\t")
    print(df.head())