import os
import pandas as pd
import json
import argparse
from utils import make_protein
# upload the output json file to the server
# https://ibs.renlab.org/#/server
# export svg and edit in inkscape or affinity designer

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument(
        "-i",
        "--input_csv",
        type=str,
        default="input/protein_domains.xlsx",
    )
    args.add_argument(
        "-s",
        "--sheet_name",
        type=str,
        default="Sheet1",
    )
    args.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default="output/protein_domains",
    )

    args = args.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Read the sheet from the input xlsx file
    df = pd.read_excel(args.input_csv, sheet_name=args.sheet_name)
    df = df.ffill()
    print(df)
    proteins_dicts = {}
    for idx, row in df.iterrows():
        protein_name = row["Protein"]
        if protein_name not in proteins_dicts:
            proteins_dicts[protein_name] = make_protein(protein_name,row["Start"], row["End"])
        child_dict = {}
        child_dict["type"] = "domain"
        child_dict["option"] = {
            "id": f"Domain_{idx}_" + str(row["Domain"]),
            "position": {
                "isLocked": False,
                "start": {
                    "site": row["Res_range"].split("-")[0],
                    "display": "hide"
                },
                "end": {
                    "site": row["Res_range"].split("-")[1],
                    "display": "hide"
                }
            },
            "displayName": f"Domain_{idx}_" + str(row["Domain"]),
            "text": "",
            "nameTextStyle": {
                "fontFamily": "Arial",
                "fontStyle": "normal",
                "fontSize": 10,
                "color": "#FFFFFF",
                "rotate": 0,
                "location": "center",
                "fontWeight": "normal"
            },
            "style": {
                "color": "#0000ff",
                "texture": {
                    "type": "none",
                    "color": "#4389b4"
                },
                "gradient": "none"
            },
            "borderStyle": {
                "color": "#000000",
                "size": 1,
                "isDash": False
            },
            "name": str(row["Domain"]),
        }
        proteins_dicts[protein_name]["data"][0]["children"].append(child_dict)

    # Write the JSON to a file
    for protein_name, protein_dict in proteins_dicts.items():
        with open(os.path.join(args.output_dir, f"{protein_name}.json"), "w") as f:
            json.dump(protein_dict, f, indent=4)