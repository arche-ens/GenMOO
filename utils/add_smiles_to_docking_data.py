#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="The docking result file from Schrödinger may lack the SMILES column. This script helps to add SMILES metadata to the docking result file based on ID matching.")
    parser.add_argument("-i", "--input", type=Path, required=True, help="Input CSV file containing the title column, ask for an additional metadata.")
    parser.add_argument("-v", "--values", type=Path, required=True, help="CSV file containing ID and SMILES columns, provide the metadata.")
    return parser.parse_args()


def force_read_csv(path: Path) -> pd.DataFrame:
    for encoding in ("utf-8", "latin-1"):
        try:
            return pd.read_csv(path, encoding=encoding, low_memory=False)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError(f"Could not decode {path} as UTF-8 or latin-1")


def merge_values(input_file: Path, values_file: Path):
    output_file = input_file
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    if not values_file.exists():
        raise FileNotFoundError(f"Values file not found: {values_file}")

    # Read CSV files (Glide/Maestro output may be latin-1 encoded)
    df = force_read_csv(input_file)
    values = force_read_csv(values_file)

    # Check required columns
    if "s_m_title" not in df.columns:
        raise ValueError("Input CSV must contain a 's_m_title' column")

    if not {"ID", "Smiles"}.issubset(values.columns):
        raise ValueError("Values CSV must contain 'ID' and 'Smiles' columns")

    # Merge values
    result = df.merge(
        values[["ID", "Smiles"]],
        how="left",
        left_on="s_m_title",
        right_on="ID"
    )

    result = result.drop(columns=["ID"])

    result.to_csv(output_file, index=False)
    print(f"Saved merged file to: {output_file}")


def main():
    args = parse_args()
    merge_values(input_file=args.input, values_file=args.values)


if __name__ == "__main__":
    main()