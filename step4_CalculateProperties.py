"""Calculate RDKit molecular descriptors for every SMILES in a CSV file.

Example
-------
python drug-like.py test.csv test_properties.csv

The input must contain an ``ID`` column and a ``Smiles`` column. Invalid or
missing SMILES are retained in the output and receive empty descriptor values.
"""

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, MolSurf
from utils.SA_Score import sascorer
from utils.calculate_qed import calculate_qed
from utils.calculate_alert import calculate_alert


def calc_properties(mol: Chem.rdchem.Mol):
    """Return a dictionary of the requested descriptors for an RDKit molecule."""
    return {
        "HBD": Lipinski.NumHDonors(mol),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "MW": round(Descriptors.MolWt(mol), 3),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RB": Lipinski.NumRotatableBonds(mol),
        "AromaticRings": Lipinski.NumAromaticRings(mol),
        "TPSA": round(MolSurf.TPSA(mol), 3),
        "MR": round(Descriptors.MolMR(mol), 3),
        "SAscore": round(sascorer.calculateScore(mol), 3),
        "Alerts": calculate_alert(mol)
    }


def add_properties(df: pd.DataFrame, smiles_column="Smiles"):
    """Calculate descriptors for each row and return an augmented DataFrame."""
    property_columns = ["HBD", "LogP", "MW", "HBA", "RB", "AromaticRings", "TPSA", "MR", "SAscore", "Alerts"]
    rule_columns = ["Lipinski_Violations", "Lipinski_Ro5", "Veber", "Ghose_Filter", "My_Rule", "QED_w"]
    records = []

    for smiles in df[smiles_column]:
        mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
        if mol is None:
            record = {column: None for column in property_columns}
        else:
            props = calc_properties(mol)
            rules = check_rules(props)
            record = {**props, **rules}
        records.append(record)

    properties_df = df.copy()
    properties_df[property_columns + rule_columns] = pd.DataFrame(records, index=properties_df.index)
    return properties_df


def check_rules(props: dict[str]):
    '''
    Parameters
    ---
        props : dict
            Dictionary or properties of an RDKit molecule.
    
    Returns
    ---
        results : dict
            Dictionary of evaluations.
    '''
    results = {}

    # Lipinski Ro5
    ro5_violations = sum([
        props["MW"] > 500,
        props["LogP"] > 5,
        props["HBD"] > 5,
        props["HBA"] > 10
    ])
    results["Lipinski_Violations"] = ro5_violations
    results["Lipinski_Ro5"] = ro5_violations <= 1  # allow 1 violation

    # Veber Rules
    results["Veber"] = (props["RB"] <= 10) and (props["TPSA"] <= 140)

    # Ghose Filter
    results["Ghose_Filter"] = (
        -0.4 <= props["LogP"] <= 5.6 and
        160 <= props["MW"] <= 480 and
        40 <= props["MR"] <= 130 and
        20 <= (props["MW"]/12) <= 70
    )

    # My Rule
    results["My_Rule"] = (calc_MyRuleScore(props))
    # Weighted QED
    # qed_weights = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95]
    qed_weights = [0.1417, 0.1108, 0.0709, 0.0639, 0.1604, 0.0757, 0.3288, 0.0479]
    results["QED_w"] = round(calculate_qed(
        values=[props[i] for i in ["MW", "LogP", "HBD", "HBA", "TPSA", "RB", "AromaticRings", "Alerts"]],
        weights=qed_weights
        ), 3)

    return results


def calc_MyRuleScore(props):
    hbd = props["HBD"]
    logp = props["LogP"]
    mw = props["MW"]
    hba = props["HBA"]
    tpsa = props["TPSA"]
    rb = props["RB"]
    arom = props["AromaticRings"]

    f_hbd = 1.0 if hbd <= 2 else 0.5 if hbd <= 5 else 0.0

    if 2 <= logp <= 5:
        f_logp = 1.0
    elif 1 <= logp < 2 or 5 < logp <= 5.5:
        f_logp = 0.5
    else:
        f_logp = 0.0

    f_mw = 1.0 if mw <= 500 else 0.5 if mw <= 600 else 0.0
    f_hba = 1.0 if hba <= 7 else 0.5 if hba <= 10 else 0.0
    f_tpsa = 1.0 if tpsa <= 120 else 0.5 if tpsa <= 140 else 0.25
    f_rb = 1.0 if rb <= 8 else 0.5 if rb <= 10 else 0.25

    if 2 <= arom <= 4:
        f_arom = 1.0
    elif arom == 1 or arom == 5:
        f_arom = 0.5
    else:  # arom == 0 or arom > 5
        f_arom = 0.25

    My_Rule_Score = round(
        0.30 * f_hbd
        + 0.25 * f_logp
        + 0.15 * f_mw
        + 0.10 * f_hba
        + 0.10 * f_tpsa
        + 0.05 * f_rb
        + 0.05 * f_arom,
        3,
    )

    return My_Rule_Score


def main():
    parser = ArgumentParser(description="Add RDKit molecular properties to a CSV file.")
    parser.add_argument("input_file", nargs="?", default="test.csv")
    parser.add_argument("--smiles-column", default="Smiles", help="Name of the SMILES column.")
    args = parser.parse_args()

    input_file = Path(args.input_file)
    if not input_file.is_file():
        raise FileNotFoundError(f"CSV file not found: {input_file}")
    output_file = input_file.parent / (input_file.stem + "_properties.csv")

    df = pd.read_csv(input_file)
    required_columns = {"ID", args.smiles_column}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required column(s): {sorted(missing_columns)}")

    properties_df = add_properties(df, args.smiles_column)
    properties_df.to_csv(output_file, index=False)
    print(f"Wrote {len(properties_df)} molecules to: {output_file}")

    invalid_count = properties_df["MW"].isna().sum()
    if invalid_count:
        print(f"Invalid or missing SMILES retained with blank properties: {invalid_count}")

    print()
    print(properties_df.head())


if __name__ == "__main__":
    main()
    # smiles = "CC(=O)Oc1ccccc1C(=O)O"
    # mol = Chem.MolFromSmiles(smiles)
    # print(sascorer.calculateScore(mol))