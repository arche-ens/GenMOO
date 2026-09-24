import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, MolStandardize
cleaner = MolStandardize.rdMolStandardize.MetalDisconnector()
fragment_chooser = MolStandardize.rdMolStandardize.LargestFragmentChooser()
normalizer = MolStandardize.rdMolStandardize.Normalizer()


def clean_molecule(mol):
    mol = cleaner.Disconnect(mol) # Disconnect metals
    mol = fragment_chooser.choose(mol) # Desalt
    mol = normalizer.normalize(mol)     # Standardize nitro, charge representation
    return mol


def deduplicate_by_inchikey(cleaned_smiles):
    """Deduplicate cleaned SMILES by full InChIKey, keeping first occurrence.


    Returns (unique_smiles, num_missing_key, num_duplicates) where the SMILES
    order follows the input order. Molecules whose InChIKey cannot be computed
    are reported and skipped.
    """
    seen = {}          # InChIKey -> SMILES kept for that molecule
    num_missing_key = 0
    num_duplicates = 0


    for smiles in cleaned_smiles:
        mol = Chem.MolFromSmiles(smiles)
        inchikey_connect = Chem.MolToInchiKey(mol).split('-')[0] if mol is not None else ""
        if not inchikey_connect:
            num_missing_key += 1
            print("[WARNING] InChIKey can not be computed!!", smiles)
            continue
        if inchikey_connect in seen:
            num_duplicates += 1
            continue
        seen[inchikey_connect] = smiles

    return list(seen.values()), num_missing_key, num_duplicates


def clean_data(input_data):
    '''
    Return
    ---
        unique_smiles : list[str]
            clean and deduplicate input data, and output a list containing unique SMILES.
    '''
    with open(input_data, "r", encoding="utf-8") as fin:
        data_smiles = fin.readlines()
        print(f"{len(data_smiles):,} molecules loaded.")

        # Clean
        cleaned_smiles = []
        for smiles in data_smiles:
            try:
                mol = Chem.MolFromSmiles(smiles.strip())
                mol = clean_molecule(mol)
                cleaned_smiles.append(Chem.MolToSmiles(mol))
            except:
                print("[WARNING] Molecule can not be cleaned!!", smiles)
        print(f"{len(cleaned_smiles):,} molecules are cleaned.")
        # breakpoint()

    # Deduplicate
    unique_smiles, num_missing_key, num_duplicates = deduplicate_by_inchikey(cleaned_smiles)
    print(f"{len(unique_smiles):,} molecules are unique.")
    # breakpoint()

    return unique_smiles

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input", default="test_data.txt")
    args = parser.parse_args()

    input_data = Path(args.input)
    if not input_data.exists():
        raise FileNotFoundError(f"Input file not found: {input_data}")
    output_data = input_data.parent / (input_data.stem + "_cleaned.txt")

    unique_smiles = clean_data(input_data)
    with open(output_data, "w") as fout:
        for smiles in unique_smiles:
            fout.write(smiles.strip() + "\n")

    print(f"{input_data} --> {output_data}")