from pathlib import Path
from rdkit import Chem
from rdkit.Chem.QED import StructuralAlertSmarts
unwanted_path = Path("utils/unwanted-groups.txt")
if unwanted_path.exists():
    with open(unwanted_path, "r") as fin:
        UNWANTED_GROUPS = [Chem.MolFromSmarts(line.strip()) for line in fin.readlines()]
else:
    print(unwanted_path, "not found, using built-in alert structures.")
    UNWANTED_GROUPS = [Chem.MolFromSmarts(smarts) for smarts in StructuralAlertSmarts]


def calculate_alert(mol: Chem.rdchem.Mol, unwanted_groups: list[Chem.rdchem.Mol]=UNWANTED_GROUPS):
    alert = 0
    for group in unwanted_groups:
        if mol.HasSubstructMatch(group):
            alert += 1

    return alert


test_mol = Chem.MolFromSmiles("CN/C(=N\CCSCc1nc[nH]c1C)NC#N") # CIMETIDINE

if __name__ == "__main__":
    print("You can modify unwanted groups at", unwanted_path.resolve())
    print()
    print("Test mol: CIMETIDINE CN/C(=N\CCSCc1nc[nH]c1C)NC#N")
    print("Alert:   ", calculate_alert(test_mol))

