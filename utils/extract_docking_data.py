#!/usr/bin/env python3
"""Extract docking scores from docking outputs.

For each class (active / inactive) the three conformation folders
``<class>_state_1`` .. ``<class>_state_3`` provide the scores
``score_1`` .. ``score_3``.  A molecule (identified by its title) may have
several rows in one file (different protonation / stereo poses); the best
(i.e. most negative) ``r_i_docking_score`` is kept.  Molecules absent from a
conformation receive an empty score.
"""

import csv
from pathlib import Path
from glob import glob

from rdkit import Chem
from rdkit.Chem import MolStandardize

SCRIPT_DIR = Path(__file__).parent
DOCKING_DIR = SCRIPT_DIR / "docking_result"

CLASSES = {
    "active": "new_molecules1_active.csv",
    "inactive": "Generate1_inactive.csv",
}

ID_CANDIDATES = ["s_sd_Title", "s_m_title"]
SMILES_CANDIDATES = ["s_sd_Original\\_SMILES", "s_sd_Original_SMILES", "Smiles"]
SCORE_COL = "r_i_docking_score"

OUTPUT_HEADER = ["ID", "Smiles", "score_1", "score_2", "score_3", "avg", "best", "n_valid"]


metal_disconnector = MolStandardize.rdMolStandardize.MetalDisconnector()
fragment_chooser = MolStandardize.rdMolStandardize.LargestFragmentChooser()
normalizer = MolStandardize.rdMolStandardize.Normalizer()

_SMILES_CACHE = {}


def clean_molecule(mol):
    mol = metal_disconnector.Disconnect(mol)  # Disconnect metals
    mol = fragment_chooser.choose(mol)        # Desalt
    mol = normalizer.normalize(mol)           # Standardize nitro, charge representation
    return mol


def clean_smiles(smiles):
    """Desalt / standardize a SMILES string; keep the original on failure."""
    if smiles in _SMILES_CACHE:
        return _SMILES_CACHE[smiles]
    cleaned = smiles
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is not None:
        cleaned_mol = clean_molecule(mol)
        if cleaned_mol is not None:
            cleaned = Chem.MolToSmiles(cleaned_mol)
    _SMILES_CACHE[smiles] = cleaned
    return cleaned


def pick_column(fieldnames, candidates, path):
    for name in candidates:
        if name in fieldnames:
            return name
    for name in fieldnames:
        for cand in candidates:
            if name.replace("\\", "") == cand.replace("\\", ""):
                return name
    raise KeyError(f"None of {candidates} found in {path}: {fieldnames}")


def read_conformation(path):
    """Return {ID: (smiles, best_score_string)} for one docking CSV."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        id_col = pick_column(reader.fieldnames, ID_CANDIDATES, path)
        smi_col = pick_column(reader.fieldnames, SMILES_CANDIDATES, path)
        if SCORE_COL not in reader.fieldnames:
            raise KeyError(f"'{SCORE_COL}' not found in {path}")

        result = {}
        for row in reader:
            mol_id = (row.get(id_col) or "").strip()
            score = (row.get(SCORE_COL) or "").strip()
            if not mol_id or not score:
                continue
            try:
                value = float(score)
            except ValueError:
                continue
            prev = result.get(mol_id)
            if prev is None or value < prev[1]:
                result[mol_id] = ((row.get(smi_col) or "").strip(), value, score)
    return result


def process_class(class_name):
    conformations = []
    for idx in (1, 2, 3):
        folder = DOCKING_DIR / f"{class_name}_state_{idx}"
        files = sorted(glob(str(folder) + "/*.csv"))
        if not files:
            raise FileNotFoundError(f"No CSV found in {folder}")
        conformations.append(read_conformation(files[0]))
        print(f"  state_{idx}: {len(conformations[-1])} molecules ({Path(files[0]).name})")

    ids = set()
    for conf in conformations:
        ids.update(conf.keys())

    rows = []
    for mol_id in ids:
        smiles = ""
        best = None
        scores = []
        for conf in conformations:
            entry = conf.get(mol_id)
            if entry is None:
                scores.append("")
                continue
            scores.append(entry[2])
            if best is None or entry[1] < best:
                best = entry[1]
                smiles = entry[0]
        valid = [float(s) for s in scores if s != ""]
        n_valid = len(valid)
        avg = str(round(sum(valid) / n_valid, 4)) if n_valid else ""
        best_str = str(round(min(valid), 4)) if n_valid else ""
        rows.append((best, mol_id, smiles, scores, avg, best_str, n_valid))

    rows.sort(key=lambda r: (r[0] is None, r[0] if r[0] is not None else 0.0, r[1]))

    out_path = SCRIPT_DIR / CLASSES[class_name]
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(OUTPUT_HEADER)
        for _, mol_id, smiles, scores, avg, best_str, n_valid in rows:
            writer.writerow([mol_id, clean_smiles(smiles), *scores, avg, best_str, n_valid])
    print(f"  -> {out_path}: {len(rows)} molecules")


def main():
    for class_name in CLASSES:
        print(f"[{class_name}]")
        process_class(class_name)


if __name__ == "__main__":
    main()
