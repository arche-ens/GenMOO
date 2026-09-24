# GenMOO2 — Multiobjective De Novo Molecule Generation (SELFIES)

A molecular-generation pipeline with the **same architecture and workflow as
`GenMOO`**, but trained and generated on **SELFIES** instead of SMILES. SELFIES
is a 100% robust molecular string representation: any sequence of SELFIES
symbols decodes to a chemically valid molecule, which nearly eliminates invalid
outputs.

```
generate  ──►  score & select (your step4)  ──►  refine  ──►  generate ...
```

## What changed vs `GenMOO`

| Aspect | `GenMOO` (SMILES) | `GenMOO2` (SELFIES) |
|--------|-------------------|---------------------|
| Representation | character-level SMILES | **token-level SELFIES** (`sf.split_selfies`) |
| Vocabulary unit | single character | one `[...]` symbol (`[C]`, `[Ring1]`, …) |
| Special tokens | `G` / `\n` / `\t` | `<start>` / `<end>` / `<pad>` |
| Validity | ~70–90% (model must learn grammar) | **~100%** (guaranteed by SELFIES) |
| Extra dependency | — | `selfies` |

The **LSTM architecture and the whole loop are unchanged** (`step2_Train.py`
and `step5_Refine.py` are identical to `GenMOO`). Only `config.py`,
`common.py`, `step1_Data_prepare.py` and `step3_Generate_molecule.py` differ.

## Boundary format stays SMILES

Step3 samples SELFIES symbols internally and **decodes them to SMILES** before
writing `new_molecules.txt`. Step5 reads SMILES from
`new_molecules_filtered.txt` and re-encodes to SELFIES. Therefore your external
scoring / Pareto-selection script (step4) works **unchanged** — it still reads
and writes SMILES.

## Pipeline

| Script | Role |
|--------|------|
| `step1_Data_prepare.py` | Clean + deduplicate SMILES, convert to SELFIES, build the symbol vocabulary and integer-encode |
| `step2_Train.py` | Train the stacked-LSTM language model from scratch (identical to GenMOO) |
| `step3_Generate_molecule.py` | Sample SELFIES symbols → decode to SMILES → write `new_molecules.txt` |
| *(external)* | Your scoring script selects the best subset → `new_molecules_filtered.txt` |
| `step5_Refine.py` | Fine-tune on the selected SMILES (re-encoded to SELFIES), then loop back to step3 |

## Model

- 3 stacked LSTM layers (`hidden_size=1024`, `dropout=0.2`) + a linear head.
- **Token-level** SELFIES modelling (one-hot input over the symbol vocabulary).
- Forget-gate bias initialized to `1`; cross-entropy with padded-position
  masking; Adam (`lr=0.001`) + gradient clipping (`3.0`).
- Generation: temperature softmax + multinomial sampling, with `<pad>` and
  `<start>` symbols masked out (`-inf`) so they are never emitted.

## Requirements

- Python 3, `numpy`, `torch`, `rdkit`, `pandas`, **`selfies`**

Recommended environment (verified: `selfies 2.1.1`, `rdkit 2024.03.5`,
`torch 2.4.1`):

```bash
/home/shanyq/miniconda3/envs/EquiScore/bin/python
```

## Quick start

```bash
cd GenMOO2

# 1. Prepare data (SMILES -> SELFIES vocabulary + encoding)
python step1_Data_prepare.py

# 2. Train the LSTM
python step2_Train.py

# 3. Generate molecules (written as SMILES)
python step3_Generate_molecule.py

# 4. (external) score new_molecules.txt -> new_molecules_filtered.txt

# 5. Refine on the selected molecules
python step5_Refine.py

# loop back to step 3 to continue the cycle
```

## Data flow

```
data/test_data.txt            (input library, one SMILES per line)
  └─ step1 ─► data/test_data_kept.txt        (cleaned SMILES)
               train/vocab.npy               (index -> SELFIES symbol)
               train/intsmiles.npz           (padded SELFIES symbol indices)
  └─ step2 ─► train/network.pth              (trained weights)
  └─ step3 ─► new_molecules.txt              (generated SMILES, decoded)
  └─ step4 ─► new_molecules_filtered.txt     (selected subset, by your script)
  └─ step5 ─► train/network.pth (updated) + train/network_refined.pth
```

## Configuration (`config.py`)

Same parameters as `GenMOO`; the differences are the special tokens and the
comment on `SEQ_LEN` (now counted in SELFIES **symbols**):

- `START_TOKEN = "<start>"`, `END_TOKEN = "<end>"`, `PAD_TOKEN = "<pad>"`.
- `SEQ_LEN = 100` (SELFIES symbols, incl. start/end).
- `NUM_MOLECULES = 20`, `TEMPERATURE = 1.0`, etc. — editable as needed.

## Notes

- SELFIES guarantees **syntactic + valence validity**, not chemical
  *meaningfulness*; the external scoring + selection loop (step4/step5) is still
  required to steer the generator toward high-quality molecules.
- The SMILES variant lives in `../GenMOO` for side-by-side comparison.
