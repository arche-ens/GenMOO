#!/usr/bin/env python
# coding: utf-8
"""Shared configuration for the GenMOO molecule-generation pipeline.

Every script in this pipeline reads its paths and hyperparameters from here so
that they stay consistent across the generate -> score -> refine loop.
"""

from pathlib import Path

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
TRAIN_DIR = BASE_DIR / "train"
RESULT_DIR = BASE_DIR / "results"
for dir in [DATA_DIR, TRAIN_DIR, RESULT_DIR]:
    if not dir.exists():
        dir.mkdir()

# --------------------------------------------------------------------------
# File paths
# --------------------------------------------------------------------------
# Input library (one SMILES per line). Provided by the user.
LIBRARY_FILE = DATA_DIR / "ChemDiv.txt"

# Cleaned SMILES (deduplicated, length-filtered). Saved by step1.
KEPT_SMILES = LIBRARY_FILE.parent / (LIBRARY_FILE.stem + "_kept.txt")

# Vocabulary (index -> character) and integer-encoded training sequences.
VOCAB_FILE = TRAIN_DIR / "vocab.npy"
INT_DATA = TRAIN_DIR / "intsmiles.npz"

# Model weights. step2 trains it; step5 refines it in place so the loop
# (generate -> score -> refine) always continues from the latest model.
MODEL_FILE = TRAIN_DIR / "network.pth"
REFINED_MODEL_FILE = TRAIN_DIR / "network_refined.pth"

# Generated molecules (output of step3) and the filtered subset produced by
# the user's external scoring script (input of step5).
GEN_FILE = RESULT_DIR / "new_molecules.txt"
FILTERED_FILE = DATA_DIR / "new_molecules_filtered.txt"

# --------------------------------------------------------------------------
# Special tokens
# --------------------------------------------------------------------------
# SELFIES is a token-level representation: every ``[Atom]``/``[Ring1]``/
# ``[Branch1_2]`` symbol is a single vocabulary unit (see common.py).
# These marker strings never occur as SELFIES symbols, so they are safe to
# use as structural markers.
START_TOKEN = "<start>"   # marks the beginning of a molecule
END_TOKEN = "<end>"       # marks the end of a molecule
PAD_TOKEN = "<pad>"       # padding for fixed-length training sequences

# --------------------------------------------------------------------------
# Model hyperparameters (mirrors the MoleculeMO LSTM)
# --------------------------------------------------------------------------
# Maximum sequence length (in SELFIES symbols, including start/end tokens).
# Molecules whose SELFIES is longer than SEQ_LEN - 2 symbols are dropped.
SEQ_LEN = 100

BATCH_SIZE = 128
HIDDEN_SIZE = 1024
NUM_LAYERS = 3
DROPOUT = 0.2

LEARNING_RATE = 0.001
GRAD_CLIP = 3.0

# --------------------------------------------------------------------------
# Training (step2)
# --------------------------------------------------------------------------
NUM_EPOCHS = 15
LOG_EVERY = 100

# --------------------------------------------------------------------------
# Generation (step3)
# --------------------------------------------------------------------------
SEED = 23333
NUM_MOLECULES = 100
TEMPERATURE = 1.0

# --------------------------------------------------------------------------
# Refinement / transfer learning (step5)
# --------------------------------------------------------------------------
REFINE_LR = 0.0001
REFINE_EPOCHS = 5
