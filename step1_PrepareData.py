#!/usr/bin/env python
# coding: utf-8
"""Step 1: Data preparation.

Reads ``data.txt`` (one SMILES per line), clean and deduplicates it, filters out
molecules that are too long, converts them to SELFIES, builds the SELFIES
symbol vocabulary and encodes every molecule into a fixed-length integer
sequence ready for LSTM training.

Outputs
-------
vocab.npy     : vocabulary (index -> SELFIES symbol)
data_kept.txt : cleaned SMILES (deduplicated, length-filtered)
intsmiles.npz : integer-encoded padded SELFIES sequences, shape (N, SEQ_LEN)
"""

import numpy as np
import selfies as sf

import config
import common
from utils.clean import clean_data


def build_vocab(selfies_list):
    """Build the vocabulary: special tokens first, then sorted SELFIES symbols."""
    alphabet = sf.get_alphabet_from_selfies(selfies_list)
    return [config.PAD_TOKEN, config.START_TOKEN, config.END_TOKEN] + sorted(alphabet)


def main():
    if not config.LIBRARY_FILE.exists():
        raise FileNotFoundError(
            f"Input library not found: {config.LIBRARY_FILE}\n"
            "Place a file with one SMILES per line at this path first."
        )


    ##### read data directly #####
    with open(config.LIBRARY_FILE, "r") as fin:
        unique_smiles = [s.strip() for s in fin.readlines()]
    print(f"{len(unique_smiles):,} molecules loaded.")
    ##############################
    # If you are not sure whether your library is clean,
    # use the following code to desalt, standardize, and deduplicate:
    ##### clean the data #####
    # unique_smiles = clean_data(config.LIBRARY_FILE)
    ##########################

    unique_selfies = []
    for s in unique_smiles:
        try:
            unique_selfies.append(sf.encoder(s))
        except:
            print("[WARNING] Cannot convert to SELFIES:", s)
    print(f"{len(unique_smiles):,} molecules are valid.")
    breakpoint()

    max_len = config.SEQ_LEN - 2  # reserve room for START + END tokens
    kept_selfies = [s for s in unique_selfies if (0 < len(list(sf.split_selfies(s))) <= max_len)]
    print(f"{len(kept_selfies):,} molecules accepted for training (SELFIES symbol length <= {max_len}).")

    with open(config.KEPT_SMILES, "w") as f:
        for s in kept_selfies:
            f.write(sf.decoder(s) + "\n")

    vocab = build_vocab(kept_selfies)
    np.save(config.VOCAB_FILE, np.array(vocab, dtype=object))

    int_kept_selfies, _ = common.encode_molecules(kept_selfies, vocab, config.SEQ_LEN)
    np.savez_compressed(config.INT_DATA, int_kept_selfies)

    print(f"\nVocabulary size: {len(vocab)}")
    print(f"Encoded data shape: {int_kept_selfies.shape}")
    print("\nOutput:")
    print(config.KEPT_SMILES.resolve())
    print(config.VOCAB_FILE.resolve())
    print(config.INT_DATA.resolve())


if __name__ == "__main__":
    main()
