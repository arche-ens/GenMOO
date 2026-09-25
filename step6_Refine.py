#!/usr/bin/env python
# coding: utf-8
"""Step 5: Refine (transfer learning) on the filtered molecules.

Continues training the LSTM on ``new_molecules_filtered.txt`` (the best
molecules selected by the external scoring script), nudging the generation
distribution toward higher-scoring regions of chemical space.

The refined weights are written back to ``network.pth`` so that the next
generate -> score -> refine round continues from the improved model. A copy is
also saved to ``network_refined.pth``.

Inputs
------
new_molecules_filtered.txt : selected molecules (one SMILES per line)
vocab.npy                  : vocabulary from step1
network.pth                : current model weights

Outputs
-------
network_refined.pth: copy of the refined weights
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import selfies as sf
import argparse
from pathlib import Path
from tqdm import tqdm
import random
import config
import common


def loader(mols_file, library):
    with open(mols_file, "r") as fin:
        unique_smiles = [s.strip() for s in fin.readlines()]
    print(f"{len(unique_smiles):,} molecules loaded.")

    with open(library, "r") as fin:
        lib_smiles = [s.strip() for s in fin.readlines()]
        rnd_lib_smiles = random.sample(lib_smiles, len(unique_smiles)//2)
    
    # merged_unique_smiles = [s for pair in zip(unique_smiles, rnd_lib_smiles) for s in pair]
    merged_unique_smiles = unique_smiles + rnd_lib_smiles
    random.shuffle(merged_unique_smiles)
    print(f"{len(merged_unique_smiles):,} molecules in total after merged.")

    unique_selfies = []
    for s in merged_unique_smiles:
        try:
            unique_selfies.append(sf.encoder(s))
        except:
            print("[WARNING] Cannot convert to SELFIES:", s)
    print(f"{len(unique_selfies):,} molecules are valid.")
    return unique_selfies


def main(args):
    device = common.get_device()
    print("Device:", device)

    if not args.mols.exists():
        raise FileNotFoundError(
            f"Filtered file not found: {args.mols}\nRun the scoring script first and point it to write {args.mols}.")

    vocab = common.load_vocab(args.vocab)
    V = len(vocab)
    pad_idx = vocab.index(config.PAD_TOKEN)
    print("Vocabulary size:", V, "| PAD index:", pad_idx)

    selfies = loader(args.mols, args.library)
    int_kept_selfies, kept = common.encode_molecules(selfies, vocab, config.SEQ_LEN)
    print(f"{len(kept):,} molecules are encoded.")

    data = torch.from_numpy(int_kept_selfies).long()
    N, L = data.shape
    print(f"Training data: {N} sequences x {L} steps")

    model = common.CharLSTM(V, config.HIDDEN_SIZE, config.NUM_LAYERS, config.DROPOUT).to(device)
    model.load_state_dict(torch.load(args.model, map_location=device, weights_only=True))

    total_params = sum(p.numel() for p in model.parameters())
    print(f"Trainable parameters: {total_params}")

    optimizer = torch.optim.Adam(model.parameters(), lr=config.REFINE_LR)
    criterion = nn.CrossEntropyLoss(ignore_index=pad_idx)

    bsize = min(config.BATCH_SIZE, N)
    n_batches = N // bsize
    print(f"Refine batch size: {bsize} | batches per epoch: {n_batches} | epochs: {config.REFINE_EPOCHS} | LR: {config.REFINE_LR}")

    for epoch in range(1, config.REFINE_EPOCHS + 1):
        perm = torch.randperm(N)
        model.train()
        running = 0.0

        for b in tqdm(range(n_batches)):
            idx = perm[b * bsize:(b + 1) * bsize]
            batch = data[idx].to(device)

            x = batch[:, :-1]
            y = batch[:, 1:]
            x_oh = F.one_hot(x, V).float()

            optimizer.zero_grad()
            logits = model(x_oh)
            loss = criterion(logits.reshape(-1, V), y.reshape(-1))
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), config.GRAD_CLIP)
            optimizer.step()

            running += loss.item()

        print(f"Refine epoch {epoch}/{config.REFINE_EPOCHS} | average loss {running / n_batches:.4f}")


    torch.save(model.state_dict(), args.output)
    print(f"Saved refined model to {args.output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--vocab", default=config.VOCAB_FILE, type=Path, 
                        help=f"Path to vocab file (npy). Default: {config.VOCAB_FILE}")
    parser.add_argument("--model", default=config.MODEL_FILE, type=Path, 
                        help=f"Path to original model file (pth). Default: {config.MODEL_FILE}")
    parser.add_argument("-o", "--output", default=config.REFINED_MODEL_FILE, type=Path, 
                        help=f"Path to refined model file (pth). Default: {config.REFINED_MODEL_FILE}")
    parser.add_argument("--mols", default=config.FILTERED_FILE, type=Path, 
                        help=f"Filtered molecules file (txt). Default: {config.FILTERED_FILE}")
    parser.add_argument("--library", default=config.LIBRARY_FILE, type=Path,
                        help=f"Path to original library that randomly mixed with filtered data. Default: {config.KEPT_SMILES}")
    args = parser.parse_args()

    main(args)
