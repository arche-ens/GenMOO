#!/usr/bin/env python
# coding: utf-8
"""Step 3: Generate molecules.

Samples new molecules from the trained LSTM (internally as SELFIES symbols,
then decoded to SMILES) and writes them to ``new_molecules.txt`` (one SMILES
per line). This file is then scored/filtered by an external script into
``new_molecules_filtered.txt``.

Inputs
------
vocab.npy   : vocabulary from step1
network.pth : trained (or refined) model weights

Outputs
-------
new_molecules.txt : NUM_MOLECULES unique generated SMILES strings
"""

import argparse
import torch
import torch.nn.functional as F
from pathlib import Path
import config
import common
from utils.calc_novelty import calc_novelty


def set_seed():
    torch.manual_seed(config.SEED)


def sample_molecule(model, vocab, device, temperature):
    """Sample a single molecule: generate SELFIES symbols, decode to SMILES."""
    V = len(vocab)
    start_idx = vocab.index(config.START_TOKEN)
    end_idx = vocab.index(config.END_TOKEN)
    pad_idx = vocab.index(config.PAD_TOKEN)
    max_steps = config.SEQ_LEN - 2

    model.eval()
    x = torch.zeros(1, 1, V, device=device)
    x[0, 0, start_idx] = 1.0
    hidden = None
    symbols = []

    with torch.no_grad():
        for _ in range(max_steps):
            out, hidden = model.lstm(x, hidden)
            logits = model.linear(out[:, -1, :]) / temperature
            logits[0, pad_idx] = float("-inf")
            logits[0, start_idx] = float("-inf")
            probs = F.softmax(logits, dim=-1)
            nxt = torch.multinomial(probs, 1).item()

            if nxt == end_idx:
                break

            symbols.append(vocab[nxt])
            x = torch.zeros(1, 1, V, device=device)
            x[0, 0, nxt] = 1.0

    return common.symbols_to_smiles(symbols)


def main(args):
    set_seed()
    device = common.get_device()
    print("Device:", device)

    if not args.model.exists():
        raise FileNotFoundError(f"Model not found: {args.model}\nRun step2_Train.py first."
        )

    vocab = common.load_vocab(args.vocab)
    V = len(vocab)

    model = common.CharLSTM(V, config.HIDDEN_SIZE, config.NUM_LAYERS, config.DROPOUT).to(device)
    model.load_state_dict(torch.load(args.model, map_location=device, weights_only=True))
    model.eval()

    seen = set()
    count = 0
    attempts = 0
    max_attempts = args.num * 10

    with open(args.output.with_suffix(".txt"), "w") as f, open(args.output.with_suffix(".csv"), "w") as f2:
        f2.write("ID,Smiles\n")
        while count < args.num and attempts < max_attempts:
            attempts += 1
            smi = sample_molecule(model, vocab, device, args.temp)
            if smi and smi not in seen:
                seen.add(smi)
                f.write(smi + "\n")
                f2.write(f"GenMOO-{attempts:05d},{smi}\n")
                count += 1
                if count % 2000 == 0:
                    print(f"Generated {count}/{args.num} unique molecules")

    print(f"Wrote {count} unique molecules to {args.output} (total attempts: {attempts})")

    try:
        novelty = calc_novelty(args.output, config.KEPT_SMILES)
        print(f"{novelty*100:.2f}% of generated molecules are novel.")
    except:
        print("[WARNING] Cannot calculate novelty")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--vocab", default=config.VOCAB_FILE, type=Path, 
                        help=f"Path to vocab file (npy). Default: {config.VOCAB_FILE}")
    parser.add_argument("--model", default=config.MODEL_FILE, type=Path, 
                        help=f"Path to original model file (pth). Default: {config.MODEL_FILE}")
    parser.add_argument("-o", "--output", default=config.GEN_FILE,type=Path, 
                        help=f"Path to generated molecules file (txt). Default: {config.GEN_FILE}")
    parser.add_argument("--num", default=config.NUM_MOLECULES, type=int, 
                        help=f"Number of generated molecules. Default: {config.NUM_MOLECULES}")
    parser.add_argument("--temp", default=config.TEMPERATURE, type=float, 
                        help=f"Temperature!! Default: {config.TEMPERATURE}")
    args = parser.parse_args()

    main(args)
