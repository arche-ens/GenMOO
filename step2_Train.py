#!/usr/bin/env python
# coding: utf-8
"""Step 2: Train the LSTM.

Trains the character-level LSTM (from scratch) on the integer-encoded SMILES
produced by step1 so that it learns to generate valid SMILES strings.

Inputs
------
vocab.npy     : vocabulary from step1
intsmiles.npz : integer-encoded sequences from step1

Outputs
-------
network.pth : trained model weights
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import config
import common


def main():
    device = common.get_device()
    print("Device:", device)

    vocab = common.load_vocab(config.VOCAB_FILE)
    V = len(vocab)
    pad_idx = vocab.index(config.PAD_TOKEN)
    print("Vocabulary size:", V, "| PAD index:", pad_idx)

    int_kept_selfies = np.load(config.INT_DATA)["arr_0"]
    data = torch.from_numpy(int_kept_selfies).long()
    N, L = data.shape
    print(f"Training data: {N} sequences x {L} steps")

    model = common.CharLSTM(V, config.HIDDEN_SIZE, config.NUM_LAYERS, config.DROPOUT).to(device)
    model.init_forget_gates()

    total_params = sum(p.numel() for p in model.parameters())
    print(f"Trainable parameters: {total_params}")

    optimizer = torch.optim.Adam(model.parameters(), lr=config.LEARNING_RATE)
    criterion = nn.CrossEntropyLoss(ignore_index=pad_idx)

    bsize = min(config.BATCH_SIZE, N)
    n_batches = N // bsize
    print(f"Batch size: {bsize} | batches per epoch: {n_batches} | epochs: {config.NUM_EPOCHS} | LR: {config.LEARNING_RATE}")

    for epoch in range(1, config.NUM_EPOCHS + 1):
        perm = torch.randperm(N)
        model.train()
        running = 0.0

        for b in range(n_batches):
            idx = perm[b * bsize:(b + 1) * bsize]
            batch = data[idx].to(device)

            x = batch[:, :-1]          # input:  chars 0..L-2
            y = batch[:, 1:]           # target: chars 1..L-1
            x_oh = F.one_hot(x, V).float()

            optimizer.zero_grad()
            logits = model(x_oh)
            loss = criterion(logits.reshape(-1, V), y.reshape(-1))
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), config.GRAD_CLIP)
            optimizer.step()

            running += loss.item()

            if (b + 1) % config.LOG_EVERY == 0:
                avg = running / (b + 1)
                print(f"  epoch {epoch}/{config.NUM_EPOCHS} | batch {b + 1}/{n_batches} | loss {loss.item():.4f} | avg {avg:.4f}")

        print(f"Epoch {epoch}/{config.NUM_EPOCHS} | average loss {running / n_batches:.4f}")
        if epoch % 10 == 0:
            torch.save(model.state_dict(), config.MODEL_FILE.parent / (config.MODEL_FILE.stem + f"_epoch{epoch}"))

    torch.save(model.state_dict(), config.MODEL_FILE)
    print("Saved model to", config.MODEL_FILE.resolve())


if __name__ == "__main__":
    main()
