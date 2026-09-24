#!/usr/bin/env python
# coding: utf-8
"""Shared building blocks for the GenMOO pipeline.

Contains the LSTM language model (faithful to the MoleculeMO architecture)
and small helpers shared by step1, step2, step3 and step5.

This version trains/generates on SELFIES symbols (token level) instead of
SMILES characters, which guarantees that any sampled sequence decodes to a
chemically valid molecule.
"""

import numpy as np
import torch
import torch.nn as nn
import selfies as sf

import config


def get_device():
    """Return the best available torch device."""
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def load_vocab(path=config.VOCAB_FILE):
    """Load the vocabulary (list of characters) saved by step1.

    The index of a character in the returned list is its integer id.
    """
    vocab = np.load(path, allow_pickle=True).tolist()
    return [str(c) for c in vocab]


def encode_molecules(selfies_list, vocab, seq_len):
    """Encode SELFIES strings into fixed-length integer sequences.

    Each molecule is first read as SELFIES, then split into symbols, and
    finally becomes ``[START] + symbols + [END] + [PAD]*...`` of length
    ``seq_len``. Molecules that fail SELFIES encoding, contain unknown symbols,
    or exceed ``seq_len`` are skipped.

    Returns:
        (array (N, seq_len) of int, list of kept SMILES)
    """
    char_to_idx = {c: i for i, c in enumerate(vocab)}
    pad_idx = char_to_idx[config.PAD_TOKEN]
    start_idx = char_to_idx[config.START_TOKEN]
    end_idx = char_to_idx[config.END_TOKEN]

    rows = []
    kept = []
    for s in selfies_list:
        s = s.strip()
        if not s:
            continue
        try:
            symbols = list(sf.split_selfies(s))
        except:
            print("Failed to split:", s)
            continue
        try:
            seq = [start_idx] + [char_to_idx[sym] for sym in symbols] + [end_idx]
        except KeyError:
            continue
        if len(seq) > seq_len:
            continue
        seq += [pad_idx] * (seq_len - len(seq))
        rows.append(seq)
        kept.append(s)

    data = np.zeros((len(rows), seq_len), dtype=np.int64)
    for i, r in enumerate(rows):
        data[i, :] = r
    return data, kept


def symbols_to_smiles(symbols):
    """Decode a list of SELFIES symbols back into a SMILES string.

    Returns ``None`` if decoding fails (should be rare, thanks to SELFIES
    robustness).
    """
    try:
        return sf.decoder("".join(symbols))
    except Exception:
        return None


class CharLSTM(nn.Module):
    """LSTM language model for SELFIES generation.

    Stacked LSTM layers followed by a linear projection to vocabulary-size logits. 
    Input is one-hot encoded.
    Each vocabulary unit is a SELFIES symbol (e.g. ``[C]``, ``[Ring1]``). 
    """

    def __init__(self, vocab_size, hidden_size, num_layers, dropout):
        super(CharLSTM, self).__init__()
        self.vocab_size = vocab_size
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.dropout = dropout

        self.lstm = nn.LSTM(
            input_size=vocab_size,
            hidden_size=hidden_size,
            num_layers=num_layers,
            dropout=dropout,
            batch_first=True,
        )
        self.linear = nn.Linear(hidden_size, vocab_size)

    def init_forget_gates(self):
        """Initialize LSTM forget-gate biases to 1.

        This biases the network to remember information by default, which
        helps capture long-range dependencies such as matching parentheses
        and ring-closure digits in SMILES strings.
        """
        with torch.no_grad():
            for names in self.lstm._all_weights:
                for name in filter(lambda n: "bias" in n, names):
                    bias = getattr(self.lstm, name)
                    n = bias.size(0)
                    bias.data[n // 4: n // 2].fill_(1.0)

    def forward(self, x):
        out, _ = self.lstm(x)
        return self.linear(out)
