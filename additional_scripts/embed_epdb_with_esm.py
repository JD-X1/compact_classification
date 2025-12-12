#!/usr/bin/env python3
"""
Embed protein sequences with ESM (esm2 / esm1b) and write:
  1) out_embeds: .npy matrix (N x D) float32
  2) out_index:  .tsv with row_idx, seq_id, seq_len, num_chunks

Key features:
- CUDA or CPU
- Optional DataParallel
- Chunking for long sequences (> max_residues), with overlap
- Writes embeddings via np.lib.format.open_memmap to avoid RAM blowups
"""

import argparse
import gzip
import os
import sys
import time
from typing import Iterator, List, Tuple, Optional

import numpy as np
import torch
from Bio import SeqIO
import esm


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def count_fasta_records(path: str) -> int:
    n = 0
    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def chunk_sequence(seq: str, max_len: int, overlap: int) -> List[str]:
    """Return list of chunks each <= max_len. Overlap is number of aa overlapping between consecutive chunks."""
    L = len(seq)
    if L <= max_len:
        return [seq]
    if overlap >= max_len:
        raise ValueError("chunk_overlap must be smaller than max_residues.")

    chunks = []
    start = 0
    step = max_len - overlap
    while start < L:
        end = min(start + max_len, L)
        chunks.append(seq[start:end])
        if end == L:
            break
        start += step
    return chunks


def effective_chunk_weights(chunk_lengths: List[int], overlap: int) -> List[float]:
    """
    Approximate "unique contribution" weights when using overlaps.
    This prevents overlapped regions from being over-counted too aggressively.

    For 1 chunk: weight = full length
    For >=2 chunks: first/last lose half overlap; middle lose full overlap.
    """
    if len(chunk_lengths) == 1:
        return [float(chunk_lengths[0])]

    w = []
    for i, L in enumerate(chunk_lengths):
        if i == 0 or i == len(chunk_lengths) - 1:
            w.append(max(1.0, float(L) - overlap / 2.0))
        else:
            w.append(max(1.0, float(L) - overlap))
    return w


def get_core_model(model: torch.nn.Module) -> torch.nn.Module:
    """Handle DataParallel wrapping."""
    return model.module if isinstance(model, torch.nn.DataParallel) else model


def mean_embed_from_tokens(
    reps: torch.Tensor,
    seq_len: int,
) -> torch.Tensor:
    """
    reps: (B, T, D) representations at a layer
    seq_len: length in amino acids (no BOS/EOS)
    We average token embeddings across residues only (exclude BOS/EOS).
    """
    # token positions: 1..seq_len inclusive (BOS at 0)
    return reps[1 : seq_len + 1].mean(dim=0)


def flush_batch(
    model: torch.nn.Module,
    device: torch.device,
    batch_converter,
    batch_items: List[Tuple[int, str, int, int]],
    # (orig_idx, chunk_seq, chunk_len, chunk_id_within_orig)
    accum: dict,
    overlap: int,
    layer: int,
):
    """
    Run one forward pass batch and accumulate per-original-seq weighted chunk means.
    """
    # Build batch for ESM converter: list[(label, seq)]
    # labels can be anything stable; we keep orig idx + chunk id
    data = [(f"{orig_idx}:{chunk_i}", chunk_seq) for (orig_idx, chunk_seq, _, chunk_i) in batch_items]
    labels, strs, toks = batch_converter(data)

    toks = toks.to(device=device, non_blocking=True)

    with torch.no_grad():
        out = model(toks, repr_layers=[layer], return_contacts=False)
        reps = out["representations"][layer]  # (B, T, D)

    # Accumulate weighted chunk means per original sequence
    for b, (orig_idx, chunk_seq, chunk_len, chunk_i) in enumerate(batch_items):
        # reps[b] shape (T, D)
        vec = mean_embed_from_tokens(reps[b], seq_len=chunk_len).detach().float().cpu()

        # We'll weight later; store chunk vec and length
        if orig_idx not in accum:
            accum[orig_idx] = {"vecs": [], "lens": []}
        accum[orig_idx]["vecs"].append(vec)
        accum[orig_idx]["lens"].append(chunk_len)


def main():
    ap = argparse.ArgumentParser(description="Embed proteins with ESM and write .npy + index TSV.")
    ap.add_argument("--fasta", required=True, help="Protein FASTA (can be .gz).")
    ap.add_argument("--model", required=True, help="ESM model name, e.g. esm2_t33_650M_UR50D")
    ap.add_argument("--batch_size", type=int, default=32, help="Number of chunks per forward pass.")
    ap.add_argument("--max_residues", type=int, default=1022, help="Max aa per chunk (<=1022 for ESM2/ESM1b).")
    ap.add_argument("--chunk_overlap", type=int, default=256, help="Overlap aa between chunks for long sequences.")
    ap.add_argument("--out_embeds", required=True, help="Output .npy embeddings matrix (N x D).")
    ap.add_argument("--out_index", required=True, help="Output TSV mapping row_idx to seq_id, seq_len, num_chunks.")
    ap.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"], help="Compute device.")
    ap.add_argument("--data_parallel", action="store_true", help="Use torch.nn.DataParallel if >1 GPU available.")
    ap.add_argument("--progress_every", type=int, default=1000, help="Progress update frequency (sequences).")
    args = ap.parse_args()

    if args.max_residues > 1022:
        print(f"[WARN] max_residues={args.max_residues} > 1022. ESM models typically max at 1022 residues.", file=sys.stderr)

    # Decide device
    if args.device == "auto":
        use_cuda = torch.cuda.is_available()
    else:
        use_cuda = (args.device == "cuda")

    device = torch.device("cuda" if use_cuda else "cpu")

    print(f"Loading model {args.model} on {device}...")
    model, alphabet = esm.pretrained.load_model_and_alphabet(args.model)
    model.eval()

    if use_cuda:
        model = model.to(device)

    if args.data_parallel and use_cuda and torch.cuda.device_count() > 1:
        print(f"Enabling DataParallel across {torch.cuda.device_count()} GPUs")
        model = torch.nn.DataParallel(model)

    core = get_core_model(model)
    # For ESM2/ESM1b this exists; if not, fall back to last layer = 33-ish only if present
    layer = getattr(core, "num_layers", None)
    if layer is None:
        raise ValueError("Model object does not have num_layers; cannot choose representation layer safely.")

    embed_dim = getattr(core, "embed_dim", None)
    if embed_dim is None:
        # Some older models may store as 'args.embed_dim'
        embed_dim = getattr(getattr(core, "args", None), "embed_dim", None)
    if embed_dim is None:
        raise ValueError("Could not determine embedding dimension (embed_dim).")

    # Batch converter with truncation length (we chunk so this should not truncate)
    batch_converter = alphabet.get_batch_converter(truncation_seq_length=args.max_residues)

    n_seqs = count_fasta_records(args.fasta)
    print(f"Found {n_seqs} sequences in {args.fasta}")

    # Prepare output embedding matrix as .npy memmap
    os.makedirs(os.path.dirname(args.out_embeds) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_index) or ".", exist_ok=True)

    embeds = np.lib.format.open_memmap(
        args.out_embeds, mode="w+", dtype=np.float32, shape=(n_seqs, int(embed_dim))
    )

    # Write index header
    with open(args.out_index, "w", encoding="utf-8") as idx_out:
        idx_out.write("row_idx\tseq_id\tseq_len\tnum_chunks\n")

    t0 = time.time()
    total_chunked = 0
    max_chunks_seen = 1

    # We'll accumulate chunk vectors per original sequence index, flush immediately per sequence
    # Since we process sequences in order, we can aggregate chunks per sequence without holding too much.
    with open(args.out_index, "a", encoding="utf-8") as idx_out:
        row_idx = 0

        # Create a small buffer for chunk batches
        batch_items: List[Tuple[int, str, int, int]] = []
        accum = {}  # orig_idx -> {vecs: [Tensor], lens: [int]}

        def finalize_sequence(orig_idx: int, seq_id: str, seq_len: int):
            """Combine chunk embeddings for one original sequence and write into embeds[row_idx]."""
            nonlocal row_idx

            info = accum.get(orig_idx, None)
            if info is None:
                return  # shouldn't happen

            vecs = info["vecs"]
            lens = info["lens"]
            num_chunks = len(vecs)
            weights = effective_chunk_weights(lens, args.chunk_overlap)
            wsum = float(sum(weights))

            # weighted mean of chunk means
            stacked = torch.stack(vecs, dim=0)  # (C, D)
            w = torch.tensor(weights, dtype=torch.float32).view(-1, 1)
            combined = (stacked * w).sum(dim=0) / wsum

            embeds[row_idx, :] = combined.numpy().astype(np.float32)

            idx_out.write(f"{row_idx}\t{seq_id}\t{seq_len}\t{num_chunks}\n")
            row_idx += 1

            # free
            del accum[orig_idx]

        # Stream fasta, one original sequence at a time
        for rec_i, rec in enumerate(SeqIO.parse(open_text(args.fasta), "fasta")):
            seq_id = rec.id
            seq = str(rec.seq)
            seq_len = len(seq)

            chunks = chunk_sequence(seq, max_len=args.max_residues, overlap=args.chunk_overlap)
            if len(chunks) > 1:
                total_chunked += 1
                max_chunks_seen = max(max_chunks_seen, len(chunks))

            # Add all chunks to batching queue, and flush batches as needed
            orig_idx = rec_i
            for ci, ch in enumerate(chunks):
                batch_items.append((orig_idx, ch, len(ch), ci))

                if len(batch_items) >= args.batch_size:
                    flush_batch(model, device, batch_converter, batch_items, accum, args.chunk_overlap, layer)
                    batch_items = []

            # Flush any remaining chunks for this sequence if they're stuck in buffer
            # (we need them embedded before we can finalize this sequence)
            if batch_items:
                flush_batch(model, device, batch_converter, batch_items, accum, args.chunk_overlap, layer)
                batch_items = []

            finalize_sequence(orig_idx, seq_id, seq_len)

            if (row_idx % args.progress_every) == 0:
                dt = time.time() - t0
                rate = row_idx / max(dt, 1e-9)
                print(f"[{row_idx}/{n_seqs}] {rate:.2f} seq/s | chunked {total_chunked} | max_chunks {max_chunks_seen}")

    dt = time.time() - t0
    print(f"Done. Wrote embeddings: {args.out_embeds}")
    print(f"Done. Wrote index:      {args.out_index}")
    print(f"Total sequences: {n_seqs} | Sequences chunked: {total_chunked} | max chunks for any sequence: {max_chunks_seen}")
    print(f"Elapsed: {dt/60:.2f} min")


if __name__ == "__main__":
    main()
