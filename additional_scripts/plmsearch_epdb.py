#!/usr/bin/env python

import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO


def normalize_rows(x: np.ndarray) -> np.ndarray:
    """L2-normalize each row, guarding against zero vectors."""
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return x / norms


def load_epdb(epdb_embeds_path: str, epdb_meta_path: str):
    emb = np.load(epdb_embeds_path)
    meta = pd.read_csv(epdb_meta_path, sep="\t")

    col_map = {}
    if "seq_id" not in meta.columns and "new_id" in meta.columns:
        col_map["new_id"] = "seq_id"
    if "marker_gene" not in meta.columns and "gene_family" in meta.columns:
        col_map["gene_family"] = "marker_gene"
    if "taxon_id" not in meta.columns and "species_code" in meta.columns:
        col_map["species_code"] = "taxon_id"
    if "clan_ids" not in meta.columns and "clans" in meta.columns:
        col_map["clans"] = "clan_ids"

    if col_map:
        meta = meta.rename(columns=col_map)

    required = {"row_idx", "seq_id", "marker_gene"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(
            f"EPDB meta file {epdb_meta_path} is missing required columns: {missing}"
        )

    meta = meta.sort_values("row_idx").reset_index(drop=True)
    if not np.array_equal(meta["row_idx"].to_numpy(), np.arange(len(meta))):
        raise ValueError(
            "EPDB meta row_idx is not contiguous 0..N-1 after sorting; "
            "make sure row_idx matches embedding rows."
        )

    if emb.shape[0] != len(meta):
        raise ValueError(
            f"EPDB embeddings rows ({emb.shape[0]}) != meta rows ({len(meta)})"
        )

    return emb.astype(np.float32), meta


def build_gene_index(epdb_embeds: np.ndarray, epdb_meta: pd.DataFrame):
    """Group EPDB rows by marker_gene and build centroids + normalized vectors."""
    gene_to_idxs = defaultdict(list)
    for row in epdb_meta.itertuples(index=False):
        gene = getattr(row, "marker_gene")
        idx = getattr(row, "row_idx")
        if pd.isna(gene):
            continue
        gene_to_idxs[gene].append(idx)

    for g in list(gene_to_idxs.keys()):
        gene_to_idxs[g] = np.array(gene_to_idxs[g], dtype=np.int32)

    emb_norm = normalize_rows(epdb_embeds)

    # Centroids per gene
    gene_centroids = {}
    for gene, idxs in gene_to_idxs.items():
        centroid = epdb_embeds[idxs].mean(axis=0, keepdims=True)
        gene_centroids[gene] = normalize_rows(centroid)[0]

    return gene_to_idxs, emb_norm, gene_centroids


def load_query(query_embeds_path: str, query_index_path: str):
    q_emb = np.load(query_embeds_path)
    idx_df = pd.read_csv(query_index_path, sep="\t")

    # normalize column names
    if "row_idx" not in idx_df.columns:
        idx_df.insert(0, "row_idx", range(len(idx_df)))

    if "seq_id" not in idx_df.columns:
        # assume first non-row_idx column is id
        candidates = [c for c in idx_df.columns if c != "row_idx"]
        if not candidates:
            raise ValueError(
                f"Cannot infer seq_id column from {query_index_path};"
                " please add a 'seq_id' column."
            )
        idx_df = idx_df.rename(columns={candidates[0]: "seq_id"})

    idx_df = idx_df.sort_values("row_idx").reset_index(drop=True)
    if not np.array_equal(idx_df["row_idx"].to_numpy(), np.arange(len(idx_df))):
        raise ValueError(
            "Query index row_idx is not contiguous 0..N-1; "
            "make sure it matches embedding rows."
        )

    if q_emb.shape[0] != len(idx_df):
        raise ValueError(
            f"Query embeddings rows ({q_emb.shape[0]}) != index rows ({len(idx_df)})"
        )

    return q_emb.astype(np.float32), idx_df


def plmsearch(
    epdb_embeds,
    epdb_meta,
    gene_to_idxs,
    epdb_norm,
    gene_centroids,
    query_embeds,
    query_meta,
    top_families=16,
    family_sim_thresh=0.15,
    top_hits=8,
    hit_sim_thresh=0.25,
):
    """Perform 2-stage PLM search and return hits DataFrame and candidate IDs."""
    # Prepare centroid matrix
    genes = sorted(gene_centroids.keys())
    G = np.stack([gene_centroids[g] for g in genes], axis=0)
    G_norm = G  # centroids already normalized

    q_norm = normalize_rows(query_embeds)
    query_ids = query_meta["seq_id"].to_numpy()

    hits = []

    for qi, qv in enumerate(q_norm):
        q_id = query_ids[qi]

        # 1) family-level similarities
        sims_fam = G_norm @ qv  # shape (G,)
        order = np.argsort(-sims_fam)

        # keep top N families and above threshold
        chosen_fams = []
        for idx in order[:top_families]:
            if sims_fam[idx] < family_sim_thresh:
                continue
            chosen_fams.append((genes[idx], float(sims_fam[idx])))

        if not chosen_fams:
            continue

        # 2) within-family hits
        for gene, fam_sim in chosen_fams:
            ref_idxs = gene_to_idxs[gene]  # indices into epdb_embeds/meta
            ref_vecs = epdb_norm[ref_idxs]  # (n_g, D)

            sims = ref_vecs @ qv  # (n_g,)
            if sims.size == 0:
                continue

            # top K hits within this family
            if sims.size <= top_hits:
                loc_order = np.argsort(-sims)
            else:
                # argpartition + sort topK
                top_loc = np.argpartition(sims, -top_hits)[-top_hits:]
                loc_order = top_loc[np.argsort(-sims[top_loc])]

            for li in loc_order:
                sim = float(sims[li])
                if sim < hit_sim_thresh:
                    continue

                ref_global_idx = int(ref_idxs[li])
                ref_row = epdb_meta.iloc[ref_global_idx]

                hits.append(
                    {
                        "query_id": q_id,
                        "query_row_idx": qi,
                        "marker_gene": ref_row["marker_gene"],
                        "ref_seq_id": ref_row["seq_id"],
                        "ref_row_idx": int(ref_row["row_idx"]),
                        "taxon_id": ref_row.get("taxon_id", None),
                        "family_sim": fam_sim,
                        "hit_sim": sim,
                    }
                )

    hits_df = pd.DataFrame(hits)
    if hits_df.empty:
        return hits_df, set()

    candidate_ids = set(hits_df["query_id"].unique())
    return hits_df, candidate_ids


def write_candidates_fasta(query_fasta_path: str, candidate_ids, out_fasta_path: str):
    keep = set(candidate_ids)
    records = []
    for rec in SeqIO.parse(query_fasta_path, "fasta"):
        if rec.id in keep or rec.name in keep or rec.description.split()[0] in keep:
            records.append(rec)
    if not records:
        print("[plmsearch_epdb] WARNING: no candidate sequences written.")
    else:
        print(f"[plmsearch_epdb] Writing {len(records)} candidate sequences to {out_fasta_path}")
    SeqIO.write(records, out_fasta_path, "fasta")


def main():
    ap = argparse.ArgumentParser(
        description="PLM-based search of query proteome against EPDB embeddings."
    )
    ap.add_argument("--epdb-embeds", required=True, help="EPDB embeddings .npy")
    ap.add_argument("--epdb-meta", required=True, help="EPDB metadata TSV")
    ap.add_argument("--query-embeds", required=True, help="Query embeddings .npy")
    ap.add_argument("--query-index", required=True, help="Query index TSV")
    ap.add_argument("--query-fasta", required=True, help="Original query proteome FASTA")
    ap.add_argument("--out-hits", required=True, help="Output TSV of hits")
    ap.add_argument("--out-candidates", required=True, help="Output FASTA of candidate query proteins")
    ap.add_argument("--top-families", type=int, default=16)
    ap.add_argument("--family-sim-thresh", type=float, default=0.15)
    ap.add_argument("--top-hits", type=int, default=8)
    ap.add_argument("--hit-sim-thresh", type=float, default=0.25)

    args = ap.parse_args()

    print("[plmsearch_epdb] Loading EPDB embeddings + meta...")
    epdb_emb, epdb_meta = load_epdb(args.epdb_embeds, args.epdb_meta)

    print("[plmsearch_epdb] Building gene index + centroids...")
    gene_to_idxs, epdb_norm, gene_centroids = build_gene_index(epdb_emb, epdb_meta)

    print("[plmsearch_epdb] Loading query embeddings + index...")
    q_emb, q_meta = load_query(args.query_embeds, args.query_index)

    print("[plmsearch_epdb] Running PLM search...")
    hits_df, candidate_ids = plmsearch(
        epdb_emb=epdb_emb,
        epdb_meta=epdb_meta,
        gene_to_idxs=gene_to_idxs,
        epdb_norm=epdb_norm,
        gene_centroids=gene_centroids,
        query_embeds=q_emb,
        query_meta=q_meta,
        top_families=args.top_families,
        family_sim_thresh=args.family_sim_thresh,
        top_hits=args.top_hits,
        hit_sim_thresh=args.hit_sim_thresh,
    )

    print(f"[plmsearch_epdb] Found {len(hits_df)} hits.")
    hits_df.to_csv(args.out_hits, sep="\t", index=False)

    print("[plmsearch_epdb] Extracting candidate FASTA...")
    write_candidates_fasta(args.query_fasta, candidate_ids, args.out_candidates)

    print("[plmsearch_epdb] Done.")


if __name__ == "__main__":
    main()
