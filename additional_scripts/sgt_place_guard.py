#!/usr/bin/env python
"""Off-ramp for markers that yield no/empty EPA-ng placement in SGT mode.

A query gene too divergent to place produces a missing, empty, or zero-placement
.jplace. Failing the run on that (`sgt_gappa` did `exit 2`) discards the
unplaceable-marker fraction, which is itself a novelty signal
(dispersion_aware_placement_calling.md sections 4-5, 7). Instead, record the gene
as unplaceable and emit a header-only profile.tsv so the DAG completes and the
marker casts zero placement votes while staying countable downstream.
"""
import argparse
import json
import os

# gappa `examine assign --best-hit` profile.tsv header. A sentinel profile
# carries this header and no data rows: sgt_summary's `grep -v LWR` empties it
# (zero votes), and is_unplaceable_profile reports True (countable).
PROFILE_HEADER = "LWR\tfract\taLWR\tafract\ttaxopath"
FLAG_NAME = "unplaceable.flag"

# 'no_placements' is the novelty-informative flavor (EPA-ng ran on a real query
# and returned nothing); the others are weaker/mechanical (doc sections 4-5).
UNPLACEABLE_REASONS = ("no_placements", "empty_file", "missing_file",
                       "unparseable", "unknown")


def classify_jplace(path):
    """Return (status, reason). status is 'placeable' or 'unplaceable'."""
    if not os.path.exists(path):
        return "unplaceable", "missing_file"
    if os.path.getsize(path) == 0:
        return "unplaceable", "empty_file"
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (json.JSONDecodeError, ValueError, UnicodeDecodeError):
        return "unplaceable", "unparseable"
    if not data.get("placements"):
        return "unplaceable", "no_placements"
    return "placeable", ""


def write_unplaceable_profile(profile_path, gene, reason):
    """Write the header-only sentinel profile and a sibling unplaceable.flag."""
    out_dir = os.path.dirname(profile_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(profile_path, "w") as fh:
        fh.write(PROFILE_HEADER + "\n")
    with open(os.path.join(out_dir, FLAG_NAME), "w") as fh:
        fh.write(f"{gene}\t{reason}\n")


def is_unplaceable_profile(profile_path):
    """True iff the profile has no placement (data) rows."""
    if not os.path.exists(profile_path):
        return True
    with open(profile_path) as fh:
        for line in fh:
            if line.strip() and "LWR" not in line:
                return False
    return True


def count_unplaceable_markers(profile_paths):
    """Return (n_unplaceable, n_total, fraction) over a MAG's per-gene profiles."""
    n_total = len(profile_paths)
    n_unplaceable = sum(1 for p in profile_paths if is_unplaceable_profile(p))
    fraction = n_unplaceable / n_total if n_total else 0.0
    return n_unplaceable, n_total, fraction


def read_unplaceable_reason(profile_path):
    """Reason for an unplaceable marker (from its sibling flag), or None if placed.

    An unplaceable profile whose flag is missing returns 'unknown' rather than
    None, so it is never silently reclassified as placed.
    """
    if not is_unplaceable_profile(profile_path):
        return None
    flag = os.path.join(os.path.dirname(profile_path), FLAG_NAME)
    if not os.path.exists(flag):
        return "unknown"
    with open(flag) as fh:
        line = fh.readline().rstrip("\n")
    parts = line.split("\t")
    reason = parts[1] if len(parts) > 1 else "unknown"
    return reason if reason in UNPLACEABLE_REASONS else "unknown"


def summarize_unplaceable(profile_paths):
    """Per-reason unplaceable breakdown for a MAG's per-gene profiles.

    Isolates the novelty-informative 'no_placements' flavor from the weaker
    empty/missing/unparseable markers.
    """
    n_total = len(profile_paths)
    counts = {reason: 0 for reason in UNPLACEABLE_REASONS}
    n_unplaceable = 0
    for path in profile_paths:
        reason = read_unplaceable_reason(path)
        if reason is None:
            continue
        n_unplaceable += 1
        counts[reason] += 1
    summary = {
        "n_markers": n_total,
        "n_unplaceable": n_unplaceable,
        "unplaceable_fraction": n_unplaceable / n_total if n_total else 0.0,
        "no_placements_fraction": counts["no_placements"] / n_total if n_total else 0.0,
    }
    for reason in UNPLACEABLE_REASONS:
        summary[f"n_{reason}"] = counts[reason]
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)

    c = sub.add_parser("classify", help="classify one jplace; sentinel if unplaceable")
    c.add_argument("--jplace", required=True)
    c.add_argument("--gene", required=True)
    c.add_argument("--profile-out", required=True,
                   help="profile.tsv to create when the gene is unplaceable")

    t = sub.add_parser("tally", help="report the unplaceable-marker fraction for a MAG")
    t.add_argument("--profiles", nargs="*", default=[],
                   help="every per-gene profile.tsv for the MAG")
    t.add_argument("--out", required=True)
    t.add_argument("--mag", default="")

    args = ap.parse_args()

    if args.command == "classify":
        status, reason = classify_jplace(args.jplace)
        if status == "placeable":
            print("PLACEABLE")
            return
        write_unplaceable_profile(args.profile_out, args.gene, reason)
        print("UNPLACEABLE")
        return

    s = summarize_unplaceable(args.profiles)
    reason_cols = [f"n_{r}" for r in UNPLACEABLE_REASONS]
    header = (["mag", "n_markers", "n_unplaceable", "unplaceable_fraction",
               "no_placements_fraction"] + reason_cols)
    values = ([args.mag, s["n_markers"], s["n_unplaceable"],
               f"{s['unplaceable_fraction']:.6f}", f"{s['no_placements_fraction']:.6f}"]
              + [s[c] for c in reason_cols])
    with open(args.out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(str(v) for v in values) + "\n")


if __name__ == "__main__":
    main()
