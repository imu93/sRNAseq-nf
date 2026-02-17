#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pysam
import polars as pl

_CATS = [
    "gene",
    "pseudogene",
    "DNA",
    "LTR",
    "LINE",
    "Unknown",
    "Other_repeat",
    "miRNA",
    "rRNA",
    "piRNA",
    "snRNA",
    "snoRNA",
    "tRNA",
    "lincRNA",
    "Other_ncRNA",
    "Unannotated",
    "Other",
]

_RULES = [
    ("gene",         re.compile(r"^gene", re.IGNORECASE)),
    ("pseudogene",   re.compile(r"pseudogene", re.IGNORECASE)),
    ("DNA",          re.compile(r"DNA|RC", re.IGNORECASE)),
    ("LTR",          re.compile(r"LTR", re.IGNORECASE)),
    ("LINE",         re.compile(r"LINE|PLE|Retro", re.IGNORECASE)),
    ("Unknown",      re.compile(r"Unknown", re.IGNORECASE)),
    ("miRNA",        re.compile(r"^miRNA", re.IGNORECASE)),
    ("piRNA",        re.compile(r"^piRNA", re.IGNORECASE)),
    ("snRNA",        re.compile(r"^snRNA", re.IGNORECASE)),
    ("snoRNA",       re.compile(r"^snoRNA", re.IGNORECASE)),
    ("lincRNA",      re.compile(r"^lincRNA", re.IGNORECASE)),
    ("rRNA",         re.compile(r"^rRNA|_rRNA|LSU", re.IGNORECASE)),
    ("tRNA",         re.compile(r"^tRNA", re.IGNORECASE)),
    ("Other_repeat", re.compile(r"Sat|Low|Simple|SINE", re.IGNORECASE)),
    ("Other_ncRNA",  re.compile(r"^yRNA|^lncRNA|other_ncRNA|circRNA", re.IGNORECASE)),
    ("Unannotated",  re.compile(r"inter|Unassigned", re.IGNORECASE)),
]


def classify_feature_id(feature_id: str) -> str:
    if not feature_id:
        return "Unannotated"
    for cat, pat in _RULES:
        if pat.search(feature_id):
            return cat
    return "Other"


def pick_best_category_from_xt(xt_value: str) -> str:
    if not xt_value:
        return "Unannotated"
    targets = [t.strip() for t in xt_value.split(",") if t.strip()]
    if not targets:
        return "Unannotated"

    best_cat = "Other"
    best_rank = 10**9
    for t in targets:
        cat = classify_feature_id(t)
        if cat == "Other":
            rank = 10**8
        else:
            rank = next((i for i, (c, _) in enumerate(_RULES) if c == cat), 10**7)
        if rank < best_rank:
            best_rank = rank
            best_cat = cat
    return best_cat


def get_rc(read: pysam.AlignedSegment, rc_tag: str) -> int:
    try:
        v = read.get_tag(rc_tag)
        if isinstance(v, int) and v > 0:
            return v
    except KeyError:
        pass
    return 1


def infer_sample_name(bam_path: str) -> str:
    name = Path(bam_path).name
    for suf in [
        ".expanded.bam.featureCounts.bam",
        ".featureCounts.bam",
        ".expanded.bam",
        ".bam",
    ]:
        if name.endswith(suf):
            name = name[: -len(suf)]
            break
    if name.endswith(".bam"):
        name = name[:-4]
    return name


def summarize_bam(
    bam_path: str,
    minlen: int,
    maxlen: int,
    all_lengths: bool,
    rc_tag: str,
    require_assigned: bool,
) -> Tuple[Dict[int, Dict[str, int]], Tuple[int, int]]:
    counts: Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    min_obs = None
    max_obs = None

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
                continue

            seq = read.query_sequence
            if not seq:
                continue

            L = len(seq)
            if not all_lengths and (L < minlen or L > maxlen):
                continue

            min_obs = L if min_obs is None else min(min_obs, L)
            max_obs = L if max_obs is None else max(max_obs, L)

            try:
                status = read.get_tag("XS")
            except KeyError:
                status = None

            if require_assigned and status is not None and str(status) != "Assigned":
                continue

            try:
                xt = str(read.get_tag("XT"))
            except KeyError:
                xt = ""

            cat = pick_best_category_from_xt(xt) if xt else "Unannotated"
            w = get_rc(read, rc_tag)
            counts[L][cat] += w

    if min_obs is None:
        if all_lengths:
            min_obs, max_obs = 0, 0
        else:
            min_obs, max_obs = minlen, maxlen

    return counts, (int(min_obs), int(max_obs))


def counts_to_tsv(
    counts: Dict[int, Dict[str, int]],
    length_range: Tuple[int, int],
    out_tsv: str,
):
    minL, maxL = length_range
    if maxL < minL:
        pl.DataFrame({"length": []}).write_csv(out_tsv, separator="\t")
        return

    lengths = list(range(minL, maxL + 1))
    rows = []
    for L in lengths:
        row = {"length": L}
        for c in _CATS:
            row[c] = int(counts.get(L, {}).get(c, 0))
        rows.append(row)

    pl.DataFrame(rows).sort("length").write_csv(out_tsv, separator="\t")


def main():
    ap = argparse.ArgumentParser(
        description="Build a matrix of counts per feature and lengths (one per input BAM file)"
    )
    ap.add_argument("--bam", required=True, nargs="+", help="input BAM")
    ap.add_argument("--minlen", type=int, default=18)
    ap.add_argument("--maxlen", type=int, default=27)
    ap.add_argument("--all-lengths", action="store_true")
    ap.add_argument("--rc-tag", default="RC")
    ap.add_argument("--require-assigned", action="store_true")
    ap.add_argument("--suffix", default=".cls_mtx.tsv")
    args = ap.parse_args()

    for bam_path in args.bam:
        sample = infer_sample_name(bam_path)

        counts, (min_obs, max_obs) = summarize_bam(
            bam_path=bam_path,
            minlen=args.minlen,
            maxlen=args.maxlen,
            all_lengths=args.all_lengths,
            rc_tag=args.rc_tag,
            require_assigned=args.require_assigned,
        )

        out_range = (min_obs, max_obs) if args.all_lengths else (args.minlen, args.maxlen)
        out_file = f"{sample}{args.suffix}"
        counts_to_tsv(counts, out_range, out_file)


if __name__ == "__main__":
    main()
