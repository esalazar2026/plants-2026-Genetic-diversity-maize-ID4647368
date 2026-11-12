#!/usr/bin/env python3
"""
01_read_based_frequencies.py
Estimate within-pool allele frequencies from the DArTseq read-count file and
validate them against the 1-row genotype calls.

Inputs
  --geno    SNP_Filtrados_plants.xlsx  (sheet 'SNPs'; DArT 1-row matrix, 20,956 SNPs x 1,157 accessions;
                                         rows 1-6 = sample annotations, codes 0 = hom. reference,
                                         1 = heterozygous, 2 = hom. alternative)
  --counts  SEQ_SNPs_counts_0_Target_Erika_Chile_For_Bio.csv  (DArT 2-row read counts)
Output
  --out     read_freq.npz  (alt-allele frequency per SNP x accession, read depth, sample DArT codes)

Technical-replicate columns of the same accession (same SAGA code) are summed. The presence/absence
of reads for each allele is checked against the genotype calls.

Python 3.12, numpy, pandas, openpyxl.
"""
import argparse, numpy as np, pandas as pd

def load_genotypes(path):
    x = pd.read_excel(path, sheet_name="SNPs", header=None)
    hdr = x.iloc[0, :36].tolist()
    markers = x.iloc[6:, :36].copy(); markers.columns = hdr; markers = markers.reset_index(drop=True)
    samples = pd.DataFrame({"dart": x.iloc[0, 36:].astype(str).values, "saga": x.iloc[1, 36:].astype(str).values,
                            "lp": x.iloc[5, 36:].astype(str).values})
    G = x.iloc[6:, 36:].apply(pd.to_numeric, errors="coerce").values.astype(float)
    return markers, samples, G

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geno", required=True); ap.add_argument("--counts", required=True)
    ap.add_argument("--out", default="read_freq.npz"); a = ap.parse_args()
    markers, samples, G = load_genotypes(a.geno)
    alt_id = {k: i for i, k in enumerate(markers["AlleleID"].astype(str))}     # e.g. 2425608|F|0-39:C>A-39:C>A
    ref_id = {k: i for i, k in enumerate(markers["AlleleID.1"].astype(str))}   # e.g. 2425608|F|0--39:C>A
    L = len(markers)
    with open(a.counts) as f:
        cols = f.readline().rstrip("\n").split(",")[1:]
        ALT = np.zeros((L, len(cols)), np.int64); REF = np.zeros((L, len(cols)), np.int64)
        for line in f:
            k, rest = line.split(",", 1)
            if k in alt_id: ALT[alt_id[k]] = np.array(rest.split(","), dtype=np.int64)
            elif k in ref_id: REF[ref_id[k]] = np.array(rest.split(","), dtype=np.int64)
    cols = np.array(cols)
    A = np.zeros((L, len(samples))); R = np.zeros_like(A)
    for j, s in enumerate(samples.saga):
        idx = np.where(cols == s)[0]
        if len(idx) == 0: raise SystemExit(f"Sample {s} not found in counts file")
        A[:, j] = ALT[:, idx].sum(1); R[:, j] = REF[:, idx].sum(1)
    called = ~np.isnan(G)
    concord = (((G > 0) == (A > 0)) & ((G < 2) == (R > 0)))[called].mean()
    print(f"Concordance between genotype calls and allele presence in reads: {concord:.4f}")
    print(f"Called cells with zero reads: {np.sum(called & (A + R == 0))}; missing cells with reads: {np.sum(~called & (A + R > 0))}")
    depth = A + R
    P = np.where(called & (depth > 0), A / np.where(depth > 0, depth, 1), np.nan)
    het = G == 1; fh = P[het]
    print(f"Heterozygous calls: {het.sum():,} ({het.sum() / called.sum():.1%} of non-missing)")
    print(f"Within-pool alt-allele frequency at heterozygous calls: mean {fh.mean():.3f}, median {np.median(fh):.3f}, "
          f"IQR {np.percentile(fh, 25):.3f}-{np.percentile(fh, 75):.3f}")
    print(f"  0.25-0.75: {np.mean((fh >= .25) & (fh < .75)):.1%}; <0.25: {np.mean(fh < .25):.1%}; >0.75: {np.mean(fh >= .75):.1%}")
    print(f"Median read depth: called cells {np.median(depth[called]):.0f}; heterozygous calls {np.median(depth[het]):.0f}")
    np.savez_compressed(a.out, P=P.astype(np.float32), depth=depth.astype(np.int32), dart=samples.dart.values)
    print("Saved", a.out)

if __name__ == "__main__":
    main()
