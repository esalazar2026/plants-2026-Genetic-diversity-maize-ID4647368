#!/usr/bin/env python3
"""
03_fst_pool_weighted.py
Pairwise pool-size-weighted FST among STRUCTURE groups (K = 2, K = 6; Table 2) and hierarchical
clusters (HC-C3, HC-C6; Table S6), using the 0/0.5/1 dosage coding and, optionally, read-based
within-pool allele frequencies (Table S11).

For each pair of groups and each SNP: group frequencies p_g are means of bulk frequencies weighted
by the number of plants in each bulk (missing calls excluded); p_bar is the weighted mean of the
two groups; Var_w(p) = sum_g a_g (p_g - p_bar)^2 with a_g the share of plants in group g.
FST = sum_loci Var_w(p) / sum_loci p_bar(1 - p_bar)   (ratio of averages across loci).

Inputs
  --geno     SNP_Filtrados_plants.xlsx
  --meta     STR_maiz_ConMetadato_80per_structure_only.xlsx (sheet 'Metadata': DART.CODE,
             BULK_PLANTS_NUM, K2(1..2), K6(1..6) membership coefficients)
  --clusters cluster_maiz_ConMetadato_80per.xlsx (sheet 'Metadata': DART.CODE, MRogers-HC-C3, MRogers-HC-C6)
  --reads    read_freq.npz (optional)
Group assignment: K = 2, Q >= 0.80; K = 6, maximum Q > 0.50; others admixed (excluded).

Partition of total gene diversity (FST-AMOVA comparison, Section 3.2): for each SNP, with bulk-size
weights w_i and bulk frequencies p_i,
  p_bar(1 - p_bar) = among groups [sum_g W_g (p_g - p_bar)^2 / W]
                   + among accessions within groups [sum_i w_i (p_i - p_g(i))^2 / W]
                   + within accessions [sum_i w_i p_i (1 - p_i) / W],
each component summed across loci and expressed as a proportion of the summed total. The
among-group proportion equals the FST computed across all groups. The script also reports the
among / (among + among-accessions) proportion and the AMOVA-type variance-component estimate
PhiST = s2_a / (s2_a + s2_w), s2_a = (MSA - MSW) / n0, from the same (unweighted) sums of squares.
"""
import argparse, numpy as np, pandas as pd

def load_genotypes(path):
    x = pd.read_excel(path, sheet_name="SNPs", header=None)
    return x.iloc[0, 36:].astype(str).values, x.iloc[6:, 36:].apply(pd.to_numeric, errors="coerce").values.astype(float)

def fst(P, w, ga, gb):
    def grp(g):
        X = P[:, g]; W = np.where(np.isnan(X), 0, w[g][None, :])
        return (W * np.nan_to_num(X)).sum(1) / W.sum(1), W.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        p1, n1 = grp(ga); p2, n2 = grp(gb); a1 = n1 / (n1 + n2)
        pb = a1 * p1 + (1 - a1) * p2
        var = a1 * (p1 - pb) ** 2 + (1 - a1) * (p2 - pb) ** 2; den = pb * (1 - pb)
    ok = np.isfinite(den) & (den > 0)
    return var[ok].sum() / den[ok].sum()

def partition(P, lab, w):
    sel = np.where(lab != "adm")[0]; X = P[:, sel]; l = lab[sel]
    Wt = np.where(np.isnan(X), 0, w[sel][None, :]); X = np.nan_to_num(X)
    W = Wt.sum(1); pb = (Wt * X).sum(1) / W; A = np.zeros(len(W)); B = np.zeros(len(W))
    for g in np.unique(l):
        m = l == g; Wg = Wt[:, m].sum(1)
        pg = np.where(Wg > 0, (Wt[:, m] * X[:, m]).sum(1) / np.where(Wg > 0, Wg, 1), 0)
        A += Wg * (pg - pb) ** 2; B += (Wt[:, m] * (X[:, m] - pg[:, None]) ** 2).sum(1)
    A /= W; B /= W; C = (Wt * X * (1 - X)).sum(1) / W; T = pb * (1 - pb); ok = T > 0
    A, B, C, T = A[ok].sum() / T[ok].sum(), B[ok].sum() / T[ok].sum(), C[ok].sum() / T[ok].sum(), 1.0
    N = len(sel); ns = np.array([np.sum(l == g) for g in np.unique(l)]); G = len(ns); n0 = (N - np.sum(ns**2) / N) / (G - 1)
    return A, B, C, A / (A + B), N, G, n0

def phist_vc(A, B, N, G, n0):
    SSA, SSW = N * A, N * B; MSW = SSW / (N - G); sa = (SSA / (G - 1) - MSW) / n0
    return sa / (sa + MSW)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geno", required=True); ap.add_argument("--meta", required=True)
    ap.add_argument("--clusters", required=True); ap.add_argument("--reads"); a = ap.parse_args()
    dart, G = load_genotypes(a.geno)
    m = pd.read_excel(a.meta, sheet_name="Metadata").set_index("DART.CODE").loc[dart]
    hc = pd.read_excel(a.clusters, sheet_name="Metadata").set_index("DART.CODE").loc[dart]
    w = m.BULK_PLANTS_NUM.values.astype(float)
    k2 = np.where(m["K2(1)"] >= 0.8, "1", np.where(m["K2(2)"] >= 0.8, "2", "adm"))
    q6 = m[[f"K6({i})" for i in range(1, 7)]].values
    k6 = np.where(q6.max(1) > 0.5, (q6.argmax(1) + 1).astype(str), "adm")
    part = {"STRUCTURE K = 2": k2, "STRUCTURE K = 6": k6,
            "HC-C3": hc["MRogers-HC-C3"].astype(str).values, "HC-C6": hc["MRogers-HC-C6"].astype(str).values}
    codings = {"dosage": G / 2}
    if a.reads:
        z = np.load(a.reads, allow_pickle=True); assert list(z["dart"]) == list(dart)
        codings["reads"] = z["P"].astype(float)
    for pname, lab in part.items():
        groups = sorted(g for g in set(lab) if g != "adm")
        print(f"\n== {pname}: " + ", ".join(f"{g} (n={np.sum(lab == g)})" for g in groups))
        for i, ga in enumerate(groups):
            for gb in groups[i + 1:]:
                vals = "  ".join(f"{c}={fst(P, w, np.where(lab == ga)[0], np.where(lab == gb)[0]):.3f}" for c, P in codings.items())
                print(f"  {ga} vs {gb}: {vals}")
        if pname.startswith("STRUCTURE"):
            for c, P in codings.items():
                for wt, ww in (("bulk-weighted", w), ("unweighted", np.ones_like(w))):
                    A, B, C, ratio, N, G, n0 = partition(P, lab, ww)
                    print(f"  Partition ({c}, {wt}): among groups {100 * A:.1f}%, among accessions within groups {100 * B:.1f}%, "
                          f"within accessions {100 * C:.1f}%; among/(among + among-acc) = {ratio:.3f}"
                          + (f"; variance-component PhiST = {phist_vc(A, B, N, G, n0):.3f}" if wt == "unweighted" else ""))

if __name__ == "__main__":
    main()
