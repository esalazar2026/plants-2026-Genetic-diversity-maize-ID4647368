#!/usr/bin/env python3
"""
08_distance_pcoa_sensitivity.py
(1) Verifies that the BIO-R distance matrix equals the squared Modified Rogers distance computed from
    the presence/absence of each allele: (1/2L) sum_loci [(ref_i - ref_j)^2 + (alt_i - alt_j)^2],
    with pairwise-complete loci.
(2) Compares Modified Rogers distances computed from allele frequencies with the 0/0.5/1 coding and
    with read-based within-pool frequencies (Table S11): matrix correlation and principal coordinates
    (classical PCoA, double-centering of -0.5 D^2).

Inputs: --geno SNP_Filtrados_plants.xlsx  --reads read_freq.npz  --dist PCOA_1157acc_MR_BIOR.xlsx
"""
import argparse, numpy as np, pandas as pd

def sqdist(X):
    """Mean squared difference over pairwise-complete loci; X = loci x samples (NaN = missing)."""
    M = (~np.isnan(X)).astype(float); Z = np.nan_to_num(X); L = M.T @ M
    S = (Z**2).T @ M + M.T @ (Z**2) - 2 * Z.T @ Z
    return np.maximum(S, 0) / L

def pcoa(D):
    n = len(D); J = np.eye(n) - 1 / n; Bc = J @ (-0.5 * D**2) @ J
    ev, V = np.linalg.eigh(Bc); o = np.argsort(ev)[::-1]; ev, V = ev[o], V[:, o]
    return V[:, :3] * np.sqrt(ev[:3]), 100 * ev[:3] / ev[ev > 0].sum()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geno", required=True); ap.add_argument("--reads", required=True); ap.add_argument("--dist", required=True)
    a = ap.parse_args()
    x = pd.read_excel(a.geno, sheet_name="SNPs", header=None)
    dart = x.iloc[0, 36:].astype(str).values; names = x.iloc[5, 36:].astype(str).tolist()
    G = x.iloc[6:, 36:].apply(pd.to_numeric, errors="coerce").values.astype(float)
    z = np.load(a.reads, allow_pickle=True); assert list(z["dart"]) == list(dart); Pr = z["P"].astype(float)
    raw = pd.read_excel(a.dist, sheet_name="Matriz de distancia MR_BIOR", header=None)
    ids = raw.iloc[3:, 0].astype(str).tolist(); Bm = raw.iloc[3:, 1:].apply(pd.to_numeric).values.astype(float)
    o = [names.index(i) for i in ids]; iu = np.triu_indices(len(o), 1)
    ref = np.where(np.isnan(G), np.nan, (G < 2) * 1.0)[:, o]; alt = np.where(np.isnan(G), np.nan, (G > 0) * 1.0)[:, o]
    Bp = (sqdist(ref) + sqdist(alt)) / 2
    print(f"(1) BIO-R matrix vs presence/absence squared MRD: max |difference| = {np.abs(Bp - Bm)[iu].max():.2e}, r = {np.corrcoef(Bp[iu], Bm[iu])[0, 1]:.6f}")
    Cb, vb = pcoa(np.sqrt(Bm)); print(f"    PCoA of the BIO-R matrix: variance explained by PC1-PC3 = {np.round(vb, 2)} %")
    Md = np.sqrt(sqdist((G / 2)[:, o])); Mr = np.sqrt(sqdist(Pr[:, o]))
    print(f"(2) MRD dosage vs MRD reads: r = {np.corrcoef(Md[iu], Mr[iu])[0, 1]:.3f}")
    Cd, vd = pcoa(Md); Cr, vr = pcoa(Mr)
    print(f"    PC1 |r| between codings = {abs(np.corrcoef(Cd[:, 0], Cr[:, 0])[0, 1]):.3f}; variance PC1: dosage {vd[0]:.1f}%, reads {vr[0]:.1f}%")

if __name__ == "__main__":
    main()
