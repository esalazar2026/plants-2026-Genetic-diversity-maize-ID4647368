#!/usr/bin/env python3
"""
07_ld_decay.py
Linkage disequilibrium (r^2) between chromosome-anchored SNPs (B73 RefGen_v4) separated by <= 1 Mb in
the landrace subset (Table S2, Figure 2), using the 0/0.5/1 dosage coding and, optionally,
read-based within-pool allele frequencies (Table S11).
r^2 = squared Pearson correlation between SNP values across landraces; missing values are replaced
by the locus mean. Mean r^2 is reported per distance interval (left-closed bins), for all anchored
SNPs and for SNPs with minor allele frequency >= --min-maf in landraces (default 0.10). The number of
pairs in the first interval that are SNPs within the same DArTseq tag (CloneID) is also reported.

Inputs: --geno SNP_Filtrados_plants.xlsx  --meta STR_maiz_ConMetadato_...xlsx  [--reads read_freq.npz]
Output: ld_pairs_<coding>.csv (chromosome, distance, r2)
"""
import argparse, numpy as np, pandas as pd

BINS = [0, 100, 1000, 2000, 3000, 4000, 5000, 7500, 10000, 15000, 20000, 30000, 50000, 100000, 250000, 500000, 1000001]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geno", required=True); ap.add_argument("--meta", required=True); ap.add_argument("--reads")
    ap.add_argument("--min-maf", type=float, default=0.10)
    a = ap.parse_args()
    x = pd.read_excel(a.geno, sheet_name="SNPs", header=None)
    hdr = x.iloc[0, :36].tolist(); mk = x.iloc[6:, :36].copy(); mk.columns = hdr
    dart = x.iloc[0, 36:].astype(str).values; G = x.iloc[6:, 36:].apply(pd.to_numeric, errors="coerce").values.astype(float)
    m = pd.read_excel(a.meta, sheet_name="Metadata").set_index("DART.CODE").loc[dart]
    lr = np.where(m.POPULATION_TYPE.astype(str).values == "LANDRACE")[0]
    chrom = mk["Chrom_Maize_v4"].astype(str).values; pos = pd.to_numeric(mk["ChromPos_Maize_v4"], errors="coerce").values
    anch = pd.Series(chrom).str.match(r"^\d+_chromosome$").values
    codings = {"dosage": G / 2}
    if a.reads:
        z = np.load(a.reads, allow_pickle=True); assert list(z["dart"]) == list(dart); codings["reads"] = z["P"].astype(float)
    clone = mk["CloneID"].astype(str).values
    p = np.nanmean(G[:, lr] / 2, 1); maf = np.minimum(p, 1 - p)
    print(f"Landraces: {len(lr)}; chromosome-anchored SNPs: {anch.sum()}; with MAF >= {a.min_maf}: {np.sum(anch & (maf >= a.min_maf))}")
    for (name, P), mmin in [(c, m) for m in (0.0, a.min_maf) for c in codings.items()]:
        rec = []
        for c in sorted(set(chrom[anch]), key=lambda s: int(s.split("_")[0])):
            idx = np.where((chrom == c) & (maf >= mmin))[0]; idx = idx[np.argsort(pos[idx], kind="stable")]
            X = P[np.ix_(idx, lr)]; X = np.where(np.isnan(X), np.nanmean(X, 1, keepdims=True), X); pp = pos[idx]
            for i in range(len(idx)):
                for j in range(i + 1, len(idx)):
                    d = pp[j] - pp[i]
                    if d > 1_000_000: break
                    if X[i].std() == 0 or X[j].std() == 0: continue
                    rec.append((int(c.split("_")[0]), d, np.corrcoef(X[i], X[j])[0, 1] ** 2, clone[idx[i]] == clone[idx[j]]))
        R = pd.DataFrame(rec, columns=["chromosome", "distance_bp", "r2", "same_tag"])
        tag = f"{name}_maf{mmin:.2f}"; R.to_csv(f"ld_pairs_{tag}.csv", index=False)
        tab = R.groupby(pd.cut(R.distance_bp, BINS, right=False), observed=False).r2.agg(["size", "mean"])
        f = R[R.distance_bp < 100]
        print(f"\n== {name}, MAF >= {mmin:.2f}: {len(R):,} pairs; overall mean r2 = {R.r2.mean():.3f}; "
              f"pairs <= 100 bp within the same tag: {int(f.same_tag.sum())} of {len(f)} (mean r2 {f[f.same_tag].r2.mean():.3f})")
        print(tab.round(3).to_string())

if __name__ == "__main__":
    main()
