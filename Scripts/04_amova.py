#!/usr/bin/env python3
"""
04_amova.py
AMOVA (Excoffier et al. 1992) for STRUCTURE groups (K = 2, Q >= 0.80; K = 6, max Q > 0.50) and for
the 22 standard landrace races, plus pairwise PhiST between races (Table S7) with permutation tests
and Benjamini-Hochberg FDR correction.

The distance matrix exported by BIO-R ('Matriz de distancia MR_BIOR') equals
(1/2L) * sum_loci [(ref_i - ref_j)^2 + (alt_i - alt_j)^2] computed from the presence/absence of each
allele, i.e. the squared Modified Rogers distance. It is therefore used directly as the squared
distance in AMOVA (it is NOT squared again).

Inputs
  --dist  PCOA_1157acc_MR_BIOR.xlsx (sheet 'Matriz de distancia MR_BIOR')
  --geno  SNP_Filtrados_plants.xlsx (only the sample annotation rows are read, to link accession
          names in the distance matrix with DArT codes)
  --meta  STR_maiz_ConMetadato_80per_structure_only.xlsx (sheet 'Metadata')
  --perm-global 999  --perm-pair 9999  --seed 1  --out pairwise_PhiST_races.csv
"""
import argparse, numpy as np, pandas as pd, openpyxl

def sample_header(path):
    ws = openpyxl.load_workbook(path, read_only=True)["SNPs"]
    rows = [r for r in ws.iter_rows(min_row=1, max_row=6, values_only=True)]
    return [str(v) for v in rows[0][36:]], [str(v) for v in rows[5][36:]]   # DArT codes, accession names

def amova(D2, lab):
    lab = np.asarray(lab); N = len(lab); groups = np.unique(lab)
    sst = D2.sum() / 2 / N
    ssw = sum(D2[np.ix_(lab == g, lab == g)].sum() / 2 / np.sum(lab == g) for g in groups)
    G = len(groups); ns = np.array([np.sum(lab == g) for g in groups]); n0 = (N - np.sum(ns**2) / N) / (G - 1)
    sw = ssw / (N - G); sa = ((sst - ssw) / (G - 1) - sw) / n0
    return sa / (sa + sw)

def perm_test(D2, lab, nperm, rng):
    obs = amova(D2, lab); c = 1
    for _ in range(nperm):
        if amova(D2, rng.permutation(lab)) >= obs - 1e-12: c += 1
    return obs, c / (nperm + 1)

def pair_test(D2, n1, nperm, rng):
    """Fast permutation test for two groups (first n1 rows/cols = group 1)."""
    N = len(D2); n2 = N - n1; lab = np.r_[np.ones(n1), np.zeros(n2)]
    sst = D2.sum() / 2 / N; n0 = N - (n1**2 + n2**2) / N
    def phi(Z):
        Zc = 1 - Z
        ssw = np.einsum("pi,ij,pj->p", Z, D2, Z) / 2 / n1 + np.einsum("pi,ij,pj->p", Zc, D2, Zc) / 2 / n2
        sw = ssw / (N - 2); sa = ((sst - ssw) - sw) / n0; return sa / (sa + sw)
    obs = phi(lab[None, :])[0]
    Z = np.array([rng.permutation(lab) for _ in range(nperm)])
    return obs, (np.sum(phi(Z) >= obs - 1e-12) + 1) / (nperm + 1)

def bh(p):
    p = np.asarray(p); m = len(p); o = np.argsort(p); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]; return np.minimum(q, 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dist", required=True); ap.add_argument("--geno", required=True); ap.add_argument("--meta", required=True)
    ap.add_argument("--perm-global", type=int, default=999); ap.add_argument("--perm-pair", type=int, default=9999)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--out", default="pairwise_PhiST_races.csv")
    a = ap.parse_args()
    raw = pd.read_excel(a.dist, sheet_name="Matriz de distancia MR_BIOR", header=None)
    ids = raw.iloc[3:, 0].astype(str).tolist(); assert ids == raw.iloc[2, 1:].astype(str).tolist()
    D2 = raw.iloc[3:, 1:].apply(pd.to_numeric).values.astype(float)
    dart, names = sample_header(a.geno); name2dart = dict(zip(names, dart))
    m = pd.read_excel(a.meta, sheet_name="Metadata").set_index("DART.CODE").loc[[name2dart[i] for i in ids]]
    k2 = np.where(m["K2(1)"] >= 0.8, 1, np.where(m["K2(2)"] >= 0.8, 2, 0))
    q6 = m[[f"K6({i})" for i in range(1, 7)]].values; k6 = np.where(q6.max(1) > 0.5, q6.argmax(1) + 1, 0)
    ptype = m.POPULATION_TYPE.astype(str).values; race = m.RACE.astype(str).values
    std = np.where((ptype == "LANDRACE") & ~np.isin(race, ["UNKNOWN", "NO CLASIFICABLE"]))[0]
    for name, sel, lab in [("K = 2", np.where(k2 > 0)[0], k2), ("K = 6", np.where(k6 > 0)[0], k6), ("22 races", std, race)]:
        phi, p = perm_test(D2[np.ix_(sel, sel)], lab[sel], a.perm_global, np.random.default_rng(a.seed))
        print(f"AMOVA {name}: n = {len(sel)}, among groups = {100 * phi:.1f}%, within = {100 - 100 * phi:.1f}%, PhiST = {phi:.4f}, p = {p:.3f}")
    races = sorted(set(race[std])); rows = []; rng = np.random.default_rng(a.seed)
    for i, ra in enumerate(races):
        for rb in races[i + 1:]:
            ia = std[race[std] == ra]; ib = std[race[std] == rb]; idx = np.r_[ia, ib]
            phi, p = pair_test(D2[np.ix_(idx, idx)], len(ia), a.perm_pair, rng)
            rows.append((ra, rb, len(ia), len(ib), phi, p))
    R = pd.DataFrame(rows, columns=["race_a", "race_b", "n_a", "n_b", "PhiST", "p"]); R["q_FDR"] = bh(R.p)
    R.to_csv(a.out, index=False)
    ns = R[R.q_FDR >= 0.05]
    print(f"Pairwise PhiST ({len(R)} comparisons): significant after FDR {np.sum(R.q_FDR < 0.05)} ({100 * np.mean(R.q_FDR < 0.05):.1f}%); "
          f"max {R.PhiST.max():.3f}; median {R.PhiST.median():.3f}; non-significant {len(ns)}, "
          f"{np.sum((ns.n_a < 15) | (ns.n_b < 15))} involving a race with < 15 accessions")
    print(ns[(ns.n_a >= 15) & (ns.n_b >= 15)].round(4).to_string(index=False))
    allsig = [r for r in races if (R[(R.race_a == r) | (R.race_b == r)].q_FDR < 0.05).all()]
    print("Races differentiated from all others:", ", ".join(allsig)); print("Saved", a.out)

if __name__ == "__main__":
    main()
