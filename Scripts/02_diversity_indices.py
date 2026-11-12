#!/usr/bin/env python3
"""
02_diversity_indices.py
Genetic diversity indices (Table 1), shared/exclusive alleles between landraces and breeding
materials, and race-private alleles. Optionally repeats the frequency-based indices with
read-based within-pool allele frequencies (output of 01_read_based_frequencies.py).

Inputs
  --geno  SNP_Filtrados_plants.xlsx          (DArT 1-row matrix; 0 = hom. ref, 1 = het, 2 = hom. alt)
  --meta  STR_maiz_ConMetadato_80per_structure_only.xlsx (sheet 'Metadata': DART.CODE, POPULATION_TYPE, RACE)
  --reads read_freq.npz (optional)

Frequency coding: operational dosage 0, 0.5, 1 (alternative allele). Per locus, p = mean across
accessions; GD = 1 - p^2 - q^2; MAF = min(p, q); H = -(p ln p + q ln q); Ho = proportion of
heterozygous calls; all averaged across loci. Population type follows the metadata file
(LP731 = landrace).
"""
import argparse, numpy as np, pandas as pd

def load_genotypes(path):
    x = pd.read_excel(path, sheet_name="SNPs", header=None)
    dart = x.iloc[0, 36:].astype(str).values
    G = x.iloc[6:, 36:].apply(pd.to_numeric, errors="coerce").values.astype(float)
    return dart, G

def indices(P, G, idx):
    X = P[:, idx]; p = np.nanmean(X, 1); q = 1 - p
    ok = ~np.isnan(p); p, q = p[ok], q[ok]
    xlogx = lambda v: np.where(v > 0, v * np.log(np.where(v > 0, v, 1)), 0)
    g = G[:, idx]
    return dict(MAF=np.minimum(p, q).mean(), GD=(1 - p**2 - q**2).mean(),
                Ho=np.nanmean(np.sum(g == 1, 1) / np.sum(~np.isnan(g), 1)),
                H=(-(xlogx(p) + xlogx(q))).mean(), missing=100 * np.isnan(g).mean())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geno", required=True); ap.add_argument("--meta", required=True); ap.add_argument("--reads")
    a = ap.parse_args()
    dart, G = load_genotypes(a.geno)
    meta = pd.read_excel(a.meta, sheet_name="Metadata").set_index("DART.CODE").loc[dart]
    ptype = meta.POPULATION_TYPE.astype(str).values; race = meta.RACE.astype(str).values
    sets = {"Overall": np.arange(len(dart)), "Landrace": np.where(ptype == "LANDRACE")[0],
            "Breeding material": np.where(ptype == "BRED MATERIAL")[0]}
    codings = {"dosage 0/0.5/1": G / 2}
    if a.reads:
        z = np.load(a.reads, allow_pickle=True)
        assert list(z["dart"]) == list(dart), "Sample order differs between files"
        codings["read-based"] = z["P"].astype(float)
    for name, P in codings.items():
        print(f"\n== Diversity indices ({name})")
        for k, idx in sets.items():
            s = indices(P, G, idx)
            print(f"  {k:18s} n={len(idx):5d}  MAF={s['MAF']:.3f}  GD={s['GD']:.3f}  Ho={s['Ho']:.3f}  H={s['H']:.3f}  missing={s['missing']:.1f}%")
    # shared / exclusive alternative alleles
    lr, bm = sets["Landrace"], sets["Breeding material"]
    alt_lr = np.nansum(G[:, lr] > 0, 1) > 0; alt_bm = np.nansum(G[:, bm] > 0, 1) > 0
    print(f"\nAlternative alleles shared: {np.sum(alt_lr & alt_bm)}; only landraces: {np.sum(alt_lr & ~alt_bm)}; "
          f"only breeding materials: {np.sum(~alt_lr & alt_bm)}; absent from both: {np.sum(~alt_lr & ~alt_bm)}")
    # race-private alleles among standard races (landraces only)
    std = np.where((ptype == "LANDRACE") & ~np.isin(race, ["UNKNOWN", "NO CLASIFICABLE"]))[0]
    races = sorted(set(race[std]))
    carr = np.column_stack([np.nansum(G[:, std[race[std] == r]] > 0, 1) > 0 for r in races])
    nr = carr.sum(1)
    print(f"\nRace-private alleles: {len(std)} landraces in {len(races)} standard races; "
          f"loci with alt allele in >=2 races: {np.sum(nr >= 2)}, in none: {np.sum(nr == 0)}, private: {np.sum(nr == 1)}")
    for i in np.where(nr == 1)[0]:
        print(f"  SNP row {i + 1}: {races[int(np.argmax(carr[i]))]}")

if __name__ == "__main__":
    main()
