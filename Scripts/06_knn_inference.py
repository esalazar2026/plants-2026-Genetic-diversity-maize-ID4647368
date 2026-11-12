#!/usr/bin/env python3
"""
06_knn_inference.py
k-nearest-neighbor inference of population type and race for accessions of unknown population type
(Tables S8 and S9), based on the BIO-R distance matrix.

For each unclassified accession, accessions of known type are ranked by increasing distance and the
category is assigned by majority vote among the k nearest (k = 3, 5, 10); ties are resolved in favor
of the closest neighbor. Race is inferred only for accessions proposed as landraces, using the same
neighbor set; neighbors with race UNKNOWN or NO CLASIFICABLE are kept but do not vote.
Confidence: High = same category for k = 3, 5, 10 and support >= 7/10; Medium = one criterion;
Low = neither, or support < 2/10. The category at k = 10 is retained.
Note: the race of LP1186 was set to Zapalote Chico by expert criterion (passport origin) in Table S9.

Inputs: --dist PCOA_1157acc_MR_BIOR.xlsx  --geno SNP_Filtrados_plants.xlsx  --meta STR_maiz_ConMetadato_...xlsx
Outputs: knn_population_type.csv, knn_race.csv
"""
import argparse, numpy as np, pandas as pd, openpyxl

def sample_header(path):
    ws = openpyxl.load_workbook(path, read_only=True)["SNPs"]
    rows = [r for r in ws.iter_rows(min_row=1, max_row=6, values_only=True)]
    return [str(v) for v in rows[0][36:]], [str(v) for v in rows[5][36:]]

def vote(labels, k):
    lab = labels[:k]; cnt = {}
    for l in lab:
        if l is not None: cnt[l] = cnt.get(l, 0) + 1
    if not cnt: return None, 0
    mx = max(cnt.values()); tied = [l for l in cnt if cnt[l] == mx]
    if len(tied) > 1: return next(l for l in lab if l in tied), mx
    return tied[0], mx

def confidence(res):
    cats = [c for c, _ in res]
    if res[2][1] < 2: return "Low"
    crit = (len(set(cats)) == 1 and cats[0] is not None) + (res[2][1] >= 7)
    return {2: "High", 1: "Medium", 0: "Low"}[crit]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dist", required=True); ap.add_argument("--geno", required=True); ap.add_argument("--meta", required=True)
    a = ap.parse_args()
    raw = pd.read_excel(a.dist, sheet_name="Matriz de distancia MR_BIOR", header=None)
    ids = raw.iloc[3:, 0].astype(str).tolist(); M = raw.iloc[3:, 1:].apply(pd.to_numeric).values.astype(float)
    dart, names = sample_header(a.geno); n2d = dict(zip(names, dart))
    m = pd.read_excel(a.meta, sheet_name="Metadata").set_index("DART.CODE").loc[[n2d[i] for i in ids]]
    pt = m.POPULATION_TYPE.astype(str).map({"LANDRACE": "Landrace", "BRED MATERIAL": "Breeding Material"}).fillna("UNKNOWN").values
    race = m.RACE.astype(str).values
    unk = np.where(pt == "UNKNOWN")[0]; known = np.where(pt != "UNKNOWN")[0]
    T, R = [], []
    for u in unk:
        order = known[np.argsort(M[u, known], kind="stable")]
        rt = [vote([pt[i] for i in order], k) for k in (3, 5, 10)]; ptype = rt[2][0]
        T.append([ids[u]] + [f"{c} ({n}/{k})" for (c, n), k in zip(rt, (3, 5, 10))] + [confidence(rt), ptype])
        rl = [None if race[i] in ("UNKNOWN", "NO CLASIFICABLE") else race[i] for i in order]
        rr = [vote(rl, k) for k in (3, 5, 10)]
        R.append([ids[u]] + [f"{c} ({n}/{k})" if c else "" for (c, n), k in zip(rr, (3, 5, 10))] +
                 ([confidence(rr), rr[2][0]] if ptype == "Landrace" else ["", "N/A (Breeding material)"]))
    T = pd.DataFrame(T, columns=["Accession", "k=3", "k=5", "k=10", "Confidence", "Proposed type"])
    R = pd.DataFrame(R, columns=["Accession", "k=3", "k=5", "k=10", "Confidence", "Proposed race"])
    T.to_csv("knn_population_type.csv", index=False); R.to_csv("knn_race.csv", index=False)
    lr = T["Proposed type"] == "Landrace"
    print(f"{len(unk)} unclassified accessions; type confidence {T.Confidence.value_counts().to_dict()}; "
          f"proposed {T['Proposed type'].value_counts().to_dict()}")
    print(f"Race confidence among proposed landraces: {R[lr].Confidence.value_counts().to_dict()}")
    print(f"Proposed races: {R[lr]['Proposed race'].value_counts().to_dict()}")

if __name__ == "__main__":
    main()
