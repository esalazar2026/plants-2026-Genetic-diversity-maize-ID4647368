"""
Pool-size-weighted pairwise FST among genetic groups (Methods 4.9, "Population Distances").

Fst_snp = Var_weighted(p)_among_groups / [ p_bar * (1 - p_bar) ]

where, for a given pair of groups (A, B) and a given SNP:
  - p_i      = bulk-level dosage allele frequency of accession i (0, 0.5, or 1)
  - w_i      = BULK_PLANTS_NUM of accession i (number of plants pooled in the bulk)
  - p_g      = w_i-weighted mean of p_i within group g (g in {A, B})
  - p_bar    = w_i-weighted mean of p_i across A U B combined
  - W_g      = sum of w_i within group g
  - Var_weighted(p)_among_groups = sum_g [ W_g * (p_g - p_bar)^2 ] / sum_g W_g

The locus-level Fst is averaged across all 20,956 quality-filtered SNPs (loci with
p_bar in {0, 1}, i.e. monomorphic for the pair, are skipped: the denominator is 0).

Reproduces Table 2 (pairwise FST among STRUCTURE groups, K = 2 and K = 6).
"""
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- load genotypes
geno = np.load('/home/claude/geno_rebuilt.npy')                       # (20956, 1157), -1 = missing, 0/1/2
saga_codes = np.load('/home/claude/saga_codes_rebuilt.npy', allow_pickle=True)
saga_codes = np.array([str(s).strip() for s in saga_codes])

dose = geno.astype(np.float64)
dose[dose == -1] = np.nan
dose = dose / 2.0                                                      # 0, 0.5, 1

# ---------------------------------------------------------------- metadata + group assignments
# IMPORTANT: join on SAGA.CODE only. LP.BD.MAIZ codes are inconsistent between source files
# (e.g. "BB12" vs "XBB12", "LP670-22" vs no counterpart) and produce silent join failures;
# SAGA.CODE is clean and matches the genotype matrix exactly (0 unmatched).
str_mtd = pd.read_excel('/mnt/user-data/uploads/STR_maiz_ConMetadato_80per_structure_only.xlsx',
                         sheet_name='Metadata')
str_mtd['SAGA.CODE'] = str_mtd['SAGA.CODE'].astype(str).str.strip()


def assign_k2(row):
    q1, q2 = row['K2(1)'], row['K2(2)']
    if pd.isna(q1) or pd.isna(q2):
        return 'NA'
    if q1 >= 0.80:
        return 'Group 1'
    elif q2 >= 0.80:
        return 'Group 2'
    else:
        return 'Admixed'


def assign_k6(row):
    qs = [row[f'K6({i})'] for i in range(1, 7)]
    if any(pd.isna(q) for q in qs):
        return 'NA'
    mx = max(qs)
    if mx > 0.50:
        return f'Group {qs.index(mx) + 1}'
    else:
        return 'Admixed'


str_mtd['K2_Group'] = str_mtd.apply(assign_k2, axis=1)
str_mtd['K6_Group'] = str_mtd.apply(assign_k6, axis=1)

saga_to_w = dict(zip(str_mtd['SAGA.CODE'], str_mtd['BULK_PLANTS_NUM']))
saga_to_k2 = dict(zip(str_mtd['SAGA.CODE'], str_mtd['K2_Group']))
saga_to_k6 = dict(zip(str_mtd['SAGA.CODE'], str_mtd['K6_Group']))

weights = np.array([saga_to_w.get(s, np.nan) for s in saga_codes], dtype=np.float64)
k2_group = np.array([saga_to_k2.get(s, 'NA') for s in saga_codes], dtype=object)
k6_group = np.array([saga_to_k6.get(s, 'NA') for s in saga_codes], dtype=object)


def weighted_fst(mask_a, mask_b, dose, weights):
    """Pool-size-weighted Fst between two groups of accessions, averaged over loci."""
    idx = np.where(mask_a | mask_b)[0]
    wA = weights[mask_a]
    wB = weights[mask_b]
    WA, WB = np.nansum(wA), np.nansum(wB)

    p = dose[:, idx]                       # (n_snps, n_acc_pair)
    w_pair = weights[idx]
    is_a = mask_a[idx]

    valid = ~np.isnan(p)
    w_mat = np.broadcast_to(w_pair, p.shape).copy()
    w_mat[~valid] = 0.0
    p_filled = np.nan_to_num(p, nan=0.0)

    # group-weighted means per locus
    wA_mat = np.where(is_a, w_mat, 0.0)
    wB_mat = np.where(~is_a, w_mat, 0.0)
    WA_loc = wA_mat.sum(axis=1)
    WB_loc = wB_mat.sum(axis=1)
    pA = np.divide(( wA_mat * p_filled).sum(axis=1), WA_loc, out=np.full(p.shape[0], np.nan), where=WA_loc > 0)
    pB = np.divide(( wB_mat * p_filled).sum(axis=1), WB_loc, out=np.full(p.shape[0], np.nan), where=WB_loc > 0)

    W_loc = WA_loc + WB_loc
    pbar = np.divide(WA_loc * pA + WB_loc * pB, W_loc, out=np.full(p.shape[0], np.nan), where=W_loc > 0)

    var_among = (WA_loc * (pA - pbar) ** 2 + WB_loc * (pB - pbar) ** 2) / W_loc
    denom = pbar * (1 - pbar)

    valid_loci = (~np.isnan(denom)) & (denom > 0) & (~np.isnan(var_among))
    fst_per_locus = var_among[valid_loci] / denom[valid_loci]

    return np.mean(fst_per_locus), valid_loci.sum(), int(WA), int(WB), int(mask_a.sum()), int(mask_b.sum())


# ---------------------------------------------------------------- K = 2
maskA = k2_group == 'Group 1'
maskB = k2_group == 'Group 2'
fst, nloci, WA, WB, nA, nB = weighted_fst(maskA, maskB, dose, weights)
print(f"K=2  Group1 (n={nA}) vs Group2 (n={nB})  Fst = {fst:.3f}   (n_loci used = {nloci})")

# ---------------------------------------------------------------- K = 6
k6_groups = ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5', 'Group 6']
print("\nK=6 pairwise Fst matrix:")
results = {}
for i in range(len(k6_groups)):
    for j in range(i + 1, len(k6_groups)):
        gi, gj = k6_groups[i], k6_groups[j]
        maskA = k6_group == gi
        maskB = k6_group == gj
        fst, nloci, WA, WB, nA, nB = weighted_fst(maskA, maskB, dose, weights)
        results[(gi, gj)] = fst
        print(f"  {gi} vs {gj} (n={nA},{nB}):  Fst = {fst:.3f}")
