#!/usr/bin/env python3
"""
05_structure_consistency.py
Consistency of STRUCTURE membership coefficients across independent runs (Table S3).

For each K, the Q matrices of all runs (*_f output files) are aligned to the reported run (highest
'Mean value of ln likelihood') by testing all K! label permutations and keeping the one that
maximizes the symmetric similarity coefficient H' (Jakobsson & Rosenberg 2007, CLUMPP):
  H'(Q1, Q2) = 1 - ||Q1 - Q2||_F / sqrt(||Q1 - W||_F * ||Q2 - W||_F),  W = matrix of 1/K.
Reports mean pairwise H', SD of Q across runs and the number of accessions with identical group
assignment in all runs (K = 2: Q >= 0.80; K > 2: maximum Q > 0.50; otherwise admixed).

Usage: python 05_structure_consistency.py --runs output_STR_K2/*_f   (repeat for K6)
"""
import argparse, itertools, re, numpy as np

def parse(path):
    txt = open(path).read().splitlines(); head = "\n".join(txt[:120])
    g = lambda pat: float(re.search(pat, head).group(1))
    i0 = [i for i, l in enumerate(txt) if l.startswith("Inferred ancestry of individuals")][0] + 2
    labs, Q = [], []
    for l in txt[i0:]:
        mm = re.match(r"\s*\d+\s+(\S+)\s+\(\d+\)\s*:\s*(.*)", l)
        if not mm: break
        labs.append(mm.group(1)); Q.append([float(v) for v in mm.group(2).split()])
    par = [l for l in txt if l.startswith("DATAFILE=")][0]
    return dict(lnp=g(r"Estimated Ln Prob of Data\s*=\s*(\S+)"), mll=g(r"Mean value of ln likelihood\s*=\s*(\S+)"),
                labs=labs, Q=np.array(Q), burnin=re.search(r"BURNIN=(\d+)", par).group(1),
                numreps=re.search(r"NUMREPS=(\d+)", par).group(1))

def Hp(a, b):
    W = np.full_like(a, 1 / a.shape[1])
    return 1 - np.linalg.norm(a - b) / np.sqrt(np.linalg.norm(a - W) * np.linalg.norm(b - W))

def assign(Q):
    K = Q.shape[1]; mx = Q.max(1)
    return np.where(mx >= 0.8, Q.argmax(1), -1) if K == 2 else np.where(mx > 0.5, Q.argmax(1), -1)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", nargs="+", required=True); a = ap.parse_args()
    runs = {f: parse(f) for f in sorted(a.runs)}
    ref = max(runs, key=lambda f: runs[f]["mll"]); K = runs[ref]["Q"].shape[1]
    assert all(r["labs"] == runs[ref]["labs"] for r in runs.values())
    print(f"K = {K}; {len(runs)} runs; burn-in {runs[ref]['burnin']}, MCMC {runs[ref]['numreps']}; reported run: {ref}")
    al = {}
    for f, r in runs.items():
        perm = max(itertools.permutations(range(K)), key=lambda p: Hp(runs[ref]["Q"], r["Q"][:, list(p)]))
        al[f] = r["Q"][:, list(perm)]
    H = [Hp(al[x], al[y]) for x, y in itertools.combinations(al, 2)]
    lnp = np.array([r["lnp"] for r in runs.values()]); Qs = np.stack(list(al.values()))
    A = np.array([assign(Q) for Q in al.values()]); ra = assign(al[ref]); same = np.all(A == ra, 0)
    print(f"Ln P(D): mean {lnp.mean():,.1f} +/- {lnp.std(ddof=1):,.1f}")
    print(f"H' mean {np.mean(H):.4f} (range {np.min(H):.4f}-{np.max(H):.4f}); mean SD of Q across runs {Qs.std(0).mean():.4f}; "
          f"max |Q - Q_reported| {np.abs(Qs - al[ref]).max():.3f}")
    print(f"Identical assignment in all runs: {same.sum()} of {len(same)} ({100 * same.mean():.1f}%)")
    labs = np.array(runs[ref]["labs"])
    for g in sorted(set(ra)):
        w = ra == g; print(f"  {'admixed' if g < 0 else 'cluster ' + str(g + 1)}: n = {w.sum()}, stable = {same[w].sum()}")
    for i in np.where(~same)[0]:
        print(f"  unstable: {labs[i]}  max Q in reported run = {al[ref][i].max():.3f}")

if __name__ == "__main__":
    main()
