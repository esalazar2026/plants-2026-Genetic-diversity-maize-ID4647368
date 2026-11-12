# Genetic diversity of maize accessions genotyped with DArTseq — data and scripts

Data and analysis scripts supporting the article submitted to *Plants* (MDPI), manuscript ID **plants-4647368**:

> **[Authors]**. **Genome-wide SNP genotyping reveals the genetic diversity and population structure of the Chilean national maize germplasm collection**. *Plants* **[year, volume, article number]**. https://doi.org/[DOI — add when published]

The repository contains the quality-filtered DArTseq SNP matrix for 1,157 maize accessions (bulked samples), the
accession metadata, and the Python scripts used for the revised analyses: read-based within-pool allele frequencies,
diversity indices, pool-weighted FST, AMOVA, consistency of STRUCTURE runs, k-NN inference of population type and race,
linkage-disequilibrium decay, and a sensitivity analysis of genetic distances and PCoA.

---

## Repository structure

```
README.md
LICENSE
requirements.txt
data/
├── SNP_Filtrados_plants.xlsx                          DArTseq 1-row SNP matrix (20,956 SNPs x 1,157 accessions)
├── SEQ_SNPs_counts_0_Target_Erika_Chile_For_Bio.csv   DArTseq read counts (2-row format)
├── STR_maiz_ConMetadato_80per_structure_only.xlsx     accession metadata + STRUCTURE Q (K = 2, K = 6)
├── cluster_maiz_ConMetadato_80per.xlsx                hierarchical-clustering assignments (HC-C3, HC-C6)
├── PCOA_1157acc_MR_BIOR.xlsx                          distance matrix exported from BIO-R
├── MRD_distance_matrix_1157x1157.csv                  Modified Rogers distance matrix (see notes below)
├── core_hunter/
│   ├── genotypes_biparental_20956SNP_1157acc.csv.gz   genotypes for Core Hunter (0/1/2, NA = missing)
│   └── strata_core_1157acc.csv                        stratum and forced race representative per accession
└── structure/
    ├── output_STR_K2/*_f                              STRUCTURE run outputs, K = 2
    └── output_STR_K6/*_f                              STRUCTURE run outputs, K = 6
scripts/
├── 01_read_based_frequencies.py
├── 02_diversity_indices.py
├── 03_fst_pool_weighted.py
├── 04_amova.py
├── 05_structure_consistency.py
├── 06_knn_inference.py
├── 07_ld_decay.py
├── 08_distance_pcoa_sensitivity.py
└── core_hunter/
    └── core_estratificada_corehunter.R                stratified core collection (Core Hunter 3, R)
results/                                               outputs written by the scripts
```

---

## Data files

### `SNP_Filtrados_plants.xlsx`
DArTseq SNP genotypes in DArT 1-row format, after quality filtering.

* **Sheet `SNPs`** — 20,956 SNPs (rows) x 1,157 accessions.
  * Columns 1–36: DArT marker annotation (`AlleleID`, `AlleleID.1`, `CloneID`, ..., `Chrom_Maize_v4`,
    `ChromPos_Maize_v4` (B73 RefGen_v4), ..., `rdepth`, `maf`).
  * Columns 37 onward: one column per accession.
  * Rows 1–6: sample annotation — (1) DArT sample code, (2) SAGA code, (3) population type, (4) race,
    (5) race typicity, (6) accession name (LP code).
  * Rows 7 onward: genotype codes **0** = homozygous reference, **1** = heterozygous,
    **2** = homozygous alternative, empty = missing (5.7 % of cells).
* **Sheet `freq`** — number of loci by allele-frequency class.
* **Sheet `otros`** — per-marker statistics (He, Ho, Ae, Shannon index, specific alleles, % missing).

### `SEQ_SNPs_counts_0_Target_Erika_Chile_For_Bio.csv`
DArTseq read counts in 2-row format (one row per allele, reference and alternative). Columns are SAGA codes;
technical replicates of the same accession appear as repeated columns and are summed by script 01.

### `STR_maiz_ConMetadato_80per_structure_only.xlsx`
Sheet `Metadata`: one row per accession with `DART.CODE`, `SAGA.CODE`, `POPULATION_TYPE` (LANDRACE, BRED MATERIAL,
or unknown), `RACE`, `BULK_PLANTS_NUM` (number of plants in each bulk), and STRUCTURE membership coefficients
`K2(1)`–`K2(2)` and `K6(1)`–`K6(6)`.

Group assignment used throughout: K = 2, Q ≥ 0.80; K = 6, maximum Q > 0.50; other accessions are treated as admixed.

### `cluster_maiz_ConMetadato_80per.xlsx`
Sheet `Metadata`: `DART.CODE` and hierarchical-cluster assignments `MRogers-HC-C3` and `MRogers-HC-C6`.

### `PCOA_1157acc_MR_BIOR.xlsx`
Sheet `Matriz de distancia MR_BIOR`: 1,157 x 1,157 distance matrix exported from BIO-R. It equals the **squared**
Modified Rogers distance computed from the presence/absence of each allele,
(1/2L) Σ [(ref_i − ref_j)² + (alt_i − alt_j)²], over pairwise-complete loci (verified by script 08).
Used by scripts 04, 06 and 08.

### `MRD_distance_matrix_1157x1157.csv`
1,157 x 1,157 Modified Rogers distance matrix (semicolon-delimited; row and column names = accession names).
It equals sqrt[(1/L) Σ (x_i − x_j)²] computed from the 0/0.5/1 dosage coding of `SNP_Filtrados_plants.xlsx`, with
missing genotypes replaced by the locus mean across all 1,157 accessions. **This is a different matrix from the BIO-R
matrix above** and is not read by any script in this repository. [Add: analysis in which it was used.]

Accession names in this file differ from those in `SNP_Filtrados_plants.xlsx` for 66 accessions. The correspondence
(confirmed by matching the distance profiles exactly) is:

| In `MRD_distance_matrix_1157x1157.csv` | In `SNP_Filtrados_plants.xlsx` |
|---|---|
| XBB1 – XBB58 | BB1 – BB58 (same number) |
| XBB59, XBB60, XBB61, XBB62, XBB63 | LP670-1, LP670-7, LP670-16, LP670-18, LP670-22 |
| XBB64, XBB65, XBB66, XBB67, XBB68 | LP670-27, LP670-46, LP670-51, LP670-52, LP670-55 |
| LP3289 | LP3290 |

### `structure/output_STR_K2/`, `structure/output_STR_K6/`
STRUCTURE output files (`*_f`) of all independent runs for K = 2 and K = 6. Used by script 05.

---

## Software requirements

* Python 3.12 — numpy, pandas, openpyxl (see `requirements.txt`)
* R with the `corehunter` package and Java 8 or later (core collection only): `install.packages("corehunter")`

```bash
pip install -r requirements.txt
```

---

## How to reproduce the analyses

Run the scripts from the `results/` folder, so that output files are written there. Script 01 must run first: it
creates `read_freq.npz`, which scripts 02, 03, 07 and 08 use (optional for 02, 03 and 07; required for 08).
Results are printed to the screen; redirect them to a log file to keep them.

```bash
cd results

python ../scripts/01_read_based_frequencies.py --geno ../data/SNP_Filtrados_plants.xlsx \
    --counts ../data/SEQ_SNPs_counts_0_Target_Erika_Chile_For_Bio.csv --out read_freq.npz > 01_read_based_frequencies.log

python ../scripts/02_diversity_indices.py --geno ../data/SNP_Filtrados_plants.xlsx \
    --meta ../data/STR_maiz_ConMetadato_80per_structure_only.xlsx --reads read_freq.npz > 02_diversity_indices.log

python ../scripts/03_fst_pool_weighted.py --geno ../data/SNP_Filtrados_plants.xlsx \
    --meta ../data/STR_maiz_ConMetadato_80per_structure_only.xlsx \
    --clusters ../data/cluster_maiz_ConMetadato_80per.xlsx --reads read_freq.npz > 03_fst_pool_weighted.log

python ../scripts/04_amova.py --dist ../data/PCOA_1157acc_MR_BIOR.xlsx --geno ../data/SNP_Filtrados_plants.xlsx \
    --meta ../data/STR_maiz_ConMetadato_80per_structure_only.xlsx \
    --perm-global 999 --perm-pair 9999 --seed 1 --out pairwise_PhiST_races.csv > 04_amova.log

python ../scripts/05_structure_consistency.py --runs ../data/structure/output_STR_K2/*_f > 05_structure_consistency_K2.log
python ../scripts/05_structure_consistency.py --runs ../data/structure/output_STR_K6/*_f > 05_structure_consistency_K6.log

python ../scripts/06_knn_inference.py --dist ../data/PCOA_1157acc_MR_BIOR.xlsx --geno ../data/SNP_Filtrados_plants.xlsx \
    --meta ../data/STR_maiz_ConMetadato_80per_structure_only.xlsx > 06_knn_inference.log

python ../scripts/07_ld_decay.py --geno ../data/SNP_Filtrados_plants.xlsx \
    --meta ../data/STR_maiz_ConMetadato_80per_structure_only.xlsx --reads read_freq.npz --min-maf 0.10 > 07_ld_decay.log

python ../scripts/08_distance_pcoa_sensitivity.py --geno ../data/SNP_Filtrados_plants.xlsx --reads read_freq.npz \
    --dist ../data/PCOA_1157acc_MR_BIOR.xlsx > 08_distance_pcoa_sensitivity.log
```

Reading the 96 MB Excel matrix takes a few minutes per script. Script 04 runs 9,999 permutations for each of the
231 pairwise race comparisons and is the slowest step.

---

## Scripts and the results they produce

| Script | Analysis | Inputs | Output files | Article |
|---|---|---|---|---|
| `01_read_based_frequencies.py` | Within-pool alternative-allele frequencies from read counts; concordance with genotype calls | SNP matrix, read counts | `read_freq.npz` | Table S11 (a) |
| `02_diversity_indices.py` | MAF, gene diversity, observed heterozygosity, Shannon index; alleles shared/exclusive between landraces and breeding materials; race-private alleles | SNP matrix, metadata, `read_freq.npz` (optional) | — | Table 1, Table S11 (b) |
| `03_fst_pool_weighted.py` | Pairwise pool-size-weighted FST among STRUCTURE groups and hierarchical clusters; partition of total gene diversity | SNP matrix, metadata, clusters, `read_freq.npz` (optional) | — | Table 2, Table S6, Table S11 (c) |
| `04_amova.py` | AMOVA for K = 2, K = 6 and the 22 standard races; pairwise PhiST between races (permutation tests, Benjamini–Hochberg FDR) | BIO-R matrix, SNP matrix, metadata | `pairwise_PhiST_races.csv` | Section 2.6, Table S7 |
| `05_structure_consistency.py` | Agreement among independent STRUCTURE runs (H′, SD of Q, stable assignments) | STRUCTURE `*_f` files | — | Table S3 |
| `06_knn_inference.py` | k-NN (k = 3, 5, 10) inference of population type and race for unclassified accessions | BIO-R matrix, SNP matrix, metadata | `knn_population_type.csv`, `knn_race.csv` | Tables S8, S9 |
| `07_ld_decay.py` | LD (r²) decay up to 1 Mb among chromosome-anchored SNPs in landraces | SNP matrix, metadata, `read_freq.npz` (optional) | `ld_pairs_<coding>_maf<x>.csv` | Table S2, Figure 2, Table S11 (e) |
| `08_distance_pcoa_sensitivity.py` | Check of the BIO-R matrix; distances and PCoA with dosage vs read-based frequencies | SNP matrix, `read_freq.npz`, BIO-R matrix | — | Table S11 (d) |
| `core_hunter/core_estratificada_corehunter.R` | Stratified core collection (Core Hunter 3), four objective configurations | Core Hunter genotypes, strata | `CORE_STRAT_<config>.csv`, `core_estratificada_log.txt` | [Section / Table] |

Notes
* Allele frequencies use the operational dosage coding 0 / 0.5 / 1 (alternative allele) unless read-based within-pool
  frequencies are requested with `--reads`.
* In Table S9, the race of LP1186 was set to Zapalote Chico by expert criterion (passport origin); script 06 reports the
  k-NN result.

---

## Core-collection selection (Core Hunter)

`scripts/core_hunter/core_estratificada_corehunter.R` selects a stratified core collection of 231 accessions (20 %)
with Core Hunter 3 (R package `corehunter`, which requires Java 8 or later).

**Inputs** (read from the working directory):
* `genotypes_biparental_20956SNP_1157acc.csv.gz` — 1,157 accessions (rows; first column = accession ID) x 20,956 SNPs,
  coded 0/1/2 (1 = heterozygous), missing = `NA`.
* `strata_core_1157acc.csv` — columns `ID`, `STRATUM` (LANDRACE or BRED MATERIAL) and `FORCED_RACE_REP` (TRUE/FALSE).

**Design**
* Two strata: LANDRACE (978) and BRED MATERIAL (179). These counts include the 47 accessions of unknown population
  type, assigned by k-NN (script 06): 40 as landraces and 7 as breeding materials.
* Proportional allocation: 195 landraces + 36 breeding materials = 231 accessions.
* One representative of each of the 22 standard races is forced into the core (`always.selected`) and counts
  within the 195 landraces.
* No accession is excluded for missing data.
* Each stratum is sampled separately with the same objectives. The search stops after 2,000 steps without
  improvement (`impr.steps = 2000`); `set.seed(20260928)` is called before each run.

**Objective configurations** (one core per configuration; CV = allele coverage, AN = accession-to-nearest-entry
distance, EN = entry-to-nearest-entry distance, both with Modified Rogers distance computed by Core Hunter):

| Configuration | Objectives and weights |
|---|---|
| `CV07_EN03` | CV 0.7, EN-MR 0.3 |
| `CV03_AN07` | CV 0.3, AN-MR 0.7 |
| `CV02_AN04_EN04` | CV 0.2, AN-MR 0.4, EN-MR 0.4 |
| `CV05_AN05` | CV 0.5, AN-MR 0.5 |

[State which configuration was retained in the article.]

**Run**
```bash
cd data/core_hunter
Rscript ../../scripts/core_hunter/core_estratificada_corehunter.R
```
Outputs: `CORE_STRAT_<configuration>.csv` (columns `ID`, `STRATUM`, `SELECTED_BY` = FORCED_RACE_REP or CORE HUNTER)
and `core_estratificada_log.txt` (date, R and `corehunter` versions, seed, stopping rule). Move them to `results/`.

---

## Data access and size

[If the large files are deposited elsewhere, e.g. Zenodo, state here: "The SNP matrix and read-count file are available
at https://doi.org/10.5281/zenodo.XXXXXXX. Download them into `data/` before running the scripts."]

---

## License

Code: [e.g. MIT License]. Data: [e.g. CC BY 4.0]. See `LICENSE`.

## Citation

If you use these data or scripts, please cite the article above. [Add the repository DOI if archived on Zenodo.]

## Contact

[Name, institution, e-mail]
