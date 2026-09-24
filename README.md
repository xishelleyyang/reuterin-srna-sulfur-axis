# reuterin-srna-sulfur-axis

Reanalysis code for: **"Public transcriptome integration reveals an
electrophile-detoxification and ISC-biased expression program in the bacterial
response to reuterin and nominates candidate sRNA regulators."**

This repository integrates four public resources -- a reuterin-exposure
microarray (GSE19760), two cross-stress RNA-seq comparators (GSE126176,
GSE165009), the Hfq CLASH sRNA-target interactome (Iosub et al. 2020, eLife
54655), and 13 RefSeq *E. coli*/*Salmonella* genomes -- to characterize the
transcriptional response to reuterin and nominate candidate post-transcriptional
regulators. No new sequencing data is generated; everything here is a
from-scratch statistical reanalysis of already-public data.

## Repository layout

```
download_inputs.py          fetches every raw input from GEO / eLife / NCBI RefSeq
run_all.py                  orchestrates the 4 analysis modules + tables + figures
config/input_manifest.json  source URL + SHA-256 for every raw input file
src/
  recompute_statistics.R      Module 1+2: microarray + RNA-seq DE, gene-set roast tests
  build_tables.py              Module 3 (CLASH network) + assembles Tables S1-S7
  conservation_analysis.py     Module 4: ryhB/iscS 5'UTR conservation across 13 genomes
  typhimurium_duplex_mfe.py    RyhB-iscS duplex free energy (MG1655 vs S. Typhimurium LT2)
  mg1655_queries.py            shared MG1655 query-sequence coordinates
  make_figures.py              Figures 1-4
  make_fig5_conservation.py    Figure 5
tests/
  test_gene_sets.py            regression guard on the curated gene-set definitions
  test_pipeline_smoke.py       output-file + headline-number smoke test (used by `run_all.py --stage validate`)
```

## Setup

**Python** (tested: 3.11, package versions in `requirements.txt`):
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```
`ViennaRNA` is the Python binding package name for the `RNA` module used by the
conservation and duplex-MFE scripts; it ships its own compiled extension (no
system ViennaRNA install required). `conservation_analysis.py` additionally
needs NCBI BLAST+ (`blastn`, `makeblastdb`) on `PATH` -- tested against 2.16.0+.

**R** (tested: 4.4.2):
```r
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"))   # pulls in statmod, locfit automatically
```
Tested versions: limma 3.62.1, edgeR 4.4.0, statmod 1.5.2, locfit 1.5-9.12.

## Quick start

```bash
python download_inputs.py --dest data          # ~35 MB, verified against GEO/eLife/RefSeq
python run_all.py --data-dir data --out-dir results --stage all --pilot   # fast sanity pass (~1-2 min)
python run_all.py --data-dir data --out-dir results --stage all           # full run (NROT=9999 rotations)
python tests/test_pipeline_smoke.py results/formal
```
`--pilot` runs the gene-set rotation tests at 199 rotations instead of 9999 and
writes to `results/pilot/` rather than `results/formal/`; pilot output is
clearly marked (`PILOT_DO_NOT_SUBMIT.txt`) and is for interface/dependency
sanity-checking only -- every q-value quoted in the manuscript requires the
formal (non-pilot) run.

## What each stage does

| stage | script(s) | inputs | key outputs |
|---|---|---|---|
| `transcriptomes` | `recompute_statistics.R DATA OUT [--pilot]` | GPL7445 SOFT + 5 GPR files; GSE126176 counts; 4 GSE165009 untreated samples; MG1655 feature table | `DE_{Reuterin,HOCl,Ferrate,Hpx}.csv`, `genesets_primary.csv` + sensitivity variants, `array_rotation_eligibility.csv`, `batch_sensitivity.csv` |
| `conservation` | `conservation_analysis.py` then `typhimurium_duplex_mfe.py` | MG1655 genome + 13 target genomes | `Table_S3_conservation_v3.csv`, `alignments_v3.json`, `Table_S3b_duplex_MFE.csv` |
| `tables` | `build_tables.py --data --stats --compgen --out` | outputs of the two stages above + CLASH xlsx | `Table_S1`-`Table_S7_..._v4.csv`, `result_summary_v4.json` |
| `figures` | `make_figures.py` then `make_fig5_conservation.py` | `stats`/`compgen`/`final` (Table S1) directories | `Fig1_..._v4.{png,svg}` through `Fig5_..._v4.{png,svg}` |
| `publish` | (inside `run_all.py`) | the two steps above | the same 7 tables / 5 figures, copied to the `_v5`-suffixed filenames actually cited in the manuscript text and captions |
| `validate` | `tests/test_pipeline_smoke.py` | `results/<formal\|pilot>/` | pass/fail against manuscript-reported numbers |

"Module 3" (the CLASH interactome Fisher tests) is not a separate CLI stage: it
is one function inside `build_tables.py` (`network_tables()`) alongside the
table-assembly logic it feeds, and that script is reused unmodified from the
verified statistics engine rather than split apart for a cosmetic 1:1 stage
mapping.

Every one of `recompute_statistics.R`, `build_tables.py`,
`conservation_analysis.py`, `make_figures.py`, and `make_fig5_conservation.py`
takes explicit `--data`/`--out`-style arguments; nothing hardcodes an absolute
path, and the whole pipeline runs correctly from any working directory
(verified by running it from a directory with none of this session's history).

## Duplex free energy for S. Typhimurium LT2

The manuscript's Discussion compares the predicted RyhB-iscS duplex MFE for
MG1655 (-54.7 kcal/mol) against S. Typhimurium LT2 (-76.0 kcal/mol).
`conservation_analysis.py` computes the MG1655 side directly (it is the BLAST
query); the LT2 side is obtained by `typhimurium_duplex_mfe.py`, which applies
the identical `RNA.duplexfold` call to LT2's own genomic ryhB-1 and iscS 5'UTR
sequences (recovered by stripping alignment gap characters from the LT2 entry
that `conservation_analysis.py` already produces in `alignments_v3.json` -- no
new alignment or method, just the missing second half of a comparison whose
first half was already computed). `alignments_v3.json` holds the same
alignment for all 13 genomes, so `--target-genome` can point at any of them;
only MG1655 and LT2 are reported because those are the only two values stated
in the manuscript text.

## Tests

```bash
python tests/test_gene_sets.py          # no data required, always runnable
python tests/test_pipeline_smoke.py results/formal   # after a full run_all.py pass
```
`test_gene_sets.py` guards against the historical class of bug this audit
found and fixed upstream of this repository: an earlier, now-retired script
used the invalid gene symbol `gorA` (canonical symbol is `gor`, b3500) and did
not exclude the Hpx-strain-deleted genes (`katE`, `katG`, `ahpC`) from
gene-set scoring for that contrast. Both checks parse the gene-set
definitions directly out of the shipped `src/` files, not a second hand-typed
copy, so they fail if the shipped definitions themselves ever regress.

## Data availability

No raw data is redistributed in this repository. `download_inputs.py` fetches
everything from its original public source (GEO, eLife, NCBI RefSeq) and
verifies each file's SHA-256 before use; see `config/input_manifest.json` for
the exact accessions, URLs, and hashes. Derived output tables (Tables S1-S7)
accompany the manuscript as supplementary files.

## License

MIT (see `LICENSE`). Analyzed data remain subject to their original sources'
terms (GEO submitter deposits; eLife CC-BY; NCBI RefSeq).
