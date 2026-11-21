# Rhodiola rosea — Decision-Support Screening of Elemental Impurities

Interactive tool and reproducible analysis for a pilot **decision-support model
(DSM)** that couples deterministic hazard indexing with Monte Carlo simulation
to screen elemental (heavy-metal) impurities in commercial *Rhodiola rosea*
dietary supplements. Six analysed lots (18 packages) serve as a case study.

**Project:** THD-2025-5312, Scientific Research Projects (BAP), Çanakkale Onsekiz Mart University
**Team:** Prof. Dr. Özlem Tonguç Yayıntaş (Principal Investigator) · Res. Asst. Dr. Berrak Damla Yağan · Assoc. Prof. Dr. Burcu Mestav

---

## Live pages and archive

| Resource | Address |
|----------|---------|
| **Interactive tool (DSM)** | https://berrakdamlayagan.github.io/rhodiola-risk-tool/ |
| **Reproducible analysis (report)** | https://berrakdamlayagan.github.io/rhodiola-risk-tool/rhodiola_monte_carlo.html |
| **Archived release (DOI)** | https://doi.org/10.5281/zenodo.24167368 |

The tool runs entirely in the browser; no data are uploaded to any server. It
computes the estimated daily intake (EDI), hazard quotient (HQ), cumulative
hazard index (HI) and the arsenic-based cancer screening index (CR) by Monte
Carlo simulation, with single-product and multi-product comparison modes and
lower-, middle- and upper-bound substitution scenarios for left-censored values.
Every result also reports the hazard index **with arsenic excluded** as a
sensitivity scenario.

## Contents

| File | Description |
|------|-------------|
| `index.html` | The interactive decision-support tool (English) |
| `rhodiola_monte_carlo.html` | Rendered reproducible analysis (code, tables and figures) |
| `rhodiola_monte_carlo.py` | Reproducibility script; writes the tables, CSVs and figures |
| `rhodiola_monte_carlo.qmd` | Quarto source of the analysis above |
| `Rhodiola_MonteCarlo_Analysis_EN.xlsx` | Workbook in English, aligned with the manuscript: constants, concentrations, live deterministic calculation, Monte Carlo P95, and the arsenic-exclusion column |
| `input_template.csv` | Column template and worked example for the comparison mode |
| `Figure1_Pb.png`, `Figure2_HI_scenarios.png`, `Figure3_As_excluded_sensitivity.png` | High-resolution figures |
| `requirements.txt`, `CITATION.cff`, `LICENSE` | Dependencies, citation metadata and licence |

## Using the tool

Open the live link in any modern browser. Load the embedded six-product example
or enter your own concentrations, choose the daily intake and body-weight range,
and run the simulation. The interface is in English.

## Reproducing the analysis

Requires Python 3 with the packages in `requirements.txt`:

```bash
pip install -r requirements.txt
python rhodiola_monte_carlo.py
```

This prints Table 4 (deterministic lower bound) and Table 5 (Monte Carlo P95),
writes `Table4_deterministic.csv` and `Table5_montecarlo.csv`, and saves the
three figures. Results are reproducible (fixed random seed 20260805). The Quarto
version renders to HTML or PDF with `quarto render rhodiola_monte_carlo.qmd`.

## Input format (comparison mode)

`input_template.csv` has one row per package. The first two columns are
`Product` and `IR` (daily product mass, g/day); the remaining columns are the
element concentrations in mg/kg (`As, Cd, Pb, Cr, Ni, Hg, Co, V, Sb, Sn`). A
value below the reporting limit is entered as `0` for the lower-bound scenario.
Rows that share a `Product` name are treated as packages of the same product.

## Method (summary)

- `EDI = C(mg/kg) × IR(g/day) × 1e-3 / BW`; `HQ = EDI / RfD`; `HI = Σ HQ`
  over As, Cd, Ni, Hg, Co, V, Sb, Sn (Pb and total Cr excluded).
- Cancer screening (inorganic arsenic only): `CR_As = EDI_As × (ED/LT) × CSF_As`
  (lifetime averaging; ED = 30 y, LT = 70 y, CSF_As = 32).
- Monte Carlo: 10,000 iterations, body weight ~ Uniform(50, 90) kg, resampling
  of the within-lot packages; left-censored values handled under lower- (0),
  middle- (LoQ/2) and upper-bound (LoQ) substitution.
- Sensitivity: the hazard index is also computed with arsenic removed from the
  summation, which isolates the influence of the reported arsenic value.

**Important.** Total arsenic is treated as 100% inorganic, a deliberately
conservative screening assumption. The arsenic-based CR is a screening
indicator, not a validated lifetime cancer risk, because the doses can exceed
the range where the 2025 U.S. EPA linear slope factor applies (~0.2 µg/kg-day).
Outputs are preliminary screening results for the single lots examined and must
not be generalised to the market. The reported arsenic signal in Sup_A requires
interference-controlled analytical confirmation before any toxicological
interpretation.

## Notes

- Product identities are anonymised as **Sup_A–Sup_F** throughout.
- Toxicity reference values (RfD, CSF) are current as of 2026 and may be revised
  by regulatory agencies; verify against current sources before use.

## Enabling the live pages (GitHub Pages)

Repository **Settings → Pages → Build and deployment → Source: Deploy from a
branch → Branch: `main` / folder: `/ (root)` → Save**. After about a minute the
tool is served at the address above and the analysis at `.../rhodiola_monte_carlo.html`.

## Citation

> Tonguç Yayıntaş Ö., Yağan B. D., Mestav B. (2026). *Rhodiola rosea
> elemental-impurity decision-support screening tool and Monte Carlo analysis*
> [Software]. Çanakkale Onsekiz Mart University, project THD-2025-5312.
> https://doi.org/10.5281/zenodo.24167368

## License

Code is released under the MIT License (see `LICENSE`). Data and figures are
released under CC-BY 4.0.
