# Rhodiola rosea — Decision Support Model (DSM) for Elemental-Impurity Screening

Software and reproducible analysis accompanying a manuscript on a **Decision
Support Model (DSM)** that couples deterministic hazard indexing with Monte
Carlo simulation to screen elemental (heavy-metal) impurities in commercial
*Rhodiola rosea* dietary supplements. The six analysed lots serve as a
**case study** demonstrating the model's behaviour under high-exposure
(worst-case) conditions.

**Project:** THD-2025-5312, Scientific Research Projects (BAP), Çanakkale Onsekiz Mart University
**Team:** Prof. Dr. Özlem Tonguç Yayıntaş (Principal Investigator) ·  Res. Asst. Dr. Berrak Damla Yağan  · Assoc. Prof. Dr. Burcu Mestav

---

## Live pages

Once GitHub Pages is enabled (see below), two pages are served:

| Page | Address | What it is |
|------|---------|------------|
| **DSM software** | `https://berrakdamlayagan.github.io/rhodiola-risk-tool/` | The interactive Decision Support Model (main link cited in the manuscript) |
| **Reproducible analysis** | `https://berrakdamlayagan.github.io/rhodiola-risk-tool/rhodiola_monte_carlo.html` | Rendered report: code + tables + figures in one page (validation evidence) |

👉 **[Open the DSM software](https://berrakdamlayagan.github.io/rhodiola-risk-tool/)**
*(this is the link to paste into the manuscript's `[Link to the DSM software]` placeholder)*

The DSM runs entirely in the browser; no data are uploaded to any server. It
computes estimated daily intake (EDI), hazard quotient (HQ), cumulative hazard
index (HI) and the arsenic-based cancer screening index (CR) by Monte Carlo
simulation, with single-product and multi-product comparison modes and
lower-/middle-/upper-bound substitution scenarios for left-censored values.

## Contents

| File | Description |
|------|-------------|
| `index.html` | The interactive Decision Support Model (served at the root link above) |
| `rhodiola_monte_carlo.html` | Rendered reproducible analysis (code + narrative + tables + figures) |
| `rhodiola_monte_carlo.qmd` | Quarto source of the analysis above |
| `rhodiola_monte_carlo.py` | Monte Carlo reproducibility script (Python); reproduces the manuscript tables and figures |
| `plot from excel.py` | Reads the Excel workbook and draws the publication figures |
| `Rhodiola_MonteCarlo_Analiz.xlsx` | Full Monte Carlo workbook (data, deterministic, brand-based MC, live engine, charts) |
| `Figure1_Pb_600dpi.png`, `Figure2_HI_scenarios_600dpi.png` | High-resolution manuscript figures |

## How to use the DSM

Open the live link above in any modern browser (Chrome, Edge, Firefox, Safari).
Use the **EN / TR** button for English or Turkish. Load the embedded example
(six supplement lots) or enter your own concentrations, then run the simulation.

## How to reproduce the analysis

Two equivalent routes:

- **Pre-rendered:** open `rhodiola_monte_carlo.html` — nothing to install.
- **From source** (Python 3 with `numpy`; `matplotlib` for figures, `openpyxl` for the Excel reader):

```bash
pip install numpy matplotlib openpyxl
python rhodiola_monte_carlo.py
```

This prints Table 4 (deterministic lower-bound) and Table 5 (Monte Carlo P95)
and saves the figures. Results are reproducible (fixed random seed 20260805).
The Quarto version re-renders the HTML/PDF with `quarto render rhodiola_monte_carlo.qmd`.

## Model (summary)

- `EDI = C(mg/kg) × IR(g/day) × 1e-3 / BW`; `HQ = EDI / RfD`; `HI = Σ HQ`
  over As, Cd, Ni, Hg, Co, V, Sb, Sn (Pb and total Cr excluded — no valid oral
  reference value / speciation needed).
- Cancer screening (inorganic arsenic only): `CR_As = EDI_As × (ED/LT) × CSF_As`
  (lifetime averaging; ED = 30 y, LT = 70 y, CSF_As = 32).
- Monte Carlo: 10,000 iterations, body weight ~ Uniform(50, 90) kg, resampling
  of same-lot packages; left-censored values handled under lower- (0),
  middle- (LoQ/2) and upper-bound (LoQ) substitution.

**Conservative-screening design.** Total arsenic is deliberately treated as
100 % inorganic. This is a built-in **worst-case / conservative screening**
assumption, not an analytical limitation: it lets the DSM flag the highest
plausible exposure before any species-specific refinement. Accordingly the
arsenic-based CR is a **screening indicator, not a validated lifetime cancer
risk** — the doses can exceed the range where the 2025 U.S. EPA linear slope
factor applies (~0.2 µg/kg-day). Outputs are preliminary screening results for
the single lots examined and must not be generalised to the market.

## Notes

- Product identities are anonymised as **Sup_A–Sup_F** throughout.
- Toxicity reference values (RfD, CSF) are current as of 2026 and may be
  revised by regulatory agencies; verify against current sources before use.

## Enabling the live pages (GitHub Pages)

Repository **Settings → Pages → Build and deployment → Source: Deploy from a
branch → Branch: `main` / folder: `/ (root)` → Save**. After a minute the DSM
is served at `https://berrakdamlayagan.github.io/rhodiola-risk-tool/` and the analysis at
`https://berrakdamlayagan.github.io/rhodiola-risk-tool/rhodiola_monte_carlo.html`.

## Citation

> Yayıntaş Ö. T., Yağan B. D., Mestav B. (2026). *Rhodiola rosea elemental
> impurity Decision Support Model (DSM) and Monte Carlo analysis* [Software].
> Çanakkale Onsekiz Mart University, project THD-2025-5312. Available at:
> https://github.com/berrakdamlayagan/rhodiola-risk-tool

*(If you archive a release on Zenodo, add the DOI here.)*

## License

Released for academic use. Suggested: MIT License for the code and
CC-BY 4.0 for the data and figures (add a `LICENSE` file to confirm).
