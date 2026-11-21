#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rhodiola rosea elemental-impurity decision-support screening
Reproducible deterministic and Monte Carlo analysis.

Reproduces:
  - Table 4  : deterministic lower-bound hazard index (HI), HI with arsenic
               excluded, and the arsenic-based cancer screening index (CR_As).
  - Table 5  : Monte Carlo 95th-percentile outputs by product and substitution
               scenario (LB/MB/UB), including HI with arsenic excluded
               (sensitivity scenario) and CR_As with the fraction of iterations
               above 1e-4.
  - Figure 1 : lead concentrations in the 18 packages.
  - Figure 2 : P95 hazard index under the LB/MB/UB substitution scenarios.
  - Figure 3 : P95 hazard index with arsenic included vs. excluded (upper-bound
               scenario), the arsenic-exclusion sensitivity analysis.

Model
  EDI = C(mg/kg) * IR(g/day) * CF / BW       (CF = 1e-3)
  HQ_e = EDI_e / RfD_e ;  HI = sum(HQ) over {As,Cd,Ni,Hg,Co,V,Sb,Sn}
  (Pb and total Cr are excluded; no valid oral RfD / speciation required.)
  Cancer screening (inorganic arsenic only):
  LADD_As = EDI_As * (ED/LT) ;  CR_As = LADD_As * CSF_As
  Left-censored values: '<MLoD' -> 0 / MLoD/2 / MLoD ; '<MLoQ' -> 0 / MLoQ/2 / MLoQ
  under the lower- (LB), middle- (MB) and upper-bound (UB) scenarios.

Requirements: python3, numpy, matplotlib
Run:          python rhodiola_monte_carlo.py
Outputs are written next to this script. Results are reproducible (fixed seed).
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
def OUT(name):
    return os.path.join(HERE, name)

# --------------------------------------------------------------------------- #
# Data (Table 3, censored).  'D' = <MLoD, 'Q' = <MLoQ.
# Order: As, Cd, Pb, Cr, Ni, Hg, Co, V, Sb, Sn
# --------------------------------------------------------------------------- #
EL = ["As", "Cd", "Pb", "Cr", "Ni", "Hg", "Co", "V", "Sb", "Sn"]

T3 = {
    "S1":  [240.7909, "D", 2.6983, 43.3909, 1.0378, "D", "D", 37.7431, 0.2185, 0.4697],
    "S2":  [283.9054, 3.0693, 7.4037, 53.1243, 4.7635, "D", 11.4081, 40.165, 0.275, 1.0231],
    "S3":  [200.7309, "Q", 2.3367, 92.4474, 1.5282, "D", 9.4352, 24.629, 0.2215, 0.6977],
    "S4":  ["Q", "D", 1.5386, 25.3708, 0.7084, "D", 0.5424, 2.1253, 0.2214, "Q"],
    "S5":  ["D", "D", 1.7837, "Q", "Q", "D", "Q", 1.1341, 0.2202, "Q"],
    "S6":  ["D", "D", 2.3413, "D", 1.5424, "D", 0.4549, 0.921, 0.1775, "Q"],
    "S7":  ["Q", "D", 2.3862, 16.4394, 14.8669, "D", 1.2866, 3.1449, 0.1869, "Q"],
    "S8":  ["Q", "Q", 2.1049, 18.9112, 16.4426, "D", 1.3335, 3.1078, 0.2094, 0.1433],
    "S9":  ["D", "D", 1.9894, 16.8279, 21.1035, "D", 1.0552, 2.7149, 0.1649, 0.4726],
    "S10": ["D", "D", 1.3051, 4.1127, 3.8824, "D", 0.3729, 1.0748, 0.1097, "Q"],
    "S11": ["Q", "Q", 1.6274, 4.1896, 5.861, "D", 0.5608, 1.3855, 0.099, "D"],
    "S12": ["Q", "D", 2.0189, 82.9713, 5.6287, "D", 0.6583, 1.6129, 0.0878, "Q"],
    "S13": ["Q", "D", 2.8387, "D", 2.0824, "D", 0.3178, 0.1754, 0.0767, "Q"],
    "S14": ["D", "D", 1.5175, 17.9459, 1.0227, "D", "D", "D", 0.044, "D"],
    "S15": ["D", "D", 1.8659, 4.2417, 50.6617, "D", 0.7702, "D", 0.0542, "D"],
    "S16": ["D", "D", 1.0437, "D", 1.7908, "D", 0.3845, "D", "Q", "D"],
    "S17": ["D", "D", 5.5885, 2.152, 0.8564, "D", 0.3074, "D", "D", "D"],
    "S18": ["D", "D", 3.1078, "Q", 9.2906, "D", "Q", "D", "Q", "D"],
}
# Three packages per product (one lot each)
PROD = {s: p for s, p in zip(T3, sum([[c] * 3 for c in "ABCDEF"], []))}
IR = {"A": 0.20, "B": 0.25, "C": 0.35, "D": 0.25, "E": 0.25, "F": 0.25}  # g/day

# Sample-specific reporting limits (mg/kg)
MLoD = {"As": 0.550, "Cd": 0.689, "Pb": 0.248, "Cr": 0.159, "Ni": 0.136,
        "Hg": 1.470, "Co": 0.0834, "V": 0.0461, "Sb": 0.0121, "Sn": 0.0351}
MLoQ = {e: MLoD[e] * 3.33 for e in EL}
MLoQ.update({"As": 1.833, "Cd": 2.293, "Pb": 0.823, "Hg": 4.894})

# Toxicity reference values (2025 US EPA IRIS for inorganic As)
RfD = {"As": 6e-5, "Cd": 1e-3, "Ni": 2e-2, "Hg": 3e-4,
       "Co": 3e-4, "V": 5e-3, "Sb": 4e-4, "Sn": 0.6}
CSF_As = 32.0

# Constants
EF, ED, LT, CF = 365, 30, 70, 1e-3
N, SEED = 10000, 20260805
HI_ELEMENTS = [e for e in RfD]          # As, Cd, Ni, Hg, Co, V, Sb, Sn

# Colorblind-safe palette (Okabe-Ito)
C_INC, C_EXC = "#0072B2", "#E69F00"
C_LB, C_MB, C_UB = "#0072B2", "#009E73", "#E69F00"


def value(e, cell, scen):
    """Resolve a (possibly censored) cell under a substitution scenario."""
    if not isinstance(cell, str):
        return cell
    if cell == "D":
        return {"LB": 0.0, "MB": MLoD[e] / 2, "UB": MLoD[e]}[scen]
    if cell == "Q":
        return {"LB": 0.0, "MB": MLoQ[e] / 2, "UB": MLoQ[e]}[scen]
    raise ValueError(cell)


def hazard_index(sample, ir, bw, exclude_as=False):
    keys = [e for e in HI_ELEMENTS if not (exclude_as and e == "As")]
    return sum((sample[e] * ir * CF / bw) / RfD[e] for e in keys)


def cr_as(sample, ir, bw):
    edi = sample["As"] * ir * CF / bw
    return edi * (ED / LT) * CSF_As


# --------------------------------------------------------------------------- #
# Table 4 - deterministic lower bound (BW = 70 kg)
# --------------------------------------------------------------------------- #
def table4():
    rows = []
    for s in T3:
        smp = {EL[j]: value(EL[j], T3[s][j], "LB") for j in range(10)}
        hi = hazard_index(smp, IR[PROD[s]], 70.0)
        hi_x = hazard_index(smp, IR[PROD[s]], 70.0, exclude_as=True)
        cr = cr_as(smp, IR[PROD[s]], 70.0)
        rows.append([f"Sup_{PROD[s]}{(int(s[1:]) - 1) % 3 + 1}",
                     round(hi, 4), round(hi_x, 4), f"{cr:.2e}"])
    return rows


# --------------------------------------------------------------------------- #
# Table 5 - Monte Carlo P95 by product and scenario
# --------------------------------------------------------------------------- #
def table5():
    rows = []
    p95_hi = {}
    for scen in ["LB", "MB", "UB"]:
        for p in "ABCDEF":
            samps = [{EL[j]: value(EL[j], T3[s][j], scen) for j in range(10)}
                     for s in T3 if PROD[s] == p]
            np.random.seed(SEED)                    # reset once per product
            bw = np.random.uniform(50, 90, N)
            idx = np.random.choice([0, 1, 2], N)
            HI = np.array([hazard_index(samps[idx[i]], IR[p], bw[i]) for i in range(N)])
            HIx = np.array([hazard_index(samps[idx[i]], IR[p], bw[i], True) for i in range(N)])
            CR = np.array([cr_as(samps[idx[i]], IR[p], bw[i]) for i in range(N)])
            p95_hi[(scen, p)] = (np.percentile(HI, 95), np.percentile(HIx, 95))
            rows.append([f"Sup_{p}", scen,
                         round(np.percentile(HI, 95), 3),
                         round(np.percentile(HIx, 95), 3),
                         f"{np.percentile(CR, 95):.2e}",
                         round(100 * (CR > 1e-4).mean(), 1)])
    return rows, p95_hi


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def style_axes(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", which="major", color="#e3e3e3", lw=0.7, zorder=0)
    ax.grid(axis="y", which="minor", color="none")
    ax.tick_params(which="minor", length=0)


def figure1_pb():
    labs = [f"{c}{k+1}" for c in "ABCDEF" for k in range(3)]
    pb = [T3[f"S{i+1}"][2] for i in range(18)]
    exc = [v > 3.0 for v in pb]
    fig, ax = plt.subplots(figsize=(9, 3.4))
    ax.bar(range(18), pb, width=0.72, zorder=3, edgecolor="white", linewidth=0.6,
           color=[C_EXC if e else C_INC for e in exc])
    ax.axhline(3.0, color="#333333", ls=(0, (6, 3)), lw=1.2, zorder=2)
    for i, v in enumerate(pb):
        if exc[i]:
            ax.text(i, v + 0.12, f"{v:.2f}", ha="center", va="bottom",
                    fontsize=8.3, fontweight="bold", color="#8a5a00")
    ax.set_xticks(range(18)); ax.set_xticklabels(labs, fontsize=8.5)
    ax.set_ylabel("Lead concentration (mg/kg)"); ax.set_ylim(0, 8.2)
    style_axes(ax)
    handles = [Patch(facecolor=C_INC, label="≤ 3.0 mg/kg"),
               Patch(facecolor=C_EXC, label="> 3.0 mg/kg comparator"),
               Line2D([0], [0], color="#333333", ls=(0, (6, 3)), lw=1.2,
                      label="3.0 mg/kg comparator")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=3, frameon=False, fontsize=9)
    plt.tight_layout()
    fig.savefig(OUT("Figure1_Pb.png"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure2_hi(p95_hi):
    sup = list("ABCDEF"); x = np.arange(6)
    style = {"LB": (C_LB, "-", "o", "Lower bound (LB)"),
             "MB": (C_MB, "--", "s", "Middle bound (MB)"),
             "UB": (C_UB, ":", "^", "Upper bound (UB)")}
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for scen in ["LB", "MB", "UB"]:
        col, ls, mk, lab = style[scen]
        ax.plot(x, [p95_hi[(scen, p)][0] for p in sup], ls, color=col, marker=mk,
                ms=6, lw=1.8, label=lab, zorder=3)
    ax.axhline(1, color="#333333", ls=(0, (6, 3)), lw=1.2, label="HI = 1 screening level")
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([f"Sup_{p}" for p in sup])
    ax.set_ylabel("Monte Carlo P95 hazard index (log scale)")
    style_axes(ax)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4,
              frameon=False, fontsize=9)
    plt.tight_layout()
    fig.savefig(OUT("Figure2_HI_scenarios.png"), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def figure3_sensitivity(p95_hi, scen="UB"):
    prods = list("ABCDEF"); x = np.arange(6); w = 0.38
    inc = [p95_hi[(scen, p)][0] for p in prods]
    exc = [p95_hi[(scen, p)][1] for p in prods]
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.bar(x - w / 2, inc, w, color=C_INC, edgecolor="#003f66", linewidth=0.8,
           zorder=3, label="Arsenic included")
    ax.bar(x + w / 2, exc, w, color=C_EXC, edgecolor="#8a5a00", linewidth=0.8,
           hatch="///", zorder=3, label="Arsenic excluded")
    ax.axhline(1.0, color="#333333", ls=(0, (6, 3)), lw=1.3, zorder=2)
    ax.set_yscale("log"); ax.set_ylim(0.004, 40)
    ax.set_xticks(x); ax.set_xticklabels([f"Sup_{p}" for p in prods])
    ax.set_ylabel("Monte Carlo P95 hazard index (log scale)")
    style_axes(ax)
    for xi, v in zip(x - w / 2, inc):
        ax.text(xi, v * 1.12, f"{v:.2f}" if v >= 0.01 else f"{v:.3f}",
                ha="center", va="bottom", fontsize=8.3, color="#003f66")
    for xi, v in zip(x + w / 2, exc):
        ax.text(xi, v * 1.12, f"{v:.3f}", ha="center", va="bottom",
                fontsize=8.3, color="#8a5a00")
    handles = [Patch(facecolor=C_INC, edgecolor="#003f66", label="Arsenic included"),
               Patch(facecolor=C_EXC, edgecolor="#8a5a00", hatch="///", label="Arsenic excluded"),
               Line2D([0], [0], color="#333333", ls=(0, (6, 3)), lw=1.3,
                      label="HI = 1 screening level")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13),
              ncol=3, frameon=False, fontsize=9.3, handlelength=1.8, columnspacing=1.6)
    plt.tight_layout()
    fig.savefig(OUT("Figure3_As_excluded_sensitivity.png"), dpi=600,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.edgecolor": "#444", "figure.facecolor": "white",
                         "axes.facecolor": "white"})

    t4 = table4()
    t5, p95_hi = table5()

    # ---- console ----
    print("\nTable 4 - deterministic lower bound (BW = 70 kg)")
    print(f"{'Sample':8} {'HI':>10} {'HI (As excl.)':>14} {'CR_As':>10}")
    for r in t4:
        print(f"{r[0]:8} {r[1]:>10} {r[2]:>14} {r[3]:>10}")

    print("\nTable 5 - Monte Carlo P95 by product and scenario")
    print(f"{'Product':8} {'Scen':4} {'HI P95':>8} {'HI P95(excl)':>12} "
          f"{'CR_As P95':>10} {'CR>1e-4 %':>9}")
    for r in t5:
        print(f"{r[0]:8} {r[1]:4} {r[2]:>8} {r[3]:>12} {r[4]:>10} {r[5]:>9}")

    # ---- CSV ----
    with open(OUT("Table4_deterministic.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Sample", "HI_LB", "HI_LB_As_excluded", "CR_As_LB"])
        w.writerows(t4)
    with open(OUT("Table5_montecarlo.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Product", "Scenario", "HI_P95", "HI_P95_As_excluded",
                    "CR_As_P95", "CR_As_gt_1e-4_pct"])
        w.writerows(t5)

    # ---- figures ----
    figure1_pb()
    figure2_hi(p95_hi)
    figure3_sensitivity(p95_hi, scen="UB")

    print("\nSaved: Table4_deterministic.csv, Table5_montecarlo.csv,")
    print("       Figure1_Pb.png, Figure2_HI_scenarios.png, "
          "Figure3_As_excluded_sensitivity.png")


if __name__ == "__main__":
    main()
