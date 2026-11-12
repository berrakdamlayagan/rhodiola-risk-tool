# -*- coding: utf-8 -*-
# =====================================================================
# Rhodiola rosea dietary supplements — elemental impurity health-risk
# Monte Carlo REPRODUCIBILITY SCRIPT  (Supplementary material)
# Project THD-2025-5312 | Statistics: B. Mestav
#
# Running this script:
#   (1) prints Table 4 (deterministic LB) and Table 5 (Monte Carlo P95),
#   (2) saves Figure 1 (lead concentrations) and Figure 2 (LB/MB/UB P95
#       hazard index) as PNG files in the current folder.
#
# Requirements : Python 3 + numpy  (matplotlib is needed only for the
#                figures; if it is missing the tables are still printed).
# Run          : python3 rhodiola_monte_carlo.py
#
# Model
#   EDI   = C(mg/kg) * IR(g/day) * 1e-3 / BW          [mg/kg-bw/day]
#   HQ_e  = EDI_e / RfD_e ;  HI = sum over {As,Cd,Ni,Hg,Co,V,Sb,Sn}
#           (Pb and total Cr are EXCLUDED from HI/CR)
#   LADD_As = EDI_As * (ED/LT) ;  CR_As = LADD_As * CSF_As
#           (inorganic-As screening; total As provisionally = 100% iAs)
# Constants : ED=30 y, LT=70 y, EF=365, AT=EF*ED, CF=1e-3,
#             BW ~ Uniform(50,90) kg, N=10000
# RNG       : numpy legacy MT19937 via np.random.seed(SEED), reset once
#             per product, then np.random.uniform / np.random.choice
# SEED      : 20260805
# Censoring : LB  <MLoD->0,        <MLoQ->0
#             MB  <MLoD->MLoD/2,   <MLoQ->MLoQ/2
#             UB  <MLoD->MLoD,     <MLoQ->MLoQ
# =====================================================================
import numpy as np

EL = ['As','Cd','Pb','Cr','Ni','Hg','Co','V','Sb','Sn']

# ---- Table-4 censored concentrations (mg/kg).  'D'=<MLoD, 'Q'=<MLoQ ----
T4 = {
 'S1' :[240.7909,'D',2.6983,43.3909,1.0378,'D','D',37.7431,0.2185,0.4697],
 'S2' :[283.9054,3.0693,7.4037,53.1243,4.7635,'D',11.4081,40.165,0.275,1.0231],
 'S3' :[200.7309,'Q',2.3367,92.4474,1.5282,'D',9.4352,24.629,0.2215,0.6977],
 'S4' :['Q','D',1.5386,25.3708,0.7084,'D',0.5424,2.1253,0.2214,'Q'],
 'S5' :['D','D',1.7837,'Q','Q','D','Q',1.1341,0.2202,'Q'],
 'S6' :['D','D',2.3413,'D',1.5424,'D',0.4549,0.921,0.1775,'Q'],
 'S7' :['Q','D',2.3862,16.4394,14.8669,'D',1.2866,3.1449,0.1869,'Q'],
 'S8' :['Q','Q',2.1049,18.9112,16.4426,'D',1.3335,3.1078,0.2094,0.1433],
 'S9' :['D','D',1.9894,16.8279,21.1035,'D',1.0552,2.7149,0.1649,0.4726],
 'S10':['D','D',1.3051,4.1127,3.8824,'D',0.3729,1.0748,0.1097,'Q'],
 'S11':['Q','Q',1.6274,4.1896,5.861,'D',0.5608,1.3855,0.099,'D'],
 'S12':['Q','D',2.0189,82.9713,5.6287,'D',0.6583,1.6129,0.0878,'Q'],
 'S13':['Q','D',2.8387,'D',2.0824,'D',0.3178,0.1754,0.0767,'Q'],
 'S14':['D','D',1.5175,17.9459,1.0227,'D','D','D',0.044,'D'],
 'S15':['D','D',1.8659,4.2417,50.6617,'D',0.7702,'D',0.0542,'D'],
 'S16':['D','D',1.0437,'D',1.7908,'D',0.3845,'D','Q','D'],
 'S17':['D','D',5.5885,2.152,0.8564,'D',0.3074,'D','D','D'],
 'S18':['D','D',3.1078,'Q',9.2906,'D','Q','D','Q','D'],
}
PROD = {s:p for s,p in zip(T4, sum([[c]*3 for c in 'ABCDEF'],[]))}  # S1-3=A ... S16-18=F
IR   = {'A':0.20,'B':0.25,'C':0.35,'D':0.25,'E':0.25,'F':0.25}       # g/day

# ---- Representative matrix MLoD/MLoQ (mg/kg) ----
MLoD = {'As':0.550,'Cd':0.689,'Pb':0.248,'Cr':0.159,'Ni':0.136,
        'Hg':1.470,'Co':0.0834,'V':0.0461,'Sb':0.0121,'Sn':0.0351}
MLoQ = {e: MLoD[e]*3.33 for e in EL}
MLoQ.update({'As':1.833,'Cd':2.293,'Pb':0.823,'Hg':4.894})

RfD = {'As':6e-5,'Cd':1e-3,'Ni':2e-2,'Hg':3e-4,'Co':3e-4,'V':5e-3,'Sb':4e-4,'Sn':0.6}
CSF_As = 32.0
EF,ED,LT,CF,N,SEED = 365,30,70,1e-3,10000,20260805
AT = EF*ED

def value(e, cell, scen):
    if not isinstance(cell,str): return cell
    if cell=='D': return {'LB':0.0,'MB':MLoD[e]/2,'UB':MLoD[e]}[scen]
    if cell=='Q': return {'LB':0.0,'MB':MLoQ[e]/2,'UB':MLoQ[e]}[scen]
    raise ValueError(cell)

def hi_cr(sample, ir, bw):
    hi = sum((sample[e]*ir*CF/bw)/RfD[e] for e in RfD)
    cr = (sample['As']*ir*CF/bw)*(ED/LT)*CSF_As
    return hi, cr

# ---------------- Table 4: deterministic LB (BW = 70 kg) ----------------
print("Table 4 - deterministic lower-bound (BW = 70 kg)")
print(f"{'Sample':7}{'Product':8}{'HI':>10}{'CR_As':>12}")
_t4 = ["Sample,Product,HI_lower_bound,CR_As_lower_bound"]
for s in T4:
    smp = {EL[j]: value(EL[j], T4[s][j], 'LB') for j in range(10)}
    hi, cr = hi_cr(smp, IR[PROD[s]], 70.0)
    print(f"{s:7}Sup_{PROD[s]:<4}{hi:10.4f}{cr:12.2e}")
    _t4.append(f"{s},Sup_{PROD[s]},{hi:.4f},{cr:.3e}")
open("Table4.csv","w",encoding="utf-8").write("\n".join(_t4))

# ---------------- Table 5: Monte Carlo P95 ----------------
print("\nTable 5 - Monte Carlo P95 (seed=%d, N=%d)" % (SEED, N))
print(f"{'Product':8}{'Scen':5}{'HI_P95':>10}{'CR_As_P95':>12}{'CR>1e-4 %':>11}")
P95_HI = {}   # (scen, product) -> HI P95  (collected for Figure 2)
_t5 = ["Product,Scenario,HI_P95,CR_As_P95,CR_As_over_1e-4_percent"]
for scen in ['LB','MB','UB']:
    for p in 'ABCDEF':
        samps = [{EL[j]: value(EL[j], T4[s][j], scen) for j in range(10)}
                 for s in T4 if PROD[s]==p]
        np.random.seed(SEED)                       # reset once per product
        bw  = np.random.uniform(50, 90, N)
        idx = np.random.choice([0,1,2], N)
        HI  = np.array([hi_cr(samps[idx[i]], IR[p], bw[i])[0] for i in range(N)])
        CR  = np.array([hi_cr(samps[idx[i]], IR[p], bw[i])[1] for i in range(N)])
        P95_HI[(scen,p)] = np.percentile(HI,95)
        print(f"Sup_{p:<4}{scen:5}{np.percentile(HI,95):10.3f}"
              f"{np.percentile(CR,95):12.2e}{100*(CR>1e-4).mean():11.1f}")
        _t5.append(f"Sup_{p},{scen},{np.percentile(HI,95):.3f},{np.percentile(CR,95):.3e},{100*(CR>1e-4).mean():.1f}")
open("Table5.csv","w",encoding="utf-8").write("\n".join(_t5))

# ===================== FIGURES =====================
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("\n[Note] matplotlib not installed - tables printed, figures skipped.")
    print("       Install with:  pip install matplotlib")
    raise SystemExit

plt.rcParams.update({"font.size":10, "savefig.dpi":600,
                     "axes.edgecolor":"#555", "axes.linewidth":0.8})

# ---- Figure 1: lead concentrations in the 18 packages ----
labs = [f"{c}{k+1}" for c in "ABCDEF" for k in range(3)]
pb   = [T4[f"S{i+1}"][2] for i in range(18)]            # Pb = index 2 (all quantified)
exc  = [v > 3.0 for v in pb]
fig, ax = plt.subplots(figsize=(9.0, 3.3))
ax.bar(range(18), pb, width=0.74, zorder=3, edgecolor="white", linewidth=0.6,
       color=["#c0702f" if e else "#5b6b7a" for e in exc])
ax.axhline(3.0, color="#a32d2d", ls="--", lw=1.2, zorder=4)
ax.text(17.4, 3.12, "3.0 mg/kg comparison level", color="#a32d2d",
        fontsize=9, ha="right", va="bottom")
for i, v in enumerate(pb):
    if exc[i]:
        ax.text(i, v+0.12, f"{v:.2f}", ha="center", va="bottom",
                fontsize=8.5, fontweight="bold")
ax.set_xticks(range(18)); ax.set_xticklabels(labs, fontsize=8.5)
ax.set_ylabel("Lead concentration (mg/kg)"); ax.set_ylim(0, 8.2)
ax.grid(axis="y", color="#e6e6e6", lw=0.6, zorder=0)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig("Figure1_Pb.png", bbox_inches="tight"); plt.close()

# ---- Figure 2: LB/MB/UB P95 hazard index by product ----
SUP = list("ABCDEF"); x = np.arange(6)
style = {'LB':("#2e5a7a","-","o","Lower-bound (LB)"),
         'MB':("#7a6a55","--","s","Middle-bound (MB)"),
         'UB':("#c0702f",":","o","Upper-bound (UB)")}
fig, ax = plt.subplots(figsize=(9.0, 3.6))
for scen in ['LB','MB','UB']:
    col, ls, mk, lab = style[scen]
    y = [P95_HI[(scen,p)] for p in SUP]
    ax.plot(x, y, ls, color=col, marker=mk, ms=5, lw=1.6, label=lab, zorder=3)
ax.axhline(1, color="#333", ls="-.", lw=1.1, label="HI = 1 screening level", zorder=2)
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([f"Sup_{p}" for p in SUP])
ax.set_ylabel("Modeled P95 hazard index (log scale)")
ax.grid(axis="y", which="both", color="#ededed", lw=0.6, zorder=0)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.legend(fontsize=8.5, ncol=2, frameon=False, loc="upper right")
plt.tight_layout(); plt.savefig("Figure2_HI_scenarios.png", bbox_inches="tight"); plt.close()

print("\nSaved files:  Table4.csv  Table5.csv  Figure1_Pb.png  Figure2_HI_scenarios.png")
