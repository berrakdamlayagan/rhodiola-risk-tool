# -*- coding: utf-8 -*-
# =====================================================================
# Rhodiola MC — Excel'den okuyup "güzel" (matplotlib) figürleri çizer.
# Reads the workbook and draws the report/manuscript figures.
#
# Gerekir / requires : openpyxl, matplotlib
#   pip install openpyxl matplotlib
# Çalıştır / run     : python grafik_ciz_exceldan.py
#   (Excel dosyası aynı klasörde olmalı.)
# Üretir / outputs   : Figure1_Pb.png, Figure2_LBMBUB.png, Sekil3_det_HI.png,
#                      Sekil4_MC_HI.png, Sekil5_MC_CR.png, Sekil6_katki.png
# =====================================================================
import os, glob
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Excel'i, çalışma klasöründen bağımsız olarak BU BETİĞİN yanında ara ---
HERE = os.path.dirname(os.path.abspath(__file__))
_cand = (glob.glob(os.path.join(HERE, "*MonteCarlo*guncel*.xlsx"))
         or glob.glob(os.path.join(HERE, "*MonteCarlo*.xlsx"))
         or glob.glob(os.path.join(HERE, "*.xlsx")))
if not _cand:
    print("Excel dosyası bu klasörde bulunamadı:", HERE)
    print("Klasördeki dosyalar:")
    for _f in os.listdir(HERE): print("  ", _f)
    raise SystemExit("Lütfen .xlsx dosyasını bu .py ile AYNI klasöre koyun.")
XLSX = _cand[0]
print("Kullanılan Excel:", os.path.basename(XLSX))

def OUT(name): return os.path.join(HERE, name)   # çıktıları da betiğin klasörüne kaydet
DPI  = 300
BL, AM, TL, RED, GRID = "#2e5a7a", "#c0702f", "#45897a", "#a32d2d", "#ededed"
plt.rcParams.update({"font.size":10,"savefig.dpi":DPI,"axes.edgecolor":"#555","axes.linewidth":0.8})

wb = openpyxl.load_workbook(XLSX, data_only=True)   # data_only: formül sonuçlarını okur

def col(ws, letter, r1, r2):
    return [ws[f"{letter}{r}"].value for r in range(r1, r2+1)]

def despine(ax):
    for s in ("top","right"): ax.spines[s].set_visible(False)

# ---------- Figure 1 — Pb konsantrasyonu (Konsantrasyon sayfası) ----------
k = wb["Konsantrasyon (ppm)"]
labs = [str(k[f"A{r}"].value).replace("Sup_","") for r in range(2,20)]
pb   = [k[f"F{r}"].value or 0 for r in range(2,20)]   # F = Pb
exc  = [v > 3.0 for v in pb]
fig, ax = plt.subplots(figsize=(9,3.3))
ax.bar(range(18), pb, width=0.74, zorder=3, edgecolor="white", linewidth=0.5,
       color=[AM if e else BL for e in exc])
ax.axhline(3.0, color=RED, ls="--", lw=1.2)
for i,v in enumerate(pb):
    if exc[i]: ax.text(i, v+0.1, f"{v:.2f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
ax.set_xticks(range(18)); ax.set_xticklabels(labs, fontsize=8.5)
ax.set_ylabel("Pb konsantrasyonu (mg/kg)"); ax.set_ylim(0,8.2)
ax.grid(axis="y", color="#e6e6e6", lw=0.6, zorder=0); despine(ax)
plt.tight_layout(); plt.savefig(OUT("Figure1_Pb.png"), bbox_inches="tight"); plt.close()

# ---------- Şekil 3 — deterministik HI (Deterministik sayfası, L sütunu) ----------
d = wb["Deterministik (nokta)"]
labs = [str(d[f"A{r}"].value).replace("Sup_","") for r in range(4,22)]
hi   = [d[f"L{r}"].value or 0 for r in range(4,22)]
fig, ax = plt.subplots(figsize=(9,3.3))
ax.bar(range(18), [max(v,1e-4) for v in hi], width=0.74, zorder=3, edgecolor="white", linewidth=0.5,
       color=[AM if v>=1 else BL for v in hi])
ax.axhline(1, color=RED, ls="--", lw=1.2); ax.set_yscale("log"); ax.set_ylim(1e-4,30)
ax.set_xticks(range(18)); ax.set_xticklabels(labs, fontsize=8.5)
ax.set_ylabel("Deterministik HI (log)")
ax.grid(axis="y", which="both", color=GRID, lw=0.6, zorder=0); despine(ax)
plt.tight_layout(); plt.savefig(OUT("Sekil3_det_HI.png"), bbox_inches="tight"); plt.close()

# ---------- Marka Bazlı MC sayfasından: ürünler, HI P95, As-CR P95, LB/MB/UB, katkı ----
m = wb["Marka Bazlı MC"]
prods = [m.cell(4, 2+j).value for j in range(6)]                 # B4:G4
hiP95 = [m.cell(14,2+j).value for j in range(6)]                 # HI P95 satırı
crP95 = [m.cell(15,2+j).value for j in range(6)]                 # As-CR P95 satırı
sensP = {m.cell(23+i,1).value: (m.cell(23+i,2).value, m.cell(23+i,3).value, m.cell(23+i,4).value)
         for i in range(6)}                                      # Ürün: (LB,MB,UB)
contrib = {m.cell(33+i,1).value: m.cell(33+i,2).value for i in range(8) if m.cell(33+i,1).value}

# ---------- Figure 2 — LB/MB/UB P95 HI ----------
import numpy as np
x = np.arange(len(prods))
fig, ax = plt.subplots(figsize=(9,3.6))
for idx,(nm,col_,ls,mk) in enumerate([("LB",BL,"-","o"),("MB","#7a6a55","--","s"),("UB",AM,":","o")]):
    y = [sensP[p][idx] for p in prods]
    ax.plot(x, y, ls, color=col_, marker=mk, ms=5, lw=1.6, label=f"{nm}")
ax.axhline(1, color="#333", ls="-.", lw=1.1, label="HI = 1")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(prods)
ax.set_ylabel("P95 HI (log)"); ax.grid(axis="y", which="both", color=GRID, lw=0.6)
ax.legend(fontsize=8.5, ncol=2, frameon=False, loc="upper right"); despine(ax)
plt.tight_layout(); plt.savefig(OUT("Figure2_LBMBUB.png"), bbox_inches="tight"); plt.close()

# ---------- Şekil 4 — MC P95 HI (marka) ----------
fig, ax = plt.subplots(figsize=(7,3.2))
ax.bar(x, hiP95, width=0.6, zorder=3, edgecolor="white",
       color=[AM if v>=1 else BL for v in hiP95])
ax.axhline(1, color=RED, ls="--", lw=1.2); ax.set_yscale("log")
ax.set_xticks(x); ax.set_xticklabels(prods); ax.set_ylabel("MC HI P95 (log)")
for i,v in enumerate(hiP95): ax.text(i,v,f"{v:.2f}" if v>=0.01 else f"{v:.0e}",ha="center",va="bottom",fontsize=8)
ax.grid(axis="y", which="both", color=GRID, lw=0.6, zorder=0); despine(ax)
plt.tight_layout(); plt.savefig(OUT("Sekil4_MC_HI.png"), bbox_inches="tight"); plt.close()

# ---------- Şekil 5 — MC P95 As-CR (marka) ----------
fig, ax = plt.subplots(figsize=(7,3.2))
ax.bar(x, [v if v and v>0 else 1e-8 for v in crP95], width=0.6, zorder=3, edgecolor="white",
       color=[AM if (v or 0)>=1e-4 else BL for v in crP95])
ax.axhline(1e-4, color=RED, ls="--", lw=1.2); ax.set_yscale("log"); ax.set_ylim(1e-6,3e-2)
ax.set_xticks(x); ax.set_xticklabels(prods); ax.set_ylabel("MC As-CR P95 (log)")
ax.grid(axis="y", which="both", color=GRID, lw=0.6, zorder=0); despine(ax)
plt.tight_layout(); plt.savefig(OUT("Sekil5_MC_CR.png"), bbox_inches="tight"); plt.close()

# ---------- Şekil 6 — element katkısı ----------
order = sorted(contrib, key=lambda e: contrib[e], reverse=True)
fig, ax = plt.subplots(figsize=(7,3.4))
ax.bar(range(len(order)), [max(contrib[e],1e-4) for e in order], width=0.6, zorder=3, edgecolor="white", color=TL)
ax.set_yscale("log"); ax.set_ylim(1e-4,300); ax.set_xticks(range(len(order))); ax.set_xticklabels(order)
ax.set_ylabel("HI'ya katkı (%, log)")
for i,e in enumerate(order): ax.text(i,max(contrib[e],1e-4),f"{contrib[e]:.2f}%" if contrib[e]>=0.01 else "0%",ha="center",va="bottom",fontsize=8)
ax.grid(axis="y", which="both", color=GRID, lw=0.6, zorder=0); despine(ax)
plt.tight_layout(); plt.savefig(OUT("Sekil6_katki.png"), bbox_inches="tight"); plt.close()

print("Bitti. Üretilen figürler:")
for f in ["Figure1_Pb","Figure2_LBMBUB","Sekil3_det_HI","Sekil4_MC_HI","Sekil5_MC_CR","Sekil6_katki"]:
    print("  "+f+".png")
