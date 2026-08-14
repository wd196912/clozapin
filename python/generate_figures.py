#!/usr/bin/env python
"""Generate publication-quality figures for clozapine FAERS manuscript."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import csv, os

OUT = "F:/clazpin/manuscript/analysis/figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
})

# ============================================================
# Figure 1 — Case Selection Flowchart (PRISMA-style)
# ============================================================
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

boxes = [
    (5, 9.2, 'FAERS ASCII database\n2022Q1 – 2025Q4\nAll reports', 'lightgray'),
    (5, 7.8, 'Clozapine primary suspect (PS) DRUG records\nN = 45,749', 'white'),
    (5, 6.4, 'Unique reports after de-duplication\nN = 44,055\n(1,694 duplicate drug records removed)', 'white'),
    (5, 5.0, 'Pulmonary infection PTs\n(17 MedDRA PTs)\nN = 1,305 unique reports', 'white'),
    (2.5, 3.4, 'Signal Detection Analysis\nClozapine vs Risperidone\n17 PTs assessed', '#e8f5e9'),
    (7.5, 3.4, 'Fatal Outcome Analysis\n309 fatal / 996 non-fatal\nMultivariable logistic regression', '#e3f2fd'),
    (2.5, 1.8, 'Excluded: duplicate drug records\nN = 1,694', '#ffebee'),
    (7.5, 1.8, 'Dose-Response Subset\n405 with valid dose data\nRCS nonlinear modeling', '#fff3e0'),
]

for x, y, text, color in boxes:
    w, h = (4.2, 1.0) if x in (2.5, 7.5) else (5.8, 1.0)
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, facecolor=color,
                                edgecolor='black', linewidth=1.5,
                                boxstyle='round,pad=0.1'))
    ax.text(x, y, text, ha='center', va='center', fontsize=9, weight='bold')

# Arrows (vertical spine)
for y1, y2 in [(8.7, 8.3), (7.3, 6.9), (5.9, 5.5)]:
    ax.annotate('', xy=(5, y1), xytext=(5, y2),
                arrowprops=dict(arrowstyle='->', lw=2, color='gray'))

# Branch arrows
ax.annotate('', xy=(2.5, 3.9), xytext=(5, 4.5),
            arrowprops=dict(arrowstyle='->', lw=1.5, color='green'))
ax.annotate('', xy=(7.5, 3.9), xytext=(5, 4.5),
            arrowprops=dict(arrowstyle='->', lw=1.5, color='blue'))
# Excluded arrow
ax.annotate('', xy=(2.5, 2.3), xytext=(5, 5.5),
            arrowprops=dict(arrowstyle='->', lw=1, color='red', ls='dashed'))

ax.set_title('Figure 1. Case Selection and Analysis Flowchart', fontsize=13, weight='bold', y=1.02)
fig.tight_layout()
fig.savefig(f'{OUT}/fig1_flowchart.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print('[OK] Fig 1: Flowchart saved')

# ============================================================
# Figure 2 — Signal Detection Forest Plot
# ============================================================
signal_data = [
    ("Lower respiratory tract infection", 81.83, 11.48, 583.38, 0.419),
    ("Pneumonia", 3.35, 2.64, 4.25, 0.283),
    ("Pneumonia aspiration", 1.80, 1.28, 2.52, 0.173),
    ("Upper respiratory tract infection", 4.64, 1.10, 19.51, 0.316),
    ("COVID-19 pneumonia", 4.29, 1.02, 18.13, 0.308),
    ("Pneumonia bacterial", 1.37, 0.39, 4.87, 0.100),
    ("Pneumonia viral", None, None, None, 0.393),
    ("Pneumonitis", None, None, None, 0.412),
    ("Empyema", None, None, None, 0.410),
    ("Pulmonary tuberculosis", None, None, None, None),
]

fig, ax = plt.subplots(1, 1, figsize=(9, 5.5))
y_positions = list(range(len(signal_data)))
pt_labels = [s[0] for s in signal_data]
rors = [s[1] for s in signal_data]
ci_low = [s[2] for s in signal_data]
ci_high = [s[3] for s in signal_data]

colors = []
for i, s in enumerate(signal_data):
    if s[1] is not None and s[2] is not None:
        colors.append('#c62828' if s[2] > 1 and s[4] is not None and s[4] > 0 else '#1565c0')
    else:
        colors.append('#9e9e9e')

# Plot ROR with CI
for i, (ror, lo, hi, col) in enumerate(zip(rors, ci_low, ci_high, colors)):
    if ror is not None and lo is not None:
        ax.errorbar(np.log(ror), i, xerr=[[np.log(ror/lo)], [np.log(hi/ror)]],
                    fmt='o', color=col, capsize=3, markersize=8, capthick=1.5, elinewidth=2)
    else:
        ax.plot(0, i, 'X', color='#9e9e9e', markersize=8, markeredgewidth=2)

# Reference line
ax.axvline(0, color='gray', linestyle='--', alpha=0.7)

# ROR labels
for i, (ror, lo, hi, col) in enumerate(zip(rors, ci_low, ci_high, colors)):
    if ror is not None and ror > 50:
        ror_label = f'ROR={ror:.0f}'
    elif ror is not None:
        ror_label = f'ROR={ror:.2f}'
    else:
        ror_label = 'ROR=Inf†'

    offset = 0.35 if ror and ror > 10 else 0.2
    ax.text(np.log(ror) + offset if ror else offset, i + 0.25,
            ror_label, fontsize=7, va='bottom', color=col, weight='bold')

# Also annotate IC for inf RORs
for i, s in enumerate(signal_data):
    if s[1] is None and s[4] is not None:
        ax.text(0.5, i - 0.25, f'IC={s[4]:.3f}', fontsize=7, color='#9e9e9e')

ax.set_yticks(y_positions)
ax.set_yticklabels(pt_labels, fontsize=9)
ax.set_xlabel('ln(ROR)', fontsize=10)
ax.set_title('Figure 2. Disproportionality Analysis: Clozapine vs Risperidone\nPulmonary Infection Preferred Terms in FAERS (2022–2025)',
             fontsize=13, weight='bold')
ax.invert_yaxis()

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#c62828', markersize=10, label='Positive Signal (IC > 0)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#1565c0', markersize=10, label='No Signal'),
    Line2D([0], [0], marker='X', color='#9e9e9e', markersize=8, label='ROR=Inf (Zero-cell in comparator)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

ax.set_xlim(-2, 8)
fig.tight_layout()
fig.savefig(f'{OUT}/fig2_signal_forest.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print('[OK] Fig 2: Signal detection forest plot saved')

# ============================================================
# Figure 4 — Fatal Risk Factors Forest Plot (Multivariable)
# ============================================================
risk_data = [
    ("Age (per 10 years)", 1.90, 1.64, 2.19, 0.0000001, True),
    ("Aspiration pneumonia", 2.40, 1.58, 3.66, 0.0000424, True),
    ("Male sex (vs female)", 0.99, 0.67, 1.45, 0.947, False),
    ("Europe region (vs other)", 0.75, 0.52, 1.10, 0.138, False),
    ("Year (per year)", 0.97, 0.90, 1.04, 0.368, False),
]

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
y_pos = list(range(len(risk_data)))

for i, (label, or_, lo, hi, p, sig) in enumerate(risk_data):
    col = '#c62828' if sig else '#455a64'
    ax.errorbar(np.log(or_), i, xerr=[[np.log(or_/lo)], [np.log(hi/or_)]],
                fmt='o', color=col, capsize=4, markersize=10, capthick=2, elinewidth=2.5)
    # OR label
    ax.text(np.log(or_) + 0.08, i + 0.28, f'OR={or_:.2f}', fontsize=9, color=col, weight='bold')
    # p-value
    p_label = 'p<0.001' if p < 0.001 else f'p={p:.3f}'
    ax.text(1.2, i - 0.28, p_label, fontsize=8, color='#757575')

ax.axvline(1, color='gray', linestyle='--', alpha=0.7, label='OR = 1 (null)')
ax.set_yticks(y_pos)
ax.set_yticklabels([r[0] for r in risk_data], fontsize=10)
ax.set_xlabel('Adjusted Odds Ratio (95% CI)', fontsize=10)
ax.set_title('Figure 4. Multivariable Logistic Regression:\nRisk Factors for Fatal Outcome in Clozapine-Associated Pulmonary Infection',
             fontsize=12, weight='bold')
ax.invert_yaxis()
ax.set_xlim(-0.5, 1.5)

# Tick labels in OR scale
or_ticks = [0.5, 1.0, 1.5, 2.0, 2.5]
ax.set_xticks([np.log(t) for t in or_ticks])
ax.set_xticklabels([f'{t:.1f}' for t in or_ticks])

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#c62828', markersize=10, label='p < 0.05'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#455a64', markersize=10, label='Not significant'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

fig.tight_layout()
fig.savefig(f'{OUT}/fig4_fatal_forest.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print('[OK] Fig 4: Fatal risk factors forest plot saved')

# ============================================================
# Figure 5 — Geographic Fatality Distribution
# ============================================================
# Countries with >=10 clozapine pulmonary infection reports (OUTC outc_cod="DE")
geo_data = [
    ("New Zealand", 56, 80, 70.0),
    ("Canada", 102, 223, 45.7),
    ("Ireland", 9, 35, 25.7),
    ("Australia", 18, 77, 23.4),
    ("France", 5, 23, 21.7),
    ("United Kingdom", 83, 460, 18.0),
    ("United States", 22, 237, 9.3),
    ("Argentina", 1, 11, 9.1),
    ("Germany", 2, 30, 6.7),
    ("Finland", 1, 28, 3.6),
    ("European Union (unspec.)", 0, 30, 0.0),
    ("Portugal", 0, 10, 0.0),
]

fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))
labels = [g[0] for g in geo_data]
pcts = [g[3] for g in geo_data]
y_pos = np.arange(len(labels))
colors5 = ['#c62828' if p >= 40 else '#ef6c00' if p >= 20 else '#1565c0' for p in pcts]
bars = ax.barh(y_pos, pcts, color=colors5, edgecolor='black', linewidth=0.6)
for i, (name, f, n, p) in enumerate(geo_data):
    ax.text(p + 0.8, i, f'{p:.1f}%  ({f}/{n})', va='center', fontsize=8)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel('Fatality (%)', fontsize=10)
ax.set_xlim(0, 78)
ax.set_title('Figure 5. Country-Level Fatality of Clozapine-Associated\nPulmonary Infections (countries with ≥10 reports)',
             fontsize=12, weight='bold')
fig.tight_layout()
fig.savefig(f'{OUT}/fig5_geographic.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print('[OK] Fig 5: Geographic fatality distribution saved')

# ============================================================
# Summary
# ============================================================
print(f'\nAll figures saved to {OUT}/')
print('  fig1_flowchart.png')
print('  fig2_signal_forest.png')
print('  fig3_dose_response.png      (from existing)')
print('  fig4_fatal_forest.png')
print('  fig5_geographic.png')
print('  figS1_subgroup_forest.png   (from existing)')
