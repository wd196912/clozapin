#!/usr/bin/env python3
"""
Recompute clozapine vs olanzapine disproportionality using ASCII-derived counts.
Clozapine: 44,055 PS reports (FAERS ASCII, same pipeline)
Olanzapine: 12,701 PS reports (FAERS ASCII, just extracted)

Replaces the openFDA-derived disproportionality_clo_vs_ola.json with
ASCII-consistent results.
"""

import json
import math
import os

DATA_DIR = "F:/clo-ola/data"
OUT_DIR = "F:/clo-ola/data"

# Clozapine ASCII totals (from the existing signal detection pipeline)
CLOZAPINE_TOTAL = 44055

# Olanzapine ASCII total (from extract_olanzapine_ascii.py)
OLANZAPINE_TOTAL = 12701

# Pulmonary infection PT list (matching the manuscript)
PULMONARY_PTS = [
    "pneumonia", "pneumonia aspiration", "lower respiratory tract infection",
    "upper respiratory tract infection", "covid-19 pneumonia",
    "respiratory tract infection", "pneumonia bacterial", "pneumonia viral",
    "pneumonitis", "empyema", "pulmonary tuberculosis",
    "pneumonia klebsiella", "respiratory tract infection viral",
    "pneumonia influenzal", "pneumonia staphylococcal",
]

def main():
    # ── Load clozapine ASCII PT counts ──
    # From the existing clozapine_lung_signal_detection.csv (clozapine vs risperidone)
    # The Clozapine_Event column gives us per-PT event counts from ASCII
    import csv

    clo_pt_counts = {}
    with open("F:/clazpin/clozapine_lung_signal_detection.csv", "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pt_lower = row["PT_Term"].strip().lower()
            event_count = int(row["Clozapine_Event"])
            clo_pt_counts[pt_lower] = event_count

    print(f"Loaded {len(clo_pt_counts)} PT counts for clozapine (ASCII, N={CLOZAPINE_TOTAL})")

    # ── Load olanzapine ASCII PT counts ──
    # Use the full extraction JSON (mixed-case PT keys) and normalize to lowercase
    with open(os.path.join(DATA_DIR, "olanzapine_ascii_extraction.json"), "r", encoding="utf-8") as f:
        ola_raw = json.load(f)
    ola_pt_counts = {k.lower(): v for k, v in ola_raw["pt_counts"].items()}

    print(f"Loaded {len(ola_pt_counts)} PT counts for olanzapine (ASCII, N={OLANZAPINE_TOTAL})")

    # ── Compute disproportionality ──
    results = []
    signals = []

    for pt in PULMONARY_PTS:
        a = clo_pt_counts.get(pt, 0)
        c = ola_pt_counts.get(pt, 0)
        b = CLOZAPINE_TOTAL - a
        d = OLANZAPINE_TOTAL - c

        # ROR with 95% CI
        if a > 0 and b > 0 and d > 0:
            if c > 0:
                ror = (a * d) / (b * c)
                se_ror = math.sqrt(1/a + 1/b + 1/c + 1/d)
                ror_ci_low = math.exp(math.log(ror) - 1.96 * se_ror)
                ror_ci_high = math.exp(math.log(ror) + 1.96 * se_ror)
                signal_ror = ror_ci_low > 1
            else:
                ror = float('inf')
                ror_ci_low = float('inf')
                ror_ci_high = float('inf')
                signal_ror = True  # positive if events only in clozapine
        elif a > 0 and c == 0:
            ror = float('inf')
            ror_ci_low = float('inf')
            ror_ci_high = float('inf')
            signal_ror = True
        else:
            ror = 0
            ror_ci_low = 0
            ror_ci_high = 0
            signal_ror = False

        # PRR
        if a > 0:
            if c > 0:
                prr = (a / CLOZAPINE_TOTAL) / (c / OLANZAPINE_TOTAL)
            else:
                prr = (a / CLOZAPINE_TOTAL) / (0.5 / OLANZAPINE_TOTAL)  # Haldane
        else:
            prr = 0

        # IC (BCPNN)
        n_pt_total = a + c
        n_all_total = CLOZAPINE_TOTAL + OLANZAPINE_TOTAL
        c_exp = (CLOZAPINE_TOTAL * n_pt_total) / n_all_total if n_all_total > 0 else 0.001

        if a > 0 and c_exp > 0:
            ic_val = math.log2((a + 0.5) / (c_exp + 0.5))
            var_ic = 1 / (math.log(2)**2) * (1/(a + 0.5) + 1/(c_exp + 0.5))
            ic025 = ic_val - 1.96 * math.sqrt(var_ic)
        else:
            ic_val, ic025 = 0, -10

        # Final signal determination (dual criterion, matching manuscript Methods)
        # For PTs with comparator events: ROR lower CI > 1 AND IC025 > 0
        # For PTs with zero comparator events: IC > 0 AND IC025 > 0
        # (ROR=∞ alone is insufficient; IC025 ≤ 0 = no detectable Bayesian signal)
        if c > 0:
            is_signal = ror_ci_low > 1 and ic025 > 0
        else:
            is_signal = ic_val > 0 and ic025 > 0

        result = {
            "pt": pt,
            "n_clozapine": a,
            "n_olanzapine": c,
            "clozapine_total": CLOZAPINE_TOTAL,
            "olanzapine_total": OLANZAPINE_TOTAL,
            "ror": round(ror, 2) if ror != float('inf') else "Inf",
            "ror_ci_95_low": round(ror_ci_low, 2) if ror_ci_low != float('inf') else "Inf",
            "ror_ci_95_high": round(ror_ci_high, 2) if ror_ci_high != float('inf') else "Inf",
            "prr": round(prr, 2) if prr != float('inf') else "Inf",
            "ic": round(ic_val, 3),
            "ic025": round(ic025, 3),
            "signal_positive": is_signal,
        }
        results.append(result)

        if is_signal:
            signals.append(result)

    # ── Print results ──
    print(f"\n{'='*100}")
    print("DISPROPORTIONALITY ANALYSIS: Clozapine vs Olanzapine (FAERS ASCII, 2022Q1–2025Q4)")
    print(f"Clozapine PS reports: {CLOZAPINE_TOTAL}  |  Olanzapine PS reports: {OLANZAPINE_TOTAL}")
    print(f"{'='*100}")
    header = f"{'PT':<42} {'CLO':>6} {'OLA':>6} {'ROR':>10} {'95% CI':>24} {'PRR':>8} {'IC':>8} {'IC025':>8} {'Signal':>8}"
    print(header)
    print("-" * 100)

    for r in results:
        ror_s = f"{r['ror']}" if isinstance(r['ror'], str) else f"{r['ror']:.2f}"
        if r['ror_ci_95_low'] == 'Inf':
            ci_s = "(Inf–Inf)"
        else:
            ci_s = f"({r['ror_ci_95_low']:.2f}–{r['ror_ci_95_high']:.2f})"
        sig = "***" if r['signal_positive'] else ""
        print(f"{r['pt']:<42} {r['n_clozapine']:>6} {r['n_olanzapine']:>6} {ror_s:>10} {ci_s:>24} {r['prr']:>8.2f} {r['ic']:>8.3f} {r['ic025']:>8.3f} {sig:>8}")

    print(f"\nTotal PTs assessed: {len(results)}")
    print(f"Positive signals:   {len(signals)}")
    for s in signals:
        ror_s = f"ROR={s['ror']}" if isinstance(s['ror'], str) else f"ROR={s['ror']:.2f}"
        print(f"  *** {s['pt']}: {ror_s}, 95%CI=({s['ror_ci_95_low']}–{s['ror_ci_95_high']}), IC025={s['ic025']}")

    # ── Save updated JSON ──
    out_path = os.path.join(OUT_DIR, "disproportionality_clo_vs_ola_ascii.json")
    output = {
        "source": "FAERS ASCII quarterly files (2022Q1–2025Q4)",
        "method": "consistent with clozapine primary analysis pipeline",
        "clozapine_total_ps": CLOZAPINE_TOTAL,
        "olanzapine_total_ps": OLANZAPINE_TOTAL,
        "comparison": "clozapine_vs_olanzapine",
        "n_pts": len(results),
        "n_signals": len(signals),
        "results": results,
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)
    print(f"\nSaved to {out_path}")

    # Also update the canonical disproportionality_clo_vs_ola.json path
    # (the generate_all_tables_figures.py script reads from this path)
    canon_path = os.path.join(OUT_DIR, "disproportionality_clo_vs_ola.json")
    # Back up the openFDA version
    backup_path = os.path.join(OUT_DIR, "disproportionality_clo_vs_ola_openfda_backup.json")
    if os.path.exists(canon_path):
        os.replace(canon_path, backup_path)
        print(f"Backed up openFDA version to {backup_path}")
    with open(canon_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"Updated canonical file: {canon_path}")

    print("\nDone!")
    return results

if __name__ == "__main__":
    main()
