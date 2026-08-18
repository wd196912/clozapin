#!/usr/bin/env python3
"""
Recompute clozapine vs olanzapine disproportionality — composite endpoint edition.
Mirrors F:/clazpin/recompute_composite_endpoint.R (new signal rules):
  - ROR/PRR reported only when comparator events >= 5 (else "NR", IC only)
  - Positive signal = IC025 > 0 AND comparator events >= 5
  - Rows: composite endpoint (primary), 3 subgroups, 17 PTs, 9-additional aggregate
Clozapine counts: F:/clazpin/composite_endpoint_signal_results.csv (ASCII-derived)
Olanzapine counts: F:/clo-ola/data/olanzapine_ascii_extraction.json (ASCII-derived)
Both denominators computed: report level (unique primaryid) and caseid level.
"""

import json
import math
import os
import csv
from collections import Counter

DATA_DIR = "F:/clo-ola/data"
CLO_RESULTS = "F:/clazpin/composite_endpoint_signal_results.csv"
OLA_EXTRACT = os.path.join(DATA_DIR, "olanzapine_ascii_extraction.json")

CLOZAPINE_TOTAL = 44055

COMPOSITE_PTS = {
    "pneumonia", "pneumonia aspiration", "lower respiratory tract infection",
    "pneumonia bacterial", "pneumonia viral", "empyema", "lung abscess",
    "pneumonia klebsiella", "pneumonia influenzal", "pneumonia staphylococcal",
}
SPECTRUM_PTS = {p for p in COMPOSITE_PTS if p != "lower respiratory tract infection"}
LRTI_PTS = {"lower respiratory tract infection"}

PT17 = [
    "pneumonia", "pneumonia aspiration", "lower respiratory tract infection",
    "upper respiratory tract infection", "covid-19 pneumonia",
    "respiratory tract infection", "pneumonia bacterial", "pneumonia viral",
    "pneumonitis", "empyema", "idiopathic interstitial pneumonia",
    "lung abscess", "pulmonary tuberculosis", "pneumonia klebsiella",
    "respiratory tract infection viral", "pneumonia influenzal",
    "pneumonia staphylococcal",
]
EXTRA9 = [
    "pneumonia necrotising", "pneumonia anthrax", "pulmonary sepsis",
    "acute interstitial pneumonitis", "atypical pneumonia",
    "eosinophilic pneumonia", "hypersensitivity pneumonitis",
    "organising pneumonia", "pleural infection",
]
ALL26 = set(PT17) | set(EXTRA9)
SPECIAL_OTHER_PTS = {p for p in ALL26 if p not in COMPOSITE_PTS}


def load_clozapine_counts():
    clo = {}
    with open(CLO_RESULTS, "r", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            clo[row["Label"]] = int(row["a"])
    return clo


def set_count(pts_by_pid, pt_set):
    return sum(1 for pts in pts_by_pid.values() if any(p.lower() in pt_set for p in pts))


def calc(a, c, n_clo, n_ola, label):
    b = n_clo - a
    d = n_ola - c
    n = n_clo + n_ola
    report_ror = (c >= 5) and (a >= 1)

    ror = (a / b) / (c / d) if (b > 0 and c > 0 and d > 0) else float("inf")
    ror_raw = ror if (b > 0 and c > 0 and d > 0 and a > 0) else None

    if report_ror and a > 0 and c > 0:
        se_ror = math.sqrt(1/a + 1/b + 1/c + 1/d)
        ror_ci_low = math.exp(math.log(ror) - 1.96 * se_ror)
        ror_ci_high = math.exp(math.log(ror) + 1.96 * se_ror)
        prr = (a / (a + b)) / (c / (c + d))
        se_prr = math.sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))
        prr_ci_low = math.exp(math.log(prr) - 1.96 * se_prr)
        prr_ci_high = math.exp(math.log(prr) + 1.96 * se_prr)
        p_val = 2 * (1 - _norm_cdf(abs(math.log(ror) / se_ror)))
    else:
        ror = ror_ci_low = ror_ci_high = prr = prr_ci_low = prr_ci_high = None
        p_val = None

    c_exp = (n_clo * (a + c)) / n
    ic_val = math.log2((a + 0.5) / (c_exp + 0.5))
    ic_sd = math.sqrt(1 / (math.log(2)**2) * (1/(a + 0.5) + 1/(c_exp + 0.5)))
    ic025 = ic_val - 1.96 * ic_sd
    ic975 = ic_val + 1.96 * ic_sd
    positive = (c >= 5) and (ic025 > 0)

    return {
        "label": label, "a": a, "b": b, "c": c, "d": d, "n_total": n,
        "ror": round(ror, 2) if ror is not None else "NR",
        "ror_ci_95_low": round(ror_ci_low, 2) if ror_ci_low is not None else "NR",
        "ror_ci_95_high": round(ror_ci_high, 2) if ror_ci_high is not None else "NR",
        "prr": round(prr, 2) if prr is not None else "NR",
        "prr_ci_95_low": round(prr_ci_low, 2) if prr_ci_low is not None else "NR",
        "prr_ci_95_high": round(prr_ci_high, 2) if prr_ci_high is not None else "NR",
        "ic": round(ic_val, 3), "ic025": round(ic025, 3), "ic975": round(ic975, 3),
        "p_value": f"{p_val:.3g}" if p_val is not None else "NR",
        "signal_positive": positive,
        "ror_raw": round(ror_raw, 4) if ror_raw is not None else None,
    }


def _norm_cdf(x):
    # Abramowitz-Stegun 7.1.26 approximation of the standard normal CDF
    t = 1 / (1 + 0.2316419 * abs(x))
    poly = 0.319381530*t - 0.356563782*t**2 + 1.781477937*t**3 \
           - 1.821255978*t**4 + 1.330274429*t**5
    pdf = 0.3989422804014327 * math.exp(-x*x/2)
    cdf = 1 - pdf * poly
    return cdf if x >= 0 else 1 - cdf


def run_variant(name, clo_counts, ola_pt_counts, ola_pts_by_pid, n_ola, total_pulm):
    rows_def = [
        ("Composite pulmonary infection endpoint (primary)", COMPOSITE_PTS, None),
        ("Pneumonia spectrum (subgroup)", SPECTRUM_PTS, None),
        ("Non-pneumonia lower respiratory tract infection (subgroup)", LRTI_PTS, None),
        ("Special/other pulmonary infection terms (subgroup)", SPECIAL_OTHER_PTS, None),
    ]
    for p in PT17:
        rows_def.append((p, {p}, p))
    rows_def.append(("Other pulmonary infection terms (9 additional PTs)", set(EXTRA9), None))

    results = []
    for label, pt_set, pt_key in rows_def:
        a = clo_counts[label]
        if pt_key:
            c = ola_pt_counts.get(pt_key.lower(), 0)
        else:
            c = set_count(ola_pts_by_pid, pt_set)
        results.append(calc(a, c, CLOZAPINE_TOTAL, n_ola, label))

    signals = [r for r in results if r["signal_positive"]]
    print(f"\n[{name}] olanzapine denominator N={n_ola} (pulmonary reports={total_pulm})")
    hdr = f"{'Row':<58} {'CLO':>6} {'OLA':>6} {'ROR':>8} {'95% CI':>22} {'PRR':>8} {'IC':>8} {'IC025':>8} {'Sig':>4}"
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        ci = f"({r['ror_ci_95_low']}–{r['ror_ci_95_high']})" if r["ror_ci_95_low"] != "NR" else "(NR)"
        print(f"{r['label']:<58} {r['a']:>6} {r['c']:>6} {str(r['ror']):>8} {ci:>22} "
              f"{str(r['prr']):>8} {r['ic']:>8} {r['ic025']:>8} {'***' if r['signal_positive'] else '':>4}")
    print(f"Positive signals: {len(signals)}")
    for s in signals:
        print(f"  *** {s['label']}: ROR={s['ror']}, IC025={s['ic025']}")
    return results


def main():
    clo_counts = load_clozapine_counts()
    print(f"Loaded {len(clo_counts)} clozapine counts from {CLO_RESULTS}")

    with open(OLA_EXTRACT, "r", encoding="utf-8") as f:
        ola = json.load(f)

    ola_pt_counts_caseid = {k.lower(): v for k, v in ola["pt_counts"].items()}
    ola_pts_by_pid_caseid = ola["pts_by_pid"]
    n_ola_caseid = ola["total_ps_reports"]
    total_pulm_caseid = ola["total_pulmonary_reports"]

    ola_pts_by_pid_all = ola["pts_by_pid_all"]
    ola_pt_counts_all = Counter()
    for pts in ola_pts_by_pid_all.values():
        for p in set(pts):
            ola_pt_counts_all[p.lower()] += 1
    n_ola_all = ola["total_ps_primaryid_level"]
    total_pulm_all = len(ola_pts_by_pid_all)

    print(f"Olanzapine extraction loaded: caseid-level N={n_ola_caseid}, "
          f"primaryid-level N={n_ola_all}")

    res_caseid = run_variant("CASEID-LEVEL", clo_counts, ola_pt_counts_caseid,
                             ola_pts_by_pid_caseid, n_ola_caseid, total_pulm_caseid)
    res_all = run_variant("PRIMARYID-LEVEL", clo_counts, ola_pt_counts_all,
                          ola_pts_by_pid_all, n_ola_all, total_pulm_all)

    out = {
        "source": "FAERS ASCII quarterly files (2022Q1–2025Q4)",
        "method": "composite endpoint + new signal rule (IC025>0 AND comparator n>=5); "
                  "ROR/PRR reported only when comparator n>=5",
        "clozapine_total_ps": CLOZAPINE_TOTAL,
        "comparison": "clozapine_vs_olanzapine",
        "caseid_level": {"olanzapine_total": n_ola_caseid, "results": res_caseid},
        "primaryid_level": {"olanzapine_total": n_ola_all, "results": res_all},
    }
    out_path = os.path.join(DATA_DIR, "disproportionality_clo_vs_ola_ascii.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"\nSaved to {out_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
