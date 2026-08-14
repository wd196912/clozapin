#!/usr/bin/env python3
"""
Extract olanzapine primary suspect reports from FAERS ASCII quarterly files.
Matches the clozapine ASCII extraction methodology: case-insensitive drug name
matching, caseid+caseversion de-duplication, primary suspect filter.
Outputs counts comparable to the clozapine ASCII-derived data.

Source: FAERS ASCII 2022Q1–2025Q4 at F:/faersdata/
"""

import csv
import json
import os
import re
import sys
import zipfile
from collections import defaultdict

DATA_DIR = "F:/faersdata"
OUT_DIR = "F:/clo-ola/data"
QUARTERS = [
    ("2022Q1", "22Q1"), ("2022Q2", "22Q2"), ("2022Q3", "22Q3"), ("2022Q4", "22Q4"),
    ("2023Q1", "23Q1"), ("2023Q2", "23Q2"), ("2023Q3", "23Q3"), ("2023Q4", "23Q4"),
    ("2024Q1", "24Q1"), ("2024Q2", "24Q2"), ("2024Q3", "24Q3"), ("2024Q4", "24Q4"),
    ("2025Q1", "25Q1"), ("2025Q2", "25Q2"), ("2025Q3", "25Q3"), ("2025Q4", "25Q4"),
]

# Olanzapine keywords — case-insensitive matching against drugname and prod_ai
OLA_KEYWORDS = re.compile(
    r"OLANZAPINE|ZYPREXA|SYMBYAX|ZYPADHERA|OLANZAPINE PAMOATE|ZYPREXA ZYDIS",
    re.IGNORECASE
)

# Exclude keywords (prevent false matches — e.g. queries that mention olanzapine
# but are actually for other drugs)
EXCLUDE_KEYWORDS = re.compile(
    r"NOT OLANZAPINE|OLANZAPINE LEVEL|OLANZAPINE CONCENTRATION",
    re.IGNORECASE
)

# Pulmonary infection PTs (matching the clozapine manuscript signal detection PT list)
PULMONARY_PTS_LOWER = {
    pt.lower() for pt in [
        "pneumonia", "pneumonia aspiration", "lower respiratory tract infection",
        "upper respiratory tract infection", "covid-19 pneumonia",
        "respiratory tract infection", "pneumonia bacterial", "pneumonia viral",
        "pneumonitis", "empyema", "pulmonary tuberculosis",
        "pneumonia klebsiella", "respiratory tract infection viral",
        "pneumonia influenzal", "pneumonia staphylococcal",
    ]
}


def open_ascii_file(quarter, table):
    """Open a FAERS ASCII table file, handling both extracted dirs and ZIPs."""
    label = quarter[1]  # e.g. "22Q1"
    fname = f"{table}{label}.txt"
    dir_path = os.path.join(DATA_DIR, f"faers_ascii_{quarter[0].lower()}")
    zip_path = os.path.join(DATA_DIR, f"faers_ascii_{quarter[0].lower()}.zip")

    # Prefer extracted directory
    extracted_file = os.path.join(dir_path, "ASCII", fname)
    if os.path.exists(extracted_file):
        return open(extracted_file, "r", encoding="latin-1", errors="replace")

    # Fall back to reading from ZIP
    if os.path.exists(zip_path):
        zf = zipfile.ZipFile(zip_path, "r")
        zip_member = f"ASCII/{fname}"
        if zip_member in zf.NameToInfo:
            return zf.open(zip_member, "r")
        # Some ZIPs use different internal structure
        for name in zf.namelist():
            if name.endswith(fname):
                return zf.open(name, "r")

    print(f"  WARNING: {fname} not found in {dir_path} or {zip_path}")
    return None


def iter_tsv(fileobj, expected_cols):
    """Iterate over a $-delimited FAERS file, yielding dicts."""
    if fileobj is None:
        return
    # Detect binary stream (from ZIP) and wrap in TextIOWrapper
    first_byte = fileobj.read(0)
    if isinstance(first_byte, bytes):
        import io
        fileobj = io.TextIOWrapper(fileobj, encoding="latin-1", errors="replace")
    reader = csv.DictReader(
        (line for line in fileobj if line.strip()),
        delimiter="$",
        quoting=csv.QUOTE_NONE,
    )
    for row in reader:
        yield row


def extract_quarter(quarter):
    """Extract olanzapine PS reports from one quarter. Returns:
    - drug_records: list of {primaryid, caseid, drugname, role_cod}
    - demo_records: dict primaryid -> {caseid, caseversion, sex, age, age_cod, reporter_country, wt, wt_cod}
    - reac_records: dict primaryid -> list of PT strings
    """
    label = quarter[0]
    print(f"  {label} ...", end=" ", flush=True)

    # ── Pass 1: read DRUG, find olanzapine PS primaryids ──
    drug_file = open_ascii_file(quarter, "DRUG")
    ola_primaryids = set()  # primaryids for olanzapine PS records
    primaryid_to_caseid = {}  # for joining to DEMO/REAC

    count_drug_total = 0
    count_ola_any = 0
    for row in iter_tsv(drug_file, 20):
        count_drug_total += 1
        drugname = row.get("drugname", "")
        prod_ai = row.get("prod_ai", "")
        combined = f"{drugname} {prod_ai}"

        if OLA_KEYWORDS.search(combined) and not EXCLUDE_KEYWORDS.search(combined):
            count_ola_any += 1
            role = row.get("role_cod", "").strip().upper()
            pid = row.get("primaryid", "").strip()
            cid = row.get("caseid", "").strip()
            if role == "PS" and pid:
                ola_primaryids.add(pid)
                if pid not in primaryid_to_caseid:
                    primaryid_to_caseid[pid] = cid

    if hasattr(drug_file, "close"):
        drug_file.close()

    if not ola_primaryids:
        print(f"0 PS records ({count_ola_any} any-role from {count_drug_total} drug rows)")
        return {}, {}, {}

    # ── Pass 2: read DEMO for olanzapine primaryids ──
    demo_file = open_ascii_file(quarter, "DEMO")
    demo_records = {}
    for row in iter_tsv(demo_file, 24):
        pid = row.get("primaryid", "").strip()
        if pid in ola_primaryids:
            demo_records[pid] = {
                "caseid": row.get("caseid", "").strip(),
                "caseversion": row.get("caseversion", "").strip(),
                "sex": row.get("sex", "").strip(),
                "age": row.get("age", "").strip(),
                "age_cod": row.get("age_cod", "").strip(),
                "reporter_country": row.get("reporter_country", "").strip(),
                "wt": row.get("wt", "").strip(),
                "wt_cod": row.get("wt_cod", "").strip(),
                "event_dt": row.get("event_dt", "").strip(),
                "occr_country": row.get("occr_country", "").strip(),
            }
    if hasattr(demo_file, "close"):
        demo_file.close()

    # ── Pass 3: read REAC for olanzapine primaryids (pulmonary PTs only) ──
    reac_file = open_ascii_file(quarter, "REAC")
    reac_records = defaultdict(list)
    for row in iter_tsv(reac_file, 4):
        pid = row.get("primaryid", "").strip()
        if pid in ola_primaryids:
            pt = row.get("pt", "").strip()
            if pt.lower() in PULMONARY_PTS_LOWER:
                reac_records[pid].append(pt)
    if hasattr(reac_file, "close"):
        reac_file.close()

    print(f"{len(ola_primaryids)} PS records ({count_ola_any} any-role from {count_drug_total} drug rows)")
    return ola_primaryids, demo_records, dict(reac_records)


def main():
    print("=" * 70)
    print("FAERS ASCII Olanzapine Extraction")
    print("Period: 2022Q1 – 2025Q4")
    print("=" * 70)

    all_demo = {}       # primaryid -> demo dict
    all_reac = {}       # primaryid -> [pt, ...]
    quarter_counts = [] # per-quarter stats

    for q in QUARTERS:
        pids, demos, reacs = extract_quarter(q)
        all_demo.update(demos)
        all_reac.update(reacs)
        quarter_counts.append((q[0], len(pids)))

    print(f"\nTotal unique primaryids before dedup: {len(all_demo)}")

    # ── Deduplication: for each caseid, keep highest caseversion ──
    print("\nDe-duplicating (caseid + caseversion)...")
    caseid_best = {}  # caseid -> (primaryid, caseversion)
    for pid, demo in all_demo.items():
        cid = demo["caseid"]
        cv_str = demo["caseversion"]
        try:
            cv = int(cv_str) if cv_str else 0
        except ValueError:
            cv = 0
        if cid not in caseid_best or cv > caseid_best[cid][1]:
            caseid_best[cid] = (pid, cv)

    dedup_pids = {v[0] for v in caseid_best.values()}
    dup_count = len(all_demo) - len(dedup_pids)
    print(f"  Unique caseids: {len(caseid_best)}")
    print(f"  Duplicates removed: {dup_count}")

    # ── Count by reporter_country ──
    country_counts = defaultdict(int)
    for pid in dedup_pids:
        demo = all_demo.get(pid, {})
        country = demo.get("reporter_country", "Unknown") or "Unknown"
        country_counts[country] += 1

    # ── Count by sex ──
    sex_counts = defaultdict(int)
    for pid in dedup_pids:
        demo = all_demo.get(pid, {})
        sex = demo.get("sex", "Unknown") or "Unknown"
        sex_counts[sex] += 1

    # ── Count pulmonary infection PTs (per primaryid, deduplicated) ──
    pt_counts = defaultdict(int)
    pulmonary_pids = set()
    for pid in dedup_pids:
        pts_for_pid = all_reac.get(pid, [])
        if pts_for_pid:
            pulmonary_pids.add(pid)
        for pt in pts_for_pid:
            pt_counts[pt] += 1

    # ── Summary ──
    total_ps = len(dedup_pids)
    total_pulm = len(pulmonary_pids)

    print(f"\n{'='*70}")
    print("RESULTS: Olanzapine FAERS ASCII Extraction (2022Q1–2025Q4)")
    print(f"{'='*70}")
    print(f"Total olanzapine PS reports (deduplicated): {total_ps}")
    print(f"Pulmonary infection reports:                {total_pulm}")
    print(f"Non-pulmonary reports:                      {total_ps - total_pulm}")
    print(f"\nPer-quarter:")
    for qname, cnt in quarter_counts:
        print(f"  {qname}: {cnt}")
    print(f"\nTop reporter countries:")
    for country, cnt in sorted(country_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {country}: {cnt} ({cnt/total_ps*100:.1f}%)")
    print(f"\nSex distribution:")
    for sex, cnt in sorted(sex_counts.items(), key=lambda x: -x[1]):
        print(f"  {sex}: {cnt} ({cnt/total_ps*100:.1f}%)")
    print(f"\nPulmonary infection PT counts:")
    for pt in sorted(pt_counts.keys()):
        print(f"  {pt}: {pt_counts[pt]}")
    print(f"\nTop 10 PTs:")
    for pt, cnt in sorted(pt_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {pt}: {cnt}")

    # ── Save results ──
    os.makedirs(OUT_DIR, exist_ok=True)

    output = {
        "source": "FAERS ASCII quarterly files (2022Q1–2025Q4)",
        "extraction_method": "case-insensitive drug name matching, PS-only, caseid+caseversion dedup",
        "total_ps_reports": total_ps,
        "total_pulmonary_reports": total_pulm,
        "duplicates_removed": dup_count,
        "quarter_counts": {q: c for q, c in quarter_counts},
        "pt_counts": dict(pt_counts),
        "country_counts": dict(country_counts),
        "sex_counts": dict(sex_counts),
        "pulmonary_pts_queried": sorted(PULMONARY_PTS_LOWER),
    }

    out_path = os.path.join(OUT_DIR, "olanzapine_ascii_extraction.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False, indent=2)
    print(f"\nSaved to {out_path}")

    # Also save a compact PT-count JSON for the disproportionality pipeline
    pt_json = {
        "source": "FAERS ASCII",
        "total_ps": total_ps,
        "pt_counts": {pt: pt_counts.get(pt, 0) for pt in PULMONARY_PTS_LOWER},
    }
    pt_path = os.path.join(OUT_DIR, "olanzapine_ascii_pt_counts.json")
    with open(pt_path, "w", encoding="utf-8") as f:
        json.dump(pt_json, f, ensure_ascii=False, indent=2)
    print(f"PT counts saved to {pt_path}")

    print("\nDone!")
    return output


if __name__ == "__main__":
    main()
