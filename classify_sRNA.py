#!/usr/bin/env python3
"""
sRNA classification script
Usage: python3 classify_sRNA.py <with_id.bed> <inter_mirna.tsv> <inter_te.tsv> <inter_tas.tsv> <output.bed> <cutoff_21> <cutoff_24> <min_support>
"""
import sys
import re
from collections import Counter, defaultdict

with_id_file    = sys.argv[1]
inter_mirna_f   = sys.argv[2]
inter_te_f      = sys.argv[3]
inter_tas_f     = sys.argv[4]
out_file        = sys.argv[5]
phase_cutoff_21 = float(sys.argv[6])
phase_cutoff_24 = float(sys.argv[7])
min_phase_support = int(sys.argv[8])

NL = chr(10)


def safe_float(x):
    try:
        if x is None:
            return None
        x = str(x).strip()
        if x in ("", ".", "NA", "NaN", "nan", "None", "N"):
            return None
        return float(x)
    except Exception:
        return None


def normalize_dicer(x):
    if x is None:
        return "N"
    s = str(x).strip()
    if s in ("", ".", "NA", "NaN", "nan", "None", "N", "NoCall"):
        return "N"
    m = re.search(r"(20|21|22|23|24|25|26)", s)
    return m.group(1) if m else "N"


def is_mirna_shortstack(x):
    if x is None:
        return False
    s = str(x).strip()
    if s in ("", ".", "N", "NA", "NaN", "nan", "None", "No", "False", "FALSE", "false", "0"):
        return False
    return True


def choose_mode(values):
    valid = [x for x in values if x not in ("", ".", "N", "NA", "None", None)]
    if not valid:
        return "N", "N:0"
    counter = Counter(valid)
    sorted_items = sorted(
        counter.items(),
        key=lambda item: (-item[1], int(item[0]) if item[0].isdigit() else 999999)
    )
    call = sorted_items[0][0]
    summary = ",".join(f"{k}:{v}" for k, v in sorted_items)
    return call, summary


# ==== Read merged records ====
records = {}
with open(with_id_file) as fh:
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 8:
            continue
        if len(cols) < 9:
            cols.append(".")

        cid = cols[3]
        samples = cols[4].split(",") if cols[4] else []
        dicers  = cols[5].split(",") if cols[5] else []
        phases  = cols[6].split(",") if cols[6] else []
        mirnas  = cols[7].split(",") if cols[7] else []
        strands = cols[8].split(",") if cols[8] else []

        n = max(len(samples), len(dicers), len(phases), len(mirnas), len(strands))
        while len(samples) < n:
            samples.append("")
        while len(dicers) < n:
            dicers.append("N")
        while len(phases) < n:
            phases.append("0")
        while len(mirnas) < n:
            mirnas.append("N")
        while len(strands) < n:
            strands.append(".")

        records[cid] = {
            "chrom": cols[0], "start": cols[1], "end": cols[2],
            "samples": samples, "dicers": dicers, "phases": phases,
            "mirnas": mirnas, "strands": strands,
        }


# ==== Read miRNA intersect ====
mirna_hits = defaultdict(set)
with open(inter_mirna_f) as fh:
    for line in fh:
        cols = line.rstrip().split("\t")
        if len(cols) < 5:
            continue
        cid = cols[3]
        if cols[4] == ".":
            continue
        attr = cols[-1]
        m = (re.search(r"Name=([^;]+)", attr) or
             re.search(r"name=([^;]+)", attr) or
             re.search(r"ID=([^;]+)", attr))
        mirna_hits[cid].add(m.group(1).strip() if m else "known_miRNA")


# ==== Read TE intersect ====
te_hits = set()
with open(inter_te_f) as fh:
    for line in fh:
        cols = line.rstrip().split("\t")
        if len(cols) < 5:
            continue
        cid = cols[3]
        if cols[4] != ".":
            te_hits.add(cid)


# ==== Read TAS intersect ====
tas_hits = set()
with open(inter_tas_f) as fh:
    for line in fh:
        cols = line.rstrip().split("\t")
        if len(cols) < 5:
            continue
        cid = cols[3]
        if cols[4] != ".":
            tas_hits.add(cid)


# ==== Classification function ====
def classify_locus(cid, info):
    samples = info["samples"]
    dicers  = info["dicers"]
    phases  = info["phases"]
    mirnas  = info["mirnas"]

    sample_to_dicers = defaultdict(list)
    sample_to_phases = defaultdict(list)
    sample_to_mirna  = defaultdict(bool)

    for s, d, p, m in zip(samples, dicers, phases, mirnas):
        s = s.strip()
        if s in ("", ".", "NA", "None"):
            continue
        d = normalize_dicer(d)
        pv = safe_float(p) or 0.0
        sample_to_dicers[s].append(d)
        sample_to_phases[s].append(pv)
        if is_mirna_shortstack(m):
            sample_to_mirna[s] = True

    sample_list = sorted(sample_to_phases.keys())
    sample_str = ",".join(sample_list) if sample_list else "."
    sample_count = len(sample_list)

    per_sample_dicer = []
    per_sample_phase = []
    shortstack_mirna = False

    for s in sample_list:
        d_mode, _ = choose_mode(sample_to_dicers[s])
        per_sample_dicer.append(d_mode)
        phase_max_s = max(sample_to_phases[s]) if sample_to_phases[s] else 0.0
        per_sample_phase.append(phase_max_s)
        if sample_to_mirna[s]:
            shortstack_mirna = True

    dicer_call, dicer_summary = choose_mode(per_sample_dicer)

    if per_sample_phase:
        phase_max  = max(per_sample_phase)
        phase_mean = sum(per_sample_phase) / len(per_sample_phase)
    else:
        phase_max = phase_mean = 0.0

    if dicer_call == "24":
        cutoff = phase_cutoff_24
    else:
        cutoff = phase_cutoff_21

    phase_support = sum(1 for x in per_sample_phase if x >= cutoff)

    ref_names = sorted(mirna_hits.get(cid, []))
    in_te  = cid in te_hits
    in_tas = cid in tas_hits
    is_phased = (phase_max >= cutoff and phase_support >= min_phase_support)

    # ========== Classification decision tree ==========
    # Priority:
    #   1. known miRNA      - overlaps reference miRNA GFF3
    #   2. tasiRNA          - 21/22-nt + phased + TAS locus
    #   3. 21-nt_phasiRNA   - 21/22-nt + phased
    #   4. 24-nt_phasiRNA   - 24-nt + phased (strict threshold)
    #   5. easiRNA          - 21/22-nt + TE region + not phased
    #   6. hc-siRNA         - 24-nt + not phased
    #   7. novel_miRNA      - ShortStack MIRNA=Y but none of above
    #   8. 21nt/22nt/23nt-siRNA - by DicerCall
    #   9. other-sRNA       - unclassified

    if ref_names:
        srna_type = "miRNA"
    elif in_tas and dicer_call in ("21", "22") and is_phased:
        srna_type = "tasiRNA"
    elif dicer_call in ("21", "22") and is_phased:
        srna_type = "21-nt_phasiRNA"
    elif dicer_call == "24" and is_phased:
        srna_type = "24-nt_phasiRNA"
    elif dicer_call in ("21", "22") and in_te:
        srna_type = "easiRNA"
    elif dicer_call == "24":
        srna_type = "hc-siRNA"
    elif shortstack_mirna:
        srna_type = "novel_miRNA"
    elif dicer_call == "21":
        srna_type = "21nt-siRNA"
    elif dicer_call == "22":
        srna_type = "22nt-siRNA"
    elif dicer_call == "23":
        srna_type = "23nt-siRNA"
    else:
        srna_type = "other-sRNA"

    if ref_names:
        srna_name = ",".join(ref_names)
    else:
        srna_name = cid

    return {
        "srna_name": srna_name,
        "srna_type": srna_type,
        "sample": sample_str,
        "sample_count": sample_count,
        "dicer_call": dicer_call,
        "dicer_summary": dicer_summary,
        "phase_max": phase_max,
        "phase_mean": phase_mean,
        "phase_support": phase_support,
        "in_te": in_te,
        "in_tas": in_tas,
    }


# ==== Write output ====
type_count = Counter()

with open(out_file, "w") as out:
    header = "\t".join([
        "#chrom", "start", "end", "sRNA_name", "sRNA_type",
        "sample", "sample_count", "DicerCall", "DicerCall_summary",
        "PhaseScore_max", "PhaseScore_mean", "PhaseScore_support",
        "overlap_TE", "overlap_TAS",
    ])
    out.write(header + NL)

    for cid in records:
        info = records[cid]
        result = classify_locus(cid, info)
        type_count[result["srna_type"]] += 1

        row = "\t".join([
            info["chrom"], info["start"], info["end"],
            result["srna_name"], result["srna_type"],
            result["sample"], str(result["sample_count"]),
            result["dicer_call"], result["dicer_summary"],
            f"{result['phase_max']:.4f}", f"{result['phase_mean']:.4f}",
            str(result["phase_support"]),
            "Y" if result["in_te"] else "N",
            "Y" if result["in_tas"] else "N",
        ])
        out.write(row + NL)

print(f"  Written: {out_file}")
print(f"  Total loci: {sum(type_count.values())}")
print("  Type distribution:")
for t, n in sorted(type_count.items(), key=lambda x: (-x[1], x[0])):
    print(f"    {t:<20} {n}")
