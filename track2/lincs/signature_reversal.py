#!/usr/bin/env python3
"""
Track 2, evidence line 3 — LINCS L1000 signature reversal (Team PHA).

Idea: find drugs whose effect on human cells is the MIRROR IMAGE of losing the
disease gene, with no assumption about any pathway.

Data (public, downloaded automatically from Enrichr / LINCS L1000):
  - LINCS_L1000_CRISPR_KO_Consensus_Sigs : genes up/down after CRISPR knockout
  - LINCS_L1000_Chem_Pert_Consensus_Sigs : genes up/down after each compound

Score for each compound d against knockout signature K:
  reverse = |K_up ∩ d_down| + |K_down ∩ d_up|
  mimic   = |K_up ∩ d_up|   + |K_down ∩ d_down|
  score   = (reverse - mimic) / (|K_up| + |K_down|)     (positive = reverses)
  p       = hypergeometric P(reverse overlap >= observed) ; BH-FDR across compounds

Pre-specified test (fixed before looking at results): do NAD+/sirtuin-axis
compounds (config "axis_compounds") rank higher than all other compounds?
One-sided Mann–Whitney U on reversal score.

Then: safety filter (same rules as repurpose.py) and clinical stage for the
top reversers, via Open Targets.

Usage:
    python signature_reversal.py --genes BUB1B,CEP57,TRIP13 \
        --config ../pipeline/config/bub1b.json --out results/
"""
import argparse, collections, json, math, os, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
import repurpose as rp

ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="
KO_LIB, CHEM_LIB = "LINCS_L1000_CRISPR_KO_Consensus_Sigs", "LINCS_L1000_Chem_Pert_Consensus_Sigs"
AXIS_DEFAULT = ["Nicotinamide", "Nadide", "SRT-1720", "Resveratrol", "SRT-2104", "Nicotinamide riboside",
                "beta-Nicotinamide mononucleotide"]


def load_lib(name):
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    path = os.path.join(HERE, "data", name + ".txt")
    if not os.path.exists(path):
        print(f"      downloading {name} …")
        urllib.request.urlretrieve(ENRICHR + name, path)
    sigs = collections.defaultdict(dict)
    for line in open(path, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        label = f[0].strip()
        if label.endswith(" Up"):
            sigs[label[:-3]]["up"] = {g for g in f[2:] if g}
        elif label.endswith(" Down"):
            sigs[label[:-5]]["down"] = {g for g in f[2:] if g}
    return {k: v for k, v in sigs.items() if "up" in v and "down" in v}


def log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hyper_sf(k, N, K, n):
    """P(X >= k), X ~ Hypergeometric(N, K, n)."""
    if k <= 0:
        return 1.0
    top = min(K, n)
    if k > top:
        return 0.0
    lc = log_comb(N, n)
    return min(1.0, sum(math.exp(log_comb(K, i) + log_comb(N - K, n - i) - lc) for i in range(k, top + 1)))


def binom_sf(k, n):
    """P(X >= k), X ~ Binomial(n, 0.5)."""
    if n == 0 or k <= 0:
        return 1.0
    return min(1.0, sum(math.exp(log_comb(n, i) - n * math.log(2)) for i in range(k, n + 1)))


def bh(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, pvals[i] * n / rank)
        q[i] = prev
    return q


def mann_whitney_greater(x, y):
    """One-sided U test (x > y), normal approximation with tie correction."""
    allv = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    ranks, i = [0.0] * len(allv), 0
    ties = 0.0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for t in range(i, j + 1):
            ranks[t] = r
        c = j - i + 1
        ties += c ** 3 - c
        i = j + 1
    n1, n2 = len(x), len(y)
    r1 = sum(r for r, (_, g) in zip(ranks, allv) if g == 0)
    u = r1 - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    n = n1 + n2
    sd = math.sqrt(n1 * n2 / 12 * ((n + 1) - ties / (n * (n - 1))))
    z = (u - mu - 0.5) / sd if sd else 0.0
    return u, z, 0.5 * math.erfc(z / math.sqrt(2))


def atc_l01_names():
    """WHO ATC class L01 (antineoplastic agents) from ChEMBL: excluded by name so a
    chemotherapy cannot slip through when a database lacks its mechanism."""
    path = os.path.join(HERE, "data", "atc_L01.json")
    if not os.path.exists(path):
        urllib.request.urlretrieve("https://www.ebi.ac.uk/chembl/api/data/atc_class.json?level2=L01&limit=1000", path)
    return {norm(x["who_name"]): x["level3_description"].lower() for x in json.load(open(path))["atc"]}


def norm(x):
    return "".join(ch for ch in x.upper() if ch.isalnum())


def ot_profiles(names, batch=15):
    """Mechanism text, targets and max clinical stage per LINCS compound name.
    LINCS spellings differ from Open Targets (hyphens, salts), so match on
    normalised names and accept a salt/ester suffix (e.g. 'Estramustine-Phosphate')."""
    out = {}
    frag = ('d{j}: search(queryString:{q}, entityNames:["drug"], page:{{index:0,size:3}}){{'
            ' hits{{ name object{{ ... on Drug {{ name maximumClinicalStage mechanismsOfAction{{ rows{{'
            ' mechanismOfAction targets{{ approvedSymbol }} }} }} }} }} }} }}')
    for i in range(0, len(names), batch):
        chunk = names[i:i + batch]
        qs = [json.dumps(n.replace("-", " ")) for n in chunk]
        d = rp.gql(rp.OT_URL, "{ " + " ".join(frag.format(j=j, q=q) for j, q in enumerate(qs)) + " }")
        for j, n in enumerate(chunk):
            for h in d[f"d{j}"]["hits"]:
                hn, nn = norm(h["name"]), norm(n)
                if hn == nn or hn.startswith(nn) or nn.startswith(hn):
                    o = h.get("object") or {}
                    rows = (o.get("mechanismsOfAction") or {}).get("rows") or []
                    out[n] = {"moa_text": "; ".join(r["mechanismOfAction"] for r in rows if r.get("mechanismOfAction")),
                              "moa_targets": sorted({t["approvedSymbol"] for r in rows for t in (r.get("targets") or [])}),
                              "stage": o.get("maximumClinicalStage") or ""}
                    break
    return out


def ot_stage(names, batch=20):
    """Max clinical stage per compound name via Open Targets search."""
    out = {}
    frag = ('d{j}: search(queryString:{q}, entityNames:["drug"], page:{{index:0,size:1}}){{'
            ' hits{{ name object{{ ... on Drug {{ maximumClinicalStage }} }} }} }}')
    for i in range(0, len(names), batch):
        chunk = names[i:i + batch]
        d = rp.gql(rp.OT_URL, "{ " + " ".join(frag.format(j=j, q=json.dumps(n)) for j, n in enumerate(chunk)) + " }")
        for j, n in enumerate(chunk):
            h = d[f"d{j}"]["hits"]
            if h and h[0]["name"].upper() == n.upper() and h[0].get("object"):
                out[n] = h[0]["object"].get("maximumClinicalStage") or ""
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, help="disease gene(s) with knockout signatures, comma-separated")
    ap.add_argument("--config", required=True)
    ap.add_argument("--axis", default=",".join(AXIS_DEFAULT), help="pre-specified compounds to test as a group")
    ap.add_argument("--top", type=int, default=100)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    with open(a.config) as f:
        cfg = dict(rp.DEFAULT_CONFIG, **json.load(f))
    if cfg.get("extends"):
        with open(os.path.join(os.path.dirname(a.config), cfg["extends"])) as f:
            base = json.load(f)
        for k in ("exclude_targets", "exclude_mechanism_keywords"):
            cfg[k] = sorted(set(base.get(k, []) + cfg.get(k, [])))

    print("[1/5] Loading LINCS consensus signatures")
    ko, chem = load_lib(KO_LIB), load_lib(CHEM_LIB)
    genes = [g.strip().upper() for g in a.genes.split(",") if g.strip()]
    found = [g for g in genes if g in ko]
    if not found:
        sys.exit(f"No knockout signature for {genes}")
    universe = set()
    for s in list(ko.values()) + list(chem.values()):
        universe |= s["up"] | s["down"]
    N = len(universe)
    print(f"      knockout signatures used: {', '.join(found)} (missing: {', '.join(set(genes) - set(found)) or 'none'})")
    print(f"      compounds: {len(chem)}; gene universe: {N}")

    lc = cfg.get("lincs", {})
    alarm = set(lc.get("alarm_genes", [])) & ko[found[0]]["up"]
    alarm_min = lc.get("alarm_suppression_min", 3)
    patterns = [[w.lower() for w in pat] for pat in cfg.get("exclude_patterns", []) + lc.get("exclude_patterns", [])]
    known = {**{k: v for k, v in cfg.get("known_mechanisms", {}).items() if not k.startswith("_")},
             **lc.get("known_mechanisms", {})}
    known = {**known, **{k.upper(): v for k, v in known.items()}}
    print(f"      alarm genes present in the knockout 'up' half: {len(alarm)}")
    print("[2/5] Scoring every compound against each knockout signature")
    rows = []
    for name, d in chem.items():
        per = {}
        for g in found:
            K = ko[g]
            sup, res = len(K["up"] & d["down"]), len(K["down"] & d["up"])
            rev, mim = sup + res, len(K["up"] & d["up"]) + len(K["down"] & d["down"])
            size = len(K["up"]) + len(K["down"])
            p = binom_sf(rev, rev + mim)          # direction test: are overlaps mostly reversing?
            p_res = hyper_sf(res, N, len(K["down"]), len(d["up"]))   # restore arm alone
            per[g] = ((rev - mim) / size, rev, mim, p, sup, res, p_res)
        main = per[found[0]]
        alarm_off = len(alarm & d["down"])
        rows.append({"compound": name, "score": main[0], "reverse": main[1], "mimic": main[2], "p": main[3],
                     "suppress_arm": main[4], "restore_arm": main[5], "p_restore": main[6],
                     "alarm_genes_switched_off": alarm_off, "alarm_suppressor": alarm_off >= alarm_min,
                     "consistent": all(v[0] > 0 for v in per.values()),
                     "mean_score_all": sum(v[0] for v in per.values()) / len(per),
                     **{f"score_{g}": round(v[0], 4) for g, v in per.items()}})
    for r, q in zip(rows, bh([r["p"] for r in rows])):
        r["fdr"] = q
    for r, q in zip(rows, bh([r["p_restore"] for r in rows])):
        r["fdr_restore"] = q
    rows.sort(key=lambda r: (-r["score"], r["p"]))
    for i, r in enumerate(rows, 1):
        r["rank"] = i

    print("[3/5] Pre-specified test: NAD+/sirtuin-axis compounds vs all others")
    axis_names = [x.strip() for x in a.axis.split(",") if x.strip()]
    axis = [r for r in rows if r["compound"].lower() in {x.lower() for x in axis_names}]
    others = [r["score"] for r in rows if r not in axis]
    u, z, p_axis = mann_whitney_greater([r["score"] for r in axis], others) if axis else (0, 0, 1)

    print("[4/5] Safety filter + clinical stage for the top reversers")
    top = [r for r in rows if r["score"] > 0 and r["fdr"] < 0.05][:a.top]
    restore = sorted([r for r in rows if r["fdr_restore"] < 0.05 and not r["alarm_suppressor"]],
                     key=lambda r: r["p_restore"])[:a.top]
    review = {r["compound"]: r for r in top + restore}
    prof = ot_profiles(list(review))
    l01 = atc_l01_names()
    for name, r in review.items():
        p = prof.get(name) or {}
        r["mechanism"] = p.get("moa_text") or known.get(name, "") or known.get(name.upper(), "")
        r["stage"] = p.get("stage", "")
        reason = rp.exclusion_reason({"moa_text": r["mechanism"], "moa_targets": p.get("moa_targets", []),
                                      "actions": []}, cfg)
        text = r["mechanism"].lower()
        nn = norm(name)
        hit = next((cls for an, cls in l01.items() if an and (nn == an or nn.startswith(an) or an.startswith(nn))), None)
        if not reason and hit:
            reason = f"antineoplastic (WHO ATC L01: {hit})"
        if not reason:
            for pat in patterns:
                if all(w in text for w in pat):
                    reason = "wrong direction or unsafe class: " + " ".join(pat)
                    break
        if not reason and r["alarm_suppressor"]:
            reason = f"switches off {r['alarm_genes_switched_off']} p53/interferon alarm genes"
        r["excluded"] = reason or ""

    print("[5/5] Writing outputs")
    cols = ["rank", "compound", "score", "reverse", "mimic", "p", "fdr", "suppress_arm", "restore_arm",
            "p_restore", "fdr_restore", "alarm_genes_switched_off", "consistent", "mean_score_all"] + \
           [f"score_{g}" for g in found]
    with open(os.path.join(a.out, "reversal_scores.csv"), "w", encoding="utf-8") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(f"\"{r[c]}\"" if c == "compound" else
                             (f"{r[c]:.3g}" if isinstance(r[c], float) else str(r[c])) for c in cols) + "\n")
    with open(os.path.join(a.out, "top_reversers_reviewed.csv"), "w", encoding="utf-8") as f:
        f.write("rank,compound,score,fdr,stage,excluded,mechanism\n")
        for r in top:
            f.write(f"{r['rank']},\"{r['compound']}\",{r['score']:.3f},{r['fdr']:.2g},{r['stage']},"
                    f"\"{r['excluded'] or ''}\",\"{r['mechanism']}\"\n")

    restore_safe = [r for r in restore if not r["excluded"]]
    with open(os.path.join(a.out, "restore_candidates.csv"), "w", encoding="utf-8") as f:
        f.write("compound,restore_arm,p_restore,fdr_restore,alarm_genes_switched_off,stage,excluded,mechanism\n")
        for r in restore:
            f.write(f"\"{r['compound']}\",{r['restore_arm']},{r['p_restore']:.2g},{r['fdr_restore']:.2g},"
                    f"{r['alarm_genes_switched_off']},{r['stage']},\"{r['excluded']}\",\"{r['mechanism']}\"\n")
    safe = [r for r in top if not r["excluded"]]
    approved_safe = [r for r in safe if r["stage"] in ("APPROVAL", "PHASE_4")]
    K = ko[found[0]]
    md = [f"# LINCS signature reversal: {', '.join(found)}\n", "> Hypotheses for investigation, not treatment advice.\n",
          f"- Knockout signature ({found[0]}): {len(K['up'])} genes up, {len(K['down'])} genes down",
          f"- Compounds scored: {len(rows)}; gene universe: {N}",
          f"- Significant reversers (score > 0, FDR < 0.05): {sum(1 for r in rows if r['score'] > 0 and r['fdr'] < 0.05)}",
          f"- Knockout signatures available in LINCS: {', '.join(found)}"
          + (f" (not available: {', '.join(g for g in genes if g not in found)})" if len(found) < len(genes) else "")
          + (f"; reversers consistent across all of them: "
             f"{sum(1 for r in rows if r['consistent'] and r['score'] > 0 and r['fdr'] < 0.05)}" if len(found) > 1 else "")
          + "\n",
          "## Pre-specified test: NAD⁺/sirtuin axis\n",
          "| Compound | Rank (of %d) | Score | FDR |" % len(rows), "|---|---|---|---|"]
    for r in axis:
        md.append(f"| {r['compound']} | {r['rank']} | {r['score']:.3f} | {r['fdr']:.2g} |")
    md += [f"\nOne-sided Mann–Whitney (axis > others): U = {u:.0f}, z = {z:.2f}, **p = {p_axis:.3g}** "
           f"(n = {len(axis)} axis compounds found of {len(axis_names)} pre-specified).\n",
           "## Top reversers after safety review\n",
           f"- Reviewed: {len(top)}; excluded by safety rules: {len(top) - len(safe)}; "
           f"remaining: {len(safe)}; of which approved: {len(approved_safe)}\n",
           "| Rank | Compound | Score | FDR | Stage | Mechanism |", "|---|---|---|---|---|---|"]
    for r in safe[:20]:
        md.append(f"| {r['rank']} | {r['compound']} | {r['score']:.3f} | {r['fdr']:.2g} | {r['stage'] or 'n/a'} | "
                  f"{(r['mechanism'] or 'n/a')[:70]} |")
    n_alarm = sum(1 for r in top if r["alarm_suppressor"])
    md += [f"\n**Alarm check:** {n_alarm} of the top {len(top)} reversers switch off ≥{alarm_min} p53/interferon alarm genes "
           f"({len(alarm)} alarm genes are in the knockout 'up' half).\n",
           "## Exploratory: restore division genes without silencing the alarm\n",
           f"Compounds that raise the knockout's 'down' genes (FDR < 0.05) and switch off < {alarm_min} alarm genes; "
           f"after safety review: {len(restore_safe)} of {len(restore)}; approved: "
           f"{sum(1 for r in restore_safe if r['stage'] in ('APPROVAL', 'PHASE_4'))}.\n",
           "| Compound | Division genes restored | FDR | Stage | Mechanism |", "|---|---|---|---|---|"]
    for r in [r for r in restore_safe if r["stage"] in ("APPROVAL", "PHASE_4")][:15]:
        md.append(f"| {r['compound']} | {r['restore_arm']} | {r['fdr_restore']:.2g} | {r['stage']} | {(r['mechanism'] or 'n/a')[:60]} |")
    def bucket(x):
        x = x.split(":")[0].split(" (")[0]
        return "switches off the p53/interferon alarm" if x.startswith("switches off") else x
    reasons = collections.Counter(bucket(r["excluded"]) for r in top if r["excluded"])
    md.append("\nExcluded among top reversers: " + ("; ".join(f"{v} × {k}" for k, v in reasons.most_common()) or "none"))
    with open(os.path.join(a.out, "lincs_summary.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    json.dump({"axis_test": {"U": u, "z": z, "p": p_axis, "n_axis": len(axis)},
               "knockouts": found, "n_compounds": len(rows), "universe": N},
              open(os.path.join(a.out, "lincs_provenance.json"), "w"), indent=2)
    print(f"Done → {a.out}/  (axis test p = {p_axis:.3g}; {len(safe)} safe top reversers, {len(approved_safe)} approved)")


if __name__ == "__main__":
    main()
