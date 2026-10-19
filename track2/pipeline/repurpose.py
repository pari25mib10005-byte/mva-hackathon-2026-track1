#!/usr/bin/env python3
"""
Track 2 — gene-agnostic drug repurposing pipeline (Team PHA).

Given a disease gene and a mechanism config, it:
  1. resolves the gene to Ensembl (Open Targets),
  2. builds a gene set = index gene + curated modifier genes + STRING neighbours,
  3. pulls drug-target evidence from DGIdb and Open Targets,
  4. applies a hard exclusion filter (harmful mechanisms) and a direction filter
     (a drug must push each modifier gene the way the mechanism needs),
  5. merges a curated literature layer for agents the databases cannot see
     (e.g. metabolic precursors),
  6. ranks candidates and writes CSV + Markdown + a provenance JSON.

Standard library only. Output is a list of HYPOTHESES FOR INVESTIGATION,
not treatment advice.

Usage:
    python repurpose.py --config config/bub1b.json --out results/
    python repurpose.py --gene TRIP13 --out results/          # no config: generic defaults
"""
import argparse, csv, datetime, json, os, sys, time, urllib.error, urllib.request

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
DGIDB_URL = "https://dgidb.org/api/graphql"

# Stage → score. Approved/marketed drugs score highest (repurposing needs them).
STAGE_SCORE = {"APPROVAL": 1.0, "PHASE_4": 1.0, "PHASE_3": 0.6, "PHASE_2_3": 0.5,
               "PHASE_2": 0.4, "PHASE_1_2": 0.3, "PHASE_1": 0.25, "EARLY_PHASE_1": 0.2,
               "PRECLINICAL": 0.1, "WITHDRAWAL": 0.0, "UNKNOWN": 0.1, None: 0.1}

UP_ACTIONS = {"AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR", "OPENER",
              "STABILISER", "PARTIAL AGONIST", "INDUCER", "POSITIVE MODULATOR"}
DOWN_ACTIONS = {"INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR", "BLOCKER",
                "DEGRADER", "INVERSE AGONIST", "ANTISENSE INHIBITOR", "RNAI INHIBITOR",
                "DISRUPTING AGENT", "NEGATIVE MODULATOR", "RELEASING AGENT"}

# Generic defaults used when no config is given (safe for chromosome-instability genes).
DEFAULT_CONFIG = {
    "mode": "loss_of_function",
    "modifiers": {},
    "network": {"max_neighbours": 20, "min_score": 0.9},
    "exclude_targets": [],
    "exclude_mechanism_keywords": [],
    "curated_candidates": [],
}


# --------------------------------------------------------------------- HTTP
def gql(url, query, variables=None, retries=3):
    body = json.dumps({"query": query, "variables": variables or {}}).encode()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=body, headers={"Content-Type": "application/json",
                                                                  "User-Agent": "teampha-track2/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.load(r)
            if data.get("errors"):
                raise RuntimeError(f"{url}: {data['errors'][0].get('message')}")
            return data["data"]
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


# ------------------------------------------------------------- Open Targets
def ot_resolve(symbol):
    d = gql(OT_URL, 'query($q:String!){ search(queryString:$q, entityNames:["target"]){ hits{ id name } } }',
            {"q": symbol})
    for h in d["search"]["hits"]:
        if h["name"].upper() == symbol.upper():
            return h["id"]
    raise SystemExit(f"Could not resolve gene symbol '{symbol}' in Open Targets")


def ot_neighbours(ensg, size, min_score):
    d = gql(OT_URL, """query($id:String!,$n:Int!){ target(ensemblId:$id){
        interactions(sourceDatabase:string, page:{index:0,size:$n}){
          rows{ score targetB{ id approvedSymbol } } } } }""", {"id": ensg, "n": size})
    rows = (d["target"]["interactions"] or {}).get("rows") or []
    return [(r["targetB"]["approvedSymbol"], r["targetB"]["id"], r["score"])
            for r in rows if r["targetB"] and r["score"] >= min_score]


def ot_drugs(ensg):
    d = gql(OT_URL, """query($id:String!){ target(ensemblId:$id){ approvedSymbol
        drugAndClinicalCandidates{ rows{ maxClinicalStage
          drug{ id name maximumClinicalStage drugType
                mechanismsOfAction{ rows{ actionType mechanismOfAction targets{ approvedSymbol } } } } } } } }""",
            {"id": ensg})
    out = []
    for r in (d["target"]["drugAndClinicalCandidates"] or {}).get("rows") or []:
        drug = r["drug"]
        moa = (drug.get("mechanismsOfAction") or {}).get("rows") or []
        out.append({
            "drug": drug["name"].upper(), "drug_id": drug["id"],
            "stage": drug.get("maximumClinicalStage") or r.get("maxClinicalStage"),
            "actions": sorted({m["actionType"] for m in moa if m.get("actionType")}),
            "moa_text": "; ".join(sorted({m["mechanismOfAction"] for m in moa if m.get("mechanismOfAction")})),
            "moa_targets": sorted({t["approvedSymbol"] for m in moa for t in (m.get("targets") or [])}),
            "source": "OpenTargets",
        })
    return out


DRUG_PROFILE_FRAGMENT = (
    'd{j}: search(queryString:{q}, entityNames:["drug"], page:{{index:0,size:1}}){{'
    ' hits{{ name object{{ ... on Drug {{ name mechanismsOfAction{{ rows{{ actionType mechanismOfAction'
    ' targets{{ approvedSymbol }} }} }} }} }} }} }}'
)


def ot_drug_profiles(names, batch=20):
    """Full mechanism profile per drug name (all targets, not just the gene it
    was found through), so a hazardous drug cannot slip in via a side record."""
    prof = {}
    for i in range(0, len(names), batch):
        chunk = names[i:i + batch]
        parts = [DRUG_PROFILE_FRAGMENT.format(j=j, q=json.dumps(n)) for j, n in enumerate(chunk)]
        d = gql(OT_URL, "{ " + " ".join(parts) + " }")
        for j, n in enumerate(chunk):
            hits = d[f"d{j}"]["hits"]
            obj = hits[0]["object"] if hits and hits[0]["name"].upper() == n.upper() else None
            if obj:
                rows = (obj.get("mechanismsOfAction") or {}).get("rows") or []
                prof[n] = {
                    "moa_text": "; ".join(r["mechanismOfAction"] for r in rows if r.get("mechanismOfAction")),
                    "moa_targets": sorted({t["approvedSymbol"] for r in rows for t in (r.get("targets") or [])}),
                }
    return prof


def _norm(x):
    return "".join(ch for ch in x.upper() if ch.isalnum())


def atc_l01():
    """WHO ATC group L01 (antineoplastic agents) via ChEMBL, cached in data/.
    Excludes chemotherapies by NAME, so one cannot slip through when a database
    record lacks a mechanism or uses a variant name (e.g. a prodrug or salt)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "atc_L01.json")
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        urllib.request.urlretrieve("https://www.ebi.ac.uk/chembl/api/data/atc_class.json?level2=L01&limit=1000", path)
    return {_norm(x["who_name"]): x["level3_description"].lower() for x in json.load(open(path))["atc"]}


def l01_reason(name, l01):
    n = _norm(name)
    for an, cls in l01.items():
        if len(an) >= 5 and (n == an or n.startswith(an)):
            return f"antineoplastic (WHO ATC L01: {cls})"
    return None


# ------------------------------------------------------------------- DGIdb
def dgidb_drugs(symbols):
    d = gql(DGIDB_URL, """query($g:[String!]){ genes(names:$g){ nodes{ name
        interactions{ drug{ name conceptId approved } interactionScore
          interactionTypes{ type directionality } sources{ sourceDbName } publications{ pmid } } } } }""",
            {"g": symbols})
    out = {}
    for node in d["genes"]["nodes"]:
        rows = []
        for i in node["interactions"]:
            types = [t["type"].upper() for t in i["interactionTypes"] if t.get("type")]
            dirs = {t.get("directionality") for t in i["interactionTypes"]}
            rows.append({
                "drug": i["drug"]["name"].upper(), "drug_id": i["drug"]["conceptId"],
                "stage": "APPROVAL" if i["drug"]["approved"] else None,
                "actions": sorted(set(types)),
                "direction_hint": "down" if "INHIBITORY" in dirs else "up" if "ACTIVATING" in dirs else None,
                "moa_text": "", "moa_targets": [node["name"]],
                "n_pubs": len(i["publications"]),
                "source": "DGIdb:" + ",".join(sorted({s["sourceDbName"] for s in i["sources"]})),
            })
        out[node["name"]] = rows
    return out


# --------------------------------------------------------------- filtering
def direction_of(rec):
    acts = set(rec.get("actions") or [])
    if acts & DOWN_ACTIONS:
        return "down"
    if acts & UP_ACTIONS:
        return "up"
    return rec.get("direction_hint")


def exclusion_reason(rec, cfg):
    ex_targets = {t.upper() for t in cfg["exclude_targets"]}
    hit = ex_targets & {t.upper() for t in rec.get("moa_targets", [])} if direction_of(rec) != "up" else set()
    if hit:
        return f"acts on excluded target(s) {','.join(sorted(hit))}"
    text = rec.get("moa_text", "").lower()
    for kw in cfg["exclude_mechanism_keywords"]:
        if kw.lower() in text:
            return f"mechanism matches excluded class '{kw}'"
    return None


def direction_check(gene, rec, cfg, index_gene):
    """Return (ok, note). For a LoF index gene, inhibiting it or its core
    complex is harmful; for modifiers the config states the wanted direction."""
    d = direction_of(rec)
    if gene in cfg["modifiers"]:
        want = cfg["modifiers"][gene]["want"]
        if d is None:
            return True, f"direction unknown (want {want})"
        return (d == want), f"{d} vs wanted {want}"
    if cfg["mode"] == "loss_of_function" and d == "down":
        if gene == index_gene:
            return False, "inhibits a loss-of-function disease gene"
        # Neighbours of a LoF checkpoint gene are mostly its mitotic partners:
        # inhibiting them risks phenocopying or worsening the defect.
        return False, "inhibits a network partner of a loss-of-function gene"
    return True, d or "direction unknown"


# ------------------------------------------------------------------ ranking
def score(rec, gene_weight, dir_ok_known):
    # Unknown direction halves the gene contribution: a binder with no stated
    # effect is weak evidence that the drug pushes the pathway the right way.
    s = 0.5 * STAGE_SCORE.get(rec.get("stage"), 0.1)
    s += 0.3 * gene_weight * (1.0 if dir_ok_known else 0.5)
    s += 0.15 if dir_ok_known else 0.0
    s += min(rec.get("n_pubs", 0), 5) * 0.01
    return round(s, 3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", help="mechanism config JSON (see config/bub1b.json)")
    ap.add_argument("--gene", help="index gene symbol (overrides config)")
    ap.add_argument("--out", default="results")
    ap.add_argument("--approved-only", action="store_true", help="keep only approved/phase-4 agents")
    args = ap.parse_args()

    cfg = dict(DEFAULT_CONFIG)
    if args.config:
        with open(args.config) as f:
            user = json.load(f)
        if user.get("extends"):  # merge a shared hazard filter (lists are concatenated)
            with open(os.path.join(os.path.dirname(args.config), user["extends"])) as f:
                base = json.load(f)
            for k in ("exclude_targets", "exclude_mechanism_keywords"):
                user[k] = sorted(set(base.get(k, []) + user.get(k, [])))
        cfg.update(user)
    gene = (args.gene or cfg.get("gene") or "").upper()
    if not gene:
        sys.exit("Give --gene or a config with 'gene'")
    os.makedirs(args.out, exist_ok=True)

    print(f"[1/6] Resolving {gene}")
    ensg = ot_resolve(gene)

    print("[2/6] Building gene set (index + modifiers + STRING neighbours)")
    genes = {gene: {"id": ensg, "weight": 1.0, "role": "index"}}
    for m, spec in cfg["modifiers"].items():
        genes[m.upper()] = {"id": ot_resolve(m), "weight": spec.get("weight", 0.8),
                            "role": f"modifier (want {spec['want']}): {spec.get('why', '')}"}
    net = cfg["network"]
    for sym, gid, sc in ot_neighbours(ensg, net["max_neighbours"], net["min_score"]):
        genes.setdefault(sym, {"id": gid, "weight": round(0.5 * sc, 3), "role": f"STRING neighbour ({sc})"})
    print(f"      {len(genes)} genes")

    print("[3/6] Querying DGIdb + Open Targets")
    dgi = dgidb_drugs(list(genes))
    records = []
    for g, info in genes.items():
        for rec in ot_drugs(info["id"]) + dgi.get(g, []):
            rec["gene"] = g
            records.append(rec)
    print(f"      {len(records)} raw drug–gene records")

    print("[4/6] Applying exclusion + direction filters")
    kept, excluded = {}, []
    known = {k.upper(): v for k, v in cfg.get("known_mechanisms", {}).items() if not k.startswith("_")}
    patterns = [[w.lower() for w in pat] for pat in cfg.get("exclude_patterns", [])]
    for rec in records:
        g = rec["gene"]
        if not rec.get("moa_text") and rec["drug"] in known:
            rec["moa_text"] = known[rec["drug"]]
        reason = exclusion_reason(rec, cfg)
        text = (rec.get("moa_text", "") + " " + rec["drug"]).lower()
        if not reason:
            for pat in patterns:
                if all(w in text for w in pat):
                    reason = "wrong direction or unsafe class: " + " ".join(pat)
                    break
        ok, note = direction_check(g, rec, cfg, gene)
        if not reason and ok and "unknown" in note and STAGE_SCORE.get(rec.get("stage"), 0) < 1.0:
            reason = "investigational compound with unknown direction (not proposable for a child)"
        if reason or not ok:
            excluded.append({"drug": rec["drug"], "gene": g, "source": rec["source"],
                             "reason": reason or f"wrong direction: {note}"})
            continue
        spec = cfg["modifiers"].get(g, {})
        if spec.get("propose") is False:
            excluded.append({"drug": rec["drug"], "gene": g, "source": rec["source"],
                             "reason": "direction consistent but not proposable: " + spec.get("no_propose_why", "")})
            continue
        if args.approved_only and STAGE_SCORE.get(rec.get("stage"), 0) < 1.0:
            continue
        s = score(rec, genes[g]["weight"], "unknown" not in note)
        key = rec["drug"]
        prev = kept.get(key)
        if prev is None or s > prev["score"]:
            kept[key] = {"drug": key, "gene": g, "role": genes[g]["role"], "stage": rec.get("stage") or "",
                         "action": ",".join(rec.get("actions") or []), "direction": note,
                         "mechanism": rec.get("moa_text", ""), "source": rec["source"],
                         "evidence_layer": "database", "score": s}

    # Whole-drug safety pass: re-check every retained drug against its full profile.
    profiles = ot_drug_profiles(list(kept))
    for name in list(kept):
        p = profiles.get(name)
        if not p:
            continue
        reason = exclusion_reason({**p, "actions": []}, cfg)
        if reason:
            excluded.append({"drug": name, "gene": kept[name]["gene"], "source": "OpenTargets drug profile",
                             "reason": "whole-drug check: " + reason})
            del kept[name]
        elif not kept[name]["mechanism"]:
            kept[name]["mechanism"] = p["moa_text"]

    # Chemotherapy-by-name pass (WHO ATC L01).
    l01 = atc_l01()
    for name in list(kept):
        reason = l01_reason(name, l01)
        if reason:
            excluded.append({"drug": name, "gene": kept[name]["gene"], "source": "WHO ATC L01 (ChEMBL)",
                             "reason": reason})
            del kept[name]

    print("[5/6] Merging curated literature layer")
    for c in cfg["curated_candidates"]:
        key = c["drug"].upper()
        row = {"drug": key, "gene": c.get("gene", ""), "role": "curated", "stage": c.get("stage", ""),
               "action": c.get("action", ""), "direction": c.get("direction", ""),
               "mechanism": c.get("mechanism", ""), "source": "; ".join(c.get("refs", [])),
               "evidence_layer": "literature (not in drug–target DBs)" if key not in kept else "database + literature",
               "score": c.get("score", 0)}
        kept[key] = row

    ranked = sorted(kept.values(), key=lambda r: -r["score"])

    print("[6/6] Writing outputs")
    cols = ["rank", "drug", "score", "gene", "role", "stage", "action", "direction",
            "mechanism", "evidence_layer", "source"]
    with open(os.path.join(args.out, f"{gene}_candidates.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for i, r in enumerate(ranked, 1):
            w.writerow({"rank": i, **r})
    with open(os.path.join(args.out, f"{gene}_excluded.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["drug", "gene", "source", "reason"])
        w.writeheader()
        w.writerows(sorted(excluded, key=lambda r: (r["reason"], r["drug"])))

    n_index = sum(1 for r in records if r["gene"] == gene)
    with open(os.path.join(args.out, f"{gene}_summary.md"), "w") as f:
        f.write(f"# Repurposing screen: {gene}\n\n")
        f.write("> Hypotheses for investigation, not treatment advice.\n\n")
        f.write(f"- Run: {datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None).isoformat(timespec='seconds')}Z\n")
        f.write(f"- Gene set: {len(genes)} genes; raw records: {len(records)}; "
                f"records on {gene} itself: {n_index}\n")
        f.write(f"- Excluded: {len(excluded)} records; retained candidates: {len(ranked)}\n\n")
        headline = [r for r in ranked if r["evidence_layer"] != "database" or "unknown" not in r["direction"]]
        f.write("## Headline candidates (curated or direction-confirmed)\n\n"
                "| # | Drug | Score | Via gene | Stage | Layer |\n|---|---|---|---|---|---|\n")
        for i, r in enumerate(headline[:15], 1):
            f.write(f"| {i} | {r['drug']} | {r['score']} | {r['gene']} | {r['stage']} | {r['evidence_layer']} |\n")
        db_known = [r for r in ranked if r["evidence_layer"] == "database" and "unknown" not in r["direction"]]
        f.write(f"\n## Database hits with a confirmed, correct direction: {len(db_known)}\n\n")
        for r in db_known[:10]:
            f.write(f"- {r['drug']} via {r['gene']} ({r['action']}, {r['stage']})\n")
        if not db_known:
            f.write("- none: no approved drug in DGIdb/Open Targets is known to push this pathway the right way; "
                    "candidates come from the curated literature layer.\n")
        unk = [r for r in ranked if r["evidence_layer"] == "database" and "unknown" in r["direction"]]
        f.write(f"\n## Approved/clinical binders with unknown direction (need manual review): {len(unk)}\n\n")
        f.write(", ".join(f"{r['drug']} ({r['gene']})" for r in unk[:25]) + ("\n" if unk else "none\n"))
        reasons = {}
        for e in excluded:
            k = e["reason"].split(" (")[0]
            reasons[k] = reasons.get(k, 0) + 1
        f.write("\n## Exclusions by reason\n\n")
        for k, v in sorted(reasons.items(), key=lambda x: -x[1]):
            f.write(f"- {v} × {k}\n")

    with open(os.path.join(args.out, f"{gene}_provenance.json"), "w") as f:
        json.dump({"gene": gene, "ensembl": ensg, "run_utc": datetime.datetime.now(datetime.timezone.utc).replace(tzinfo=None).isoformat() + "Z",
                   "sources": {"open_targets": OT_URL, "dgidb": DGIDB_URL},
                   "config": cfg, "gene_set": genes}, f, indent=2)
    print(f"Done → {args.out}/  ({len(ranked)} candidates, {len(excluded)} excluded)")


if __name__ == "__main__":
    main()
