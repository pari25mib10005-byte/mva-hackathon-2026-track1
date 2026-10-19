#!/usr/bin/env python3
"""
Track 2, evidence line 2 — phenotype-driven analysis (Team PHA).

Input: the patient's HPO terms in a LOCAL text file (one HP:xxxxxxx per line;
lines ending in "# family" are family history and are excluded from the
patient's phenotype). The file never needs to leave your machine.

Steps
  1. Score every gene by how well it explains the patient's phenotype
     (sum of information content of matched HPO terms; specific terms count more).
  2. Validation: report where the known causal gene ranks (phenotype only).
  3. Coverage: which patient terms the causal gene(s) explain, and which they
     do not (flag for clinicians: second process, prematurity, or under-annotation).
  4. Network: STRING (score >= 0.7) among the top phenotype genes; the causal
     gene's direct partners form the "shared biology" module.
  5. Pathways: Reactome / KEGG / GO enrichment of the module (STRING API).
  6. Drugs: Open Targets + DGIdb on module genes, filtered with the SAME safety
     and direction rules as repurpose.py (config/bub1b.json).

Usage:
    python phenotype.py --hpo my_patient_hpo.txt --gene BUB1B --related CEP57,TRIP13 \
                        --config ../pipeline/config/bub1b.json --out results/
"""
import argparse, collections, json, math, os, sys, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
import repurpose as rp  # reuse database clients and filters

HPO_P2G = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype_to_genes.txt"
STRING = "https://string-db.org/api/json/"


def fetch(url, path):
    if not os.path.exists(path):
        print(f"      downloading {os.path.basename(path)} …")
        urllib.request.urlretrieve(url, path)
    return path


def string_post(endpoint, ids):
    data = urllib.parse.urlencode({"identifiers": "\r".join(ids), "species": 9606,
                                   "required_score": 700, "caller_identity": "teampha"}).encode()
    with urllib.request.urlopen(urllib.request.Request(STRING + endpoint, data=data), timeout=120) as r:
        return json.load(r)


def hyper_sf(k, N, K, n):
    """P(X >= k), X ~ Hypergeometric(N, K, n)."""
    lc = lambda a, b: math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    if k <= 0:
        return 1.0
    return min(1.0, sum(math.exp(lc(K, i) + lc(N - K, n - i) - lc(N, n)) for i in range(k, min(K, n) + 1)))


def read_hpo(path):
    own, family = [], []
    for line in open(path):
        raw = line.strip()
        if not raw or raw.startswith("#"):
            continue
        hp = raw.split()[0].upper().replace("HP", "HP:").replace("HP::", "HP:")
        (family if "# family" in raw.lower() else own).append(hp)
    return sorted(set(own)), sorted(set(family))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hpo", required=True)
    ap.add_argument("--gene", required=True, help="causal gene from Track 1")
    ap.add_argument("--related", default="", help="comma-separated genes causing the same disease")
    ap.add_argument("--config", required=True)
    ap.add_argument("--top", type=int, default=200)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    gene = a.gene.upper()
    related = [g.strip().upper() for g in a.related.split(",") if g.strip()]

    with open(a.config) as f:
        cfg = dict(rp.DEFAULT_CONFIG, **json.load(f))
    if cfg.get("extends"):
        with open(os.path.join(os.path.dirname(a.config), cfg["extends"])) as f:
            base = json.load(f)
        for k in ("exclude_targets", "exclude_mechanism_keywords"):
            cfg[k] = sorted(set(base.get(k, []) + cfg.get(k, [])))

    print("[1/6] Reading HPO terms and annotations")
    own, family = read_hpo(a.hpo)
    p2g = fetch(HPO_P2G, os.path.join(HERE, "data", "phenotype_to_genes.txt"))
    term_genes, names, all_genes = collections.defaultdict(set), {}, set()
    for line in open(p2g):
        f = line.rstrip("\n").split("\t")
        if f[0] == "hpo_id":
            continue
        all_genes.add(f[3])
        if f[0] in own or f[0] in family:
            term_genes[f[0]].add(f[3]); names[f[0]] = f[1]
    missing = [t for t in own if t not in term_genes]
    own = [t for t in own if t in term_genes]
    N = len(all_genes)
    ic = {t: math.log(N / len(term_genes[t])) for t in own}

    print("[2/6] Scoring genes by phenotype match")
    score, hits = collections.Counter(), collections.defaultdict(list)
    for t in own:
        for g in term_genes[t]:
            score[g] += ic[t]; hits[g].append(t)
    ranked = [g for g, _ in score.most_common()]
    rank_of = {g: ranked.index(g) + 1 for g in [gene] + related if g in ranked}

    print("[3/6] Coverage of the patient's phenotype")
    disease_genes = [gene] + related
    coverage = {t: sorted(g for g in disease_genes if g in term_genes[t]) for t in own}

    print("[4/6] Network module (STRING >= 0.7)")
    top = ranked[:a.top]
    for g in disease_genes:
        if g not in top:
            top.append(g)
    edges = string_post("network", top)
    adj = collections.defaultdict(set)
    for e in edges:
        adj[e["preferredName_A"]].add(e["preferredName_B"]); adj[e["preferredName_B"]].add(e["preferredName_A"])
    module = sorted(adj[gene] | {gene} | (set(related) & set(score)))
    second = sorted(set().union(*[adj[x] for x in adj[gene]]) - set(module)) if adj[gene] else []

    print("[5/6] Pathway enrichment of the module")
    enr = [x for x in string_post("enrichment", module)
           if x["category"] in ("RCTM", "KEGG", "Process") and x["fdr"] < 1e-3]
    enr.sort(key=lambda x: x["fdr"])

    print("      feature test: are BUB1B's own STRING partners enriched for each patient feature?")
    q = urllib.parse.urlencode({"identifiers": gene, "species": 9606, "required_score": 700, "limit": 100,
                                "caller_identity": "teampha"})
    with urllib.request.urlopen(STRING + "interaction_partners?" + q, timeout=120) as r:
        partners = {x["preferredName_B"] for x in json.load(r)} & all_genes
    feature_rows = []
    for t in own:
        K, k = len(term_genes[t] & all_genes), len(term_genes[t] & partners)
        pval = hyper_sf(k, N, K, len(partners))
        if coverage[t]:
            cat = "genetic: disease gene annotated" + (" + network-enriched" if pval < 0.01 else "")
        elif pval < 0.01:
            cat = "genetic: network-enriched (disease-gene partners)"
        else:
            cat = "not explained by the disease-gene network: secondary, acquired, or a second cause"
        feature_rows.append((t, k, len(partners) * K / N, pval, cat, sorted(term_genes[t] & partners)))

    print("[6/6] Drugs acting on the module (with safety + direction filters)")
    dgi = rp.dgidb_drugs(module)
    kept, excluded, unclear = {}, [], set()
    for g in module:
        try:
            ensg = rp.ot_resolve(g)
        except SystemExit:
            continue
        for rec in rp.ot_drugs(ensg) + dgi.get(g, []):
            reason = rp.exclusion_reason(rec, cfg)
            ok, note = rp.direction_check(g, rec, cfg, gene)
            if reason or not ok:
                excluded.append((rec["drug"], g, reason or f"wrong direction: {note}"))
                continue
            if "unknown" in note:   # association only (e.g. TP53 as a trial biomarker): not evidence of action
                unclear.add((rec["drug"], g, rec.get("source", "")))
                continue
            kept.setdefault(rec["drug"], {"drug": rec["drug"], "gene": g, "stage": rec.get("stage") or "",
                                          "direction": note, "symptoms": "; ".join(names[t] for t in hits.get(g, []))})
    profiles = rp.ot_drug_profiles(list(kept))
    for name in list(kept):
        p = profiles.get(name)
        if p and rp.exclusion_reason({**p, "actions": []}, cfg):
            excluded.append((name, kept[name]["gene"], "whole-drug check: " + rp.exclusion_reason({**p, "actions": []}, cfg)))
            del kept[name]

    # ---------------- outputs
    with open(os.path.join(a.out, "gene_scores.csv"), "w") as f:
        f.write("rank,gene,score,n_terms,terms\n")
        for i, g in enumerate(ranked[:500], 1):
            f.write(f"{i},{g},{score[g]:.2f},{len(hits[g])},\"{'; '.join(names[t] for t in hits[g])}\"\n")
    with open(os.path.join(a.out, "module_drugs.csv"), "w") as f:
        f.write("drug,via_gene,stage,direction,symptoms_linked_to_gene\n")
        for r in kept.values():
            f.write(f"{r['drug']},{r['gene']},{r['stage']},{r['direction']},\"{r['symptoms']}\"\n")
    with open(os.path.join(a.out, "module_excluded.csv"), "w") as f:
        f.write("drug,via_gene,reason\n")
        for d, g, r in sorted(set(excluded)):
            f.write(f"{d},{g},\"{r}\"\n")

    with open(os.path.join(a.out, "module_unclear_direction.csv"), "w") as f:
        f.write("drug,via_gene,source\n")
        for d, g, src in sorted(unclear):
            f.write(f"{d},{g},{src}\n")

    md = [f"# Phenotype-driven analysis: {gene}\n", "> Hypotheses for investigation, not treatment advice.\n",
          f"- Patient terms used: {len(own)} (family-history terms excluded: {len(family)}; unknown IDs: {len(missing)})",
          f"- Genes with any HPO annotation: {N}; genes matching ≥1 patient term: {len(score)}\n",
          "## 1. Validation: does the phenotype alone point to the causal gene?\n"]
    for g, r in rank_of.items():
        md.append(f"- **{g}** ranks **{r} of {len(score)}** from phenotype alone (top {100*r/len(score):.1f}%)")
    md += ["\n## 2. Phenotype coverage by the disease gene(s)\n", "| Patient term | Explained by |", "|---|---|"]
    for t in own:
        md.append(f"| {names[t]} ({t}) | {', '.join(coverage[t]) or '**not explained — flag for clinicians**'} |")
    md += ["\n### Genetic or not: network test per feature\n",
           f"{gene}'s {len(partners)} strongest STRING partners (score ≥ 0.7) with HPO annotations, tested per feature "
           f"(hypergeometric, universe {N} genes):\n",
           "| Feature | Partners linked | Expected | p | Category |", "|---|---|---|---|---|"]
    for t, k, exp, pval, cat, gs in feature_rows:
        md.append(f"| {names[t]} ({t}) | {k} | {exp:.1f} | {pval:.2g} | {cat} |")
    if family:
        md.append("\nFamily-history terms (not used for the patient's phenotype): " +
                  ", ".join(f"{names.get(t, t)} ({t})" for t in family))
    md += ["\n## 3. Shared biology: the causal gene's network module\n",
           f"Direct STRING partners (≥0.7) of {gene} among the top {a.top} phenotype genes: **{', '.join(sorted(adj[gene]))}**",
           f"\nSecond-step neighbours (context only, not proposed as targets): {', '.join(second[:40])}",
           "\n| Pathway / process | Source | Genes | FDR |", "|---|---|---|---|"]
    for x in enr[:12]:
        md.append(f"| {x['description']} | {x['category']} | {x['number_of_genes']} | {x['fdr']:.1e} |")
    md += [f"\n## 4. Drugs acting on the module\n",
           f"- Retained (known action, correct direction, passed safety): **{len(kept)}**",
           f"- Association-only links with no stated action (e.g. a gene used as a trial biomarker), set aside: **{len(unclear)}**",
           f"- Excluded: **{len(set(excluded))}**\n"]
    for r in kept.values():
        md.append(f"- {r['drug']} via {r['gene']} ({r['stage'] or 'stage n/a'}; {r['direction']})")
    reasons = collections.Counter(r.split(":")[0].split(" (")[0] for _, _, r in set(excluded))
    md.append("\nExclusion reasons: " + "; ".join(f"{v} × {k}" for k, v in reasons.most_common()))
    with open(os.path.join(a.out, "phenotype_summary.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"Done → {a.out}/  (module {len(module)} genes, {len(kept)} kept, {len(unclear)} association-only, {len(set(excluded))} excluded)")


if __name__ == "__main__":
    main()
