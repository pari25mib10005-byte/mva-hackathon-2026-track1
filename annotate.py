import subprocess, json, requests, csv, time

VCF = "candidates.norm.vcf.gz"
out = subprocess.run(["bcftools","query","-f","%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD\t%DP]\n",VCF],
                     capture_output=True,text=True).stdout.strip().split("\n")
variants = []
for line in out:
    chrom,pos,ref,alt,gt,ad,dp = line.split("\t")
    variants.append(dict(chrom=chrom,pos=int(pos),ref=ref,alt=alt,gt=gt,ad=ad,dp=dp))

url = "https://rest.ensembl.org/vep/human/region"
hdr = {"Content-Type":"application/json","Accept":"application/json"}
params = {"canonical":1,"pick":1,"af_gnomadg":1,"sift":"b","polyphen":"b","hgvs":1}

def call(batch):
    body = {"variants":[f'{v["chrom"]} {v["pos"]} . {v["ref"]} {v["alt"]} . . .' for v in batch]}
    for attempt in range(4):
        r = requests.post(url, headers=hdr, params=params, data=json.dumps(body))
        if r.status_code == 200: return r.json()
        print(f"  HTTP {r.status_code}, retry {attempt+1}..."); time.sleep(5*(attempt+1))
    raise SystemExit(f"Failed batch starting at {batch[0]['chrom']}:{batch[0]['pos']}: {r.text[:300]}")

rows = []
for i in range(0,len(variants),50):
    batch = variants[i:i+50]
    print(f"annotating {i+1}-{i+len(batch)} of {len(variants)}")
    res_list = call(batch)
    by_key = {(str(x["seq_region_name"]),x["start"]):x for x in res_list}
    for v in batch:
        res = next((x for x in res_list if x.get("input","").startswith(f'{v["chrom"]} {v["pos"]} ')), {})
        tc = (res.get("transcript_consequences") or [{}])[0]
        cv = (res.get("colocated_variants") or [{}])[0]
        rows.append({
            "chrom":v["chrom"],"pos":v["pos"],"ref":v["ref"],"alt":v["alt"],
            "gt":v["gt"],"ad":v["ad"],"dp":v["dp"],
            "gene":tc.get("gene_symbol"),"impact":tc.get("impact"),
            "consequence":res.get("most_severe_consequence"),
            "hgvsc":tc.get("hgvsc"),"hgvsp":tc.get("hgvsp"),
            "sift":tc.get("sift_prediction"),"polyphen":tc.get("polyphen_prediction"),
            "gnomad_af":cv.get("frequencies",{}).get(v["alt"],{}).get("gnomadg"),
            "clinvar":";".join(cv.get("clin_sig",[])),
            "rsid":cv.get("id",""),
        })
    time.sleep(1)

order = {"HIGH":0,"MODERATE":1,"LOW":2,"MODIFIER":3}
rows.sort(key=lambda r:(order.get(r["impact"],4), r["gnomad_af"] if r["gnomad_af"] is not None else -1))
with open("annotated.csv","w",newline="") as f:
    w = csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(f"\n{len(rows)} variants annotated -> annotated.csv\n")
for r in rows[:15]:
    print(r["gene"],r["chrom"],r["pos"],r["ref"],r["alt"],r["gt"],r["impact"],r["consequence"],r["hgvsp"],r["gnomad_af"],r["clinvar"])
