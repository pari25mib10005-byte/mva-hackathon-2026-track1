import subprocess, sys, numpy as np
VCF = sys.argv[1]
cmd = ["bcftools","query","-r","1,20,21,22",
       "-i",'TYPE="snp" && GT="het" && FMT/DP>=30 && FMT/GQ>=30 && QUAL>=50',
       "-f","%CHROM\t[%AD]\n",VCF]
data = {}
for line in subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True).stdout:
    c, ad = line.rstrip().split("\t"); a = ad.split(",")
    if len(a)!=2: continue
    r, al = int(a[0]), int(a[1]); data.setdefault(c,[]).append(al/(r+al))
bins = np.arange(0.30, 0.71, 0.02)
for c in ["1","20","21","22"]:
    b = np.array(data[c]); h,_ = np.histogram(b, bins); h = h/h.max()
    print(f"\nchr{c}  n={len(b):,}  mean={b.mean():.3f}")
    for lo,v in zip(bins[:-1],h): print(f"  {lo:.2f} {'#'*int(v*50)}")
