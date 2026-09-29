import subprocess, sys, numpy as np, pandas as pd
VCF = sys.argv[1]; BIN = 5_000_000
cmd = ["bcftools","query","-r","1,20,21,22",
       "-i",'TYPE="snp" && GT="het" && FMT/DP>=20 && FMT/GQ>=30 && QUAL>=50',
       "-f","%CHROM\t%POS\t[%AD]\t[%DP]\n",VCF]
rows=[]
for line in subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True).stdout:
    c,p,ad,dp=line.rstrip().split("\t"); a=ad.split(",")
    if len(a)!=2: continue
    r,al=int(a[0]),int(a[1]); rows.append((c,int(p)//BIN*BIN//1_000_000,al/(r+al),int(dp)))
df=pd.DataFrame(rows,columns=["chrom","bin_Mb","baf","dp"])
g=df.groupby(["chrom","bin_Mb"]).agg(n=("baf","size"),mean_baf=("baf","mean"),
    dev=("baf",lambda x:(x-0.5).abs().mean()),med_dp=("dp","median")).reset_index()
g=g[g.n>=500]; g["dp_ratio"]=g.med_dp/df[df.chrom=="1"].dp.median()
pd.set_option("display.width",160)
for c in ["1","20","21","22"]:
    print(f"\nchr{c}"); print(g[g.chrom==c].drop(columns="chrom").round(3).to_string(index=False))
g.round(4).to_csv("mosaic_bins.csv",index=False)
