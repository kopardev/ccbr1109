import gtfparse
import sys

if len(sys.argv)!=2:
	print("1 argument required! Need path to GTF file!")
	exit()

df = gtfparse.read_gtf(sys.argv[1])
df = df[df["gene_type"]=="protein_coding"]
genes_only = df[df["feature"]=="gene"]
gene_boundaries=dict()
for i,row in genes_only.iterrows():
	gene_name=row["gene_name"]
	if not gene_name in gene_boundaries:
		gene_boundaries[gene_name]=dict()
	if row["strand"]=="+":
		gene_boundaries[gene_name]["start"]=row["start"]
		gene_boundaries[gene_name]["end"]=row["end"]
	else:
		gene_boundaries[gene_name]["start"]=row["end"]
		gene_boundaries[gene_name]["end"]=row["start"]
out=dict()
for k,v in gene_boundaries.items():
	x=df[df["feature"]=="UTR"]
	x=x[x["gene_name"]==k]
	for i,row in x.iterrows():
		if row["strand"]=="+":
			end=row["end"]
		else:
			end=row["start"]
		if end==v["end"]:
			r=[row["seqname"],str(row["start"]-1),str(row["end"]),k,".",row["strand"]]
			if k in out:
				utrlen=row["end"]-row["start"]
				if utrlen > out[k]['utrlen']:
					out[k]['outstr']="\t".join(r)
			else:
				out[k]=dict()
				out[k]['utrlen']=row["end"]-row["start"]
				out[k]['outstr']="\t".join(r)
for k,v in out.items():
	print(v['outstr'])
