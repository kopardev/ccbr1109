import pandas
from functools import reduce

#parameters
file_dir = '/data/CCBR/projects/ccbr1109/rnaseq/'
output_dir = '/data/CCBR/projects/ccbr1109/analysis/deg_data/'
date_in = '20210811'
DEG_list = ['HNPK_SH1-SCR','HNPK_SH2-SCR','GCN1_SH2-SCR']

#read in file
def readdf(filename):
    data=pandas.read_csv(filename,sep="\t",header=0,usecols=["ensid_gene","fc","fdr","gsea_ranking_score"])
    new = data["ensid_gene"].str.split("|", n = 0, expand = True)
    data['ensid']=new[0]
    data.drop(columns=['ensid_gene'],inplace=True)
    data=data[["ensid","fc","fdr","gsea_ranking_score"]]
    return data

#create df with fc and fdr only
def merge_fc(type_in):
    flag=1

    for comps in DEG_list:
        tmp = readdf(file_dir + "DEG_" + comps + "_1_1/" + type_in + "_DEG_" + comps + "_all_genes.txt")
        tmp.drop(columns=['gsea_ranking_score'],inplace=True)
        tmp.columns = ["gene",comps + ".fc", comps + ".fdr"]

        if flag==1:
            mergeddf = tmp
            flag=2
        else:
            #merge df's
            mergeddf=reduce(lambda a,b:pandas.merge(a,b,how="outer",on="gene"),[tmp,mergeddf])   

    #clean merge df's
    mergeddf.fillna(' ',inplace=True)
    mergeddf.drop_duplicates(inplace=True)
    mergeddf.to_csv(output_dir + "merged_deg_fc_" + type_in + "_" + date_in + ".txt",sep="\t",index=False)

#create df with fc fdr and gsea
def merge_gsea(type_in):
    flag = 1

    for comps in DEG_list:
        tmp = readdf(file_dir + "DEG_" + comps + "_1_1/" + type_in + "_DEG_" + comps + "_all_genes.txt")
        tmp.columns = ["gene", comps + ".fc", comps + ".fdr", comps + ".gsea"]

        if flag==1:
            mergeddf = tmp
            flag=2
        else:
            #merge df's
            mergeddf=reduce(lambda a,b:pandas.merge(a,b,how="outer",on="gene"),[tmp,mergeddf])   

    #clean merge df's
    mergeddf.fillna(' ',inplace=True)
    mergeddf.drop_duplicates(inplace=True)
    mergeddf.to_csv(output_dir + "merged_deg_gsea_" + type_in + "_" + date_in + ".txt",sep="\t",index=False)

#generate df for deg analysis
merge_fc("DESeq2")
merge_fc("limma")
    
#create gsea df
merge_gsea("DESeq2")
merge_gsea("limma")