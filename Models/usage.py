import pandas as pd
import numpy as np
import re
import pickle
import sklearn
from collections import Counter
from Bio import motifs
from Bio.Seq import Seq
from Bio.motifs.matrix import PositionWeightMatrix,PositionSpecificScoringMatrix
import warnings

warnings.filterwarnings("ignore")


### load models
with open('RFs_kmersFinal.pkl', 'rb') as handle:
    RFs_kmers = pickle.load(handle)

with open('RFs_flanksFinal.pkl', 'rb') as handle:
    RFs_flanks = pickle.load(handle)

with open('RFs_metaFinal.pkl', 'rb') as handle:
    meta_rfc = pickle.load(handle)



## columns of meta model
colsmeta=list()
with open("meta_colsFinal.txt","r") as file:
    for line in file:
        colsmeta.append(line.strip())

## columns of flanks model
colsflanks=list()
with open("Flanks_colsFinal.txt","r") as file:
    for line in file:
        colsflanks.append(line.strip())

## columns of kmers model
colskmers=list()
with open("Kmers_colsFinal.txt","r") as file:
    for line in file:
        colskmers.append(line.strip())

with open('FinalRBps.pkl', 'rb') as f:
    rbps = pickle.load(f)


for i in rbps: ## convert all dataframes to numeric
    rbps[i]=rbps[i].apply(pd.to_numeric)

for i in rbps:
    data=rbps[i]
    newdata=dict()
    newdata['A']=list(data.iloc[0,:])
    newdata['C']=list(data.iloc[1,:])
    newdata['G']=list(data.iloc[2,:])
    newdata['T']=list(data.iloc[3,:])
    rbps[i]=PositionSpecificScoringMatrix("ACGT",newdata)



#### Save the efficiencies
nearcogEff=pd.read_csv("NearCognateEff.txt",sep="\t")
atgEff=pd.read_csv("ATGEff.txt",sep="\t")
ncEff=dict()
atEff=dict()
totnearcogeff=sum(nearcogEff["TIS Efficiency"])
totalatgeff=sum(atgEff["efficiency"])
for i in range(nearcogEff.shape[0]):
    seq=nearcogEff.loc[i,"TIS Sequence"].replace("U","T")
    ncEff[seq]=nearcogEff.loc[i,"TIS Efficiency"]/totnearcogeff

for i in range(atgEff.shape[0]):
    seq=atgEff.loc[i,"sequence"].replace("U","T")
    atEff[seq]=atgEff.loc[i,"efficiency"]/totalatgeff

def count_kmers(sequence, k):
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return Counter(kmers)

def normalize_kmer_counts(sequence, k):
    kmer_counts = count_kmers(sequence, k)
    total_kmers = len(sequence) - k + 1
    normalized_counts = {kmer: count / total_kmers for kmer, count in kmer_counts.items()}
    return normalized_counts

def substring_occurance(string,cod):
    '''Get occurances of start codons in cdna'''
    arr=[m.start() for m in re.finditer(cod, string)]
    return arr

def getFeatures(seq,codons,rbps,atgeffs,nceffs):
    '''Retrieve required features'''
    datapos=list()
    kmer3_order=['CTC', 'TCA', 'CAC', 'ACA', 'CAG', 'AGT', 'GTG', 'TGA', 'GAC', 'ACC', 'CCC', 'CCT', 'CTG', 'GAT', 'ATC', 'TCT', 'TGG', 'GGT', 'GTA', 'TAA', 'AAA', 'AAG', 'AGC', 'GCT', 'TCC', 'CCA', 'CAT', 'TGC', 'GCC', 'ACT', 'TGT', 'GTC', 'ATG', 'GGG', 'GGC', 'GCA', 'AGG', 'AGA', 'GGA', 'GAG', 'AAC', 'CAA', 'GAA', 'TAG', 'AAT', 'ATT', 'TTT', 'TTC', 'ATA', 'TAC', 'CTT', 'TTA', 'CCG', 'CGA', 'GTT', 'CTA', 'ACG', 'CGG', 'GCG', 'TCG', 'TTG', 'CGT', 'CGC', 'TAT']
    kmer2_order=['CT','TG','GC','CG','CC','CA','AT','TC','GG','AG','GT','TT','AC','GA','AA','TA']
    kmer1_order=['C','T','G','A']
    for startcodon in codons:
        positions=substring_occurance(seq,startcodon)
        for z in positions:
            ## skip the positions at the right edge
            if z>len(seq)-20 or z<20:
                continue
            arr=list()
            ## upstream sequences
            k3u=normalize_kmer_counts(seq[0:z],3)
            k2u=normalize_kmer_counts(seq[0:z],2)
            k1u=normalize_kmer_counts(seq[0:z],1)
            ### use ordering to keep it
            for j in kmer3_order:
                if j in k3u:
                    arr.append(k3u[j])
                else:
                    arr.append(0)
            for j in kmer2_order:
                if j in k2u:
                    arr.append(k2u[j])
                else:
                    arr.append(0)
            for j in kmer1_order:
                if j in k1u:
                    arr.append(k1u[j])
                else:
                    arr.append(0)
            ## downstream sequences
            k3d=normalize_kmer_counts(seq[z+3:len(seq)],3)
            k2d=normalize_kmer_counts(seq[z+3:len(seq)],2)
            k1d=normalize_kmer_counts(seq[z+3:len(seq)],1)
            ### use ordering to keep it
            for j in kmer3_order:
                if j in k3d:
                    arr.append(k3d[j])
                else:
                    arr.append(0)
            for j in kmer2_order:
                if j in k2d:
                    arr.append(k2d[j])
                else:
                    arr.append(0)
            for j in kmer1_order:
                if j in k1d:
                    arr.append(k1d[j])
                else:
                    arr.append(0)
            ## append start codon
            startcod=seq[z:z+3]
            arr.append(startcod)
            ## append efficiencies
            if startcod=="ATG":
                arr.append(atgeffs[seq[z-6:z+5]])
            else:
                arr.append(nceffs[seq[z-4:z+4]])
            
            temp=list(seq[z-20:z])
            arr.extend(temp)
            temp=list(seq[z+3:z+23])
            arr.extend(temp)
            ## add rbp info
            avg = []
            for rbp_id in rbps:
                # Precompute data and k-mers for the current RBP
                data = rbps[rbp_id]
                seqs = np.array(data.calculate(seq[0:z]))
                seqs=seqs/data.calculate(data.consensus)
                temp = []
                temp=[z-i for i, value in enumerate(seqs) if value >= 0.8]
                # Calculate average distance if matches are found
                if temp:
                    avg.append(np.mean(temp))
                else:
                    avg.append(0)
            arr.extend(avg)
            arr.append(z)
            aa=translate(seq[z:len(seq)])
            if ("_3Extension" in aa) | ("NoStop" in aa):
                arr.append("No Stop Codon")
                arr.append("No Stop Codon")
            else:
                arr.append(z+(len(aa)*3)-1)
                arr.append(len(aa))
            arr.append(aa)
            datapos.append(arr)
        datapos=pd.DataFrame(datapos)
        colnames=list()
        orders=list()
        orders.extend(kmer3_order)
        orders.extend(kmer2_order)
        orders.extend(kmer1_order)
        for i in range(2):
            for j in orders:
                if i==0:
                    colnames.append(j+"_U")
                else:
                    colnames.append(j+"_D")
        colnames.append("start_codon")
        colnames.append("efficiency")
        ## add flank column names
        ## upstream
        for i in range(20,0,-1): 
            colnames.append(f'{i}_U')
        ## downstream
        for i in range(20):
            colnames.append(f'{i+1}_D')
        colnames.extend(list(rbps.keys()))
        colnames.append("Start")
        colnames.append("Stop")
        colnames.append("ORF Length")
        colnames.append("Amino Acid")
        datapos.columns=colnames
        if datapos.shape[0]==0:
            return None
        else:
            return datapos


def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""

    for i in range(0, len(seq), 3):

        codon = seq[i:i + 3]
        if len(codon)<3:
            protein+="_3Extension"
            return protein

        if table[codon]=="_":
            return protein
        protein += table[codon]
    protein+="NoStop"
    return protein


def metaPred(metamodel,flanks,kmers,colsmeta,colsflanks,colskmers,data):
    '''Predict TIS'''
    X_comparison=data.iloc[:,0:data.shape[1]-4]
    comp=pd.get_dummies(X_comparison)
    for i in colsmeta:
        if i in comp.columns:
            continue
        else:
            comp[i]=0
    #### KmersERF
    kmerscomp=comp[colskmers]
    pred_probas=list()
    for i in kmers:
        arr=kmers[i].predict_proba(kmerscomp)[:,1]
        pred_probas.append(list(arr))

    pred_probas=pd.DataFrame(pred_probas)
    meanpreds_kmers=pred_probas.apply(np.mean)
    finalpreds_kmers=list()
    for i in meanpreds_kmers:
        if i>=0.5:
            finalpreds_kmers.append(1)
        else:
            finalpreds_kmers.append(0)

    ### FlanksERF
    flankscomp=comp[colsflanks]
    pred_probas=list()
    for i in flanks:
        arr=flanks[i].predict_proba(flankscomp)[:,1]
        pred_probas.append(list(arr))

    pred_probas=pd.DataFrame(pred_probas)
    meanpreds_flanks=pred_probas.apply(np.mean)
    finalpreds_flanks=list()
    for i in meanpreds_flanks:
        if i>=0.5:
            finalpreds_flanks.append(1)
        else:
            finalpreds_flanks.append(0)

    ## meta model 
    comp["kmers_pred"]=meanpreds_kmers
    comp["flanks_pred"]=meanpreds_flanks
    X_comparison=comp[colsmeta]

    finalpredsproba=metamodel.predict_proba(X_comparison)[:,1]
    #finalpreds=meta_rfc.predict(X_comparison)
    data2=data.iloc[:,data.shape[1]-4:data.shape[1]]
    data2["KmersERF"]=meanpreds_kmers
    data2["FlanksERF"]=meanpreds_flanks
    data2["MetaTIS"]=finalpredsproba
    data2["Start Codon"]=data["start_codon"]
    finalcols=["Start","Stop","ORF Length","Amino Acid","Start Codon","KmersERF","FlanksERF","MetaTIS"]
    data2=data2[finalcols]
    return data2



