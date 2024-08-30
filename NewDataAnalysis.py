#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 13:45:15 2021

@author: tingtingzhao
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 16:17:08 2021

@author: tingtingzhao
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:04:57 2021

@author: tingtingzhao
"""

import os
import sys
wd = '/Users/tingtingzhao/Documents/Research/YKnock/Deep-YKnock_Tingting'
#wd
#wd = '/Users/guangyu/OneDrive - University of Florida/Research/Projects/YKnock'
os.chdir(wd)


sys.path.append(wd)
sys.path.append(wd+'/codes')

## read csv in python
import pandas as pd
from sklearn import preprocessing
from knockpy.knockoffs import GaussianSampler
from knockpy.knockoff_stats import data_dependent_threshhold
import numpy as np
from stat_modelY_classification_coef import stat_modelY_classification_coef
from sklearn.model_selection import train_test_split
from rpy2.robjects.packages import importr
knockoff = importr('knockoff')
import csv
from readGCTX import readGCTX
from create_Yknockoff_para import create_Yknockoff_para
from RealDataAnalysisUtil import findSubset, scaleYCatFeatures, getMuSigma, generateMultiKnockff,scaleYCatFeatures
import pickle
from DataAnalysisUtil import obtainYkMean, getXnum, selectVar, obtainBetaMtx,obtainSummaryBeta,sortVarFreq,writeSelectGenes
from rpy2.robjects.packages import importr
coin = importr('coin')


# Read the observations in the two cluster dataset
df = pd.read_csv("/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/TwoClusterGenes.csv")
class_0_data = df[df['cluster'] == 0]
class_1_data = df[df['cluster'] == 1]

# Remove the cluster column for 0 and 1 labels
class0 = class_0_data.drop("cluster", axis=1)
# 1370, 978
class1 = class_1_data.drop("cluster", axis=1)
# 1737, 978

with open("/Users/tingtingzhao/Documents/Research/data/rubicinAndOther/readcombinedDatavorinostat.pkl", 'rb') as f:
    GeneMtx1, combined_df1 = pickle.load(f)
#with open(r'C:\Users\student\OneDrive - Bryant University\Desktop\Honors Thesis\readcombinedDatavorinostat.pkl', 'rb') as f:
#    GeneMtx1, combined_df1 = pickle.load(f)

DMSO = GeneMtx1[3107:,:]
# combine DMSO with class0 and class 1 respectively
num1 = class0.shape[0]

dmso1 = DMSO[0:num1, :]
dmso2 = DMSO[num1:, :]

dmso1 = pd.DataFrame(dmso1, columns=class0.columns)
dmso2 = pd.DataFrame(dmso2, columns=class0.columns) 


# concatenate DMSO1 with class0, and concatenate DMSO2 with class1
exp1 = pd.concat([class0, dmso1], axis=0)
exp2 = pd.concat([class1, dmso2], axis=0)

y1 = [1] * num1 + [0] * num1
y2 =[1] * (3107-num1) + [0] * (3107-num1)

nRep = 100 
mu1,Sigma1 = getMuSigma(exp1.values)
seed = np.arange(0, nRep, 1)
mu2,Sigma2 = getMuSigma(exp2.values)

YkList1 = generateMultiKnockff(Y=exp1,mu=mu1, Sigma=Sigma1, n=nRep,seed=seed) 
YkList2 = generateMultiKnockff(Y=exp2,mu=mu2, Sigma=Sigma2, n=nRep,seed=seed) 

with open('/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/YkListCluster1.pkl', 'wb') as f:
    pickle.dump(YkList1, f) 

with open('/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/YkListCluster1.pkl', 'wb') as f:
    pickle.dump(YkList2, f) 


path1 = '/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/'
path2 ='/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/'
YkMean1 = obtainYkMean(YkList1, path1)
YkMean2 = obtainYkMean(YkList2, path2)


np.savetxt("/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/KnockoffAverage1.csv", YkMean1, delimiter=",")
np.savetxt("/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/KnockoffAverage2.csv", YkMean2, delimiter=",")


def selectVarNew(YkList, y1, GeneMtx):
    nRep = len(YkList)
    Z_max_list = []
    W_list = []
    S_max01_list = []
    S_max02_list = []
    X_num = y1
    r = GeneMtx.shape[1]
    nSam = GeneMtx.shape[0]
    for i in range(nRep):
        print(i)
        Yk = YkList[i]
        scaleYfeatures = scaleYCatFeatures(GeneMtx,Yk)
        Z_max = stat_modelY_classification_coef(np.array(X_num), scaleYfeatures)
        Z_max_list.append(Z_max)
        Z=np.abs(Z_max)
        Z=Z/Z.max()
        W = Z[0:r] - Z[r:]
        W_list.append(W)
    
        tau_max01 = data_dependent_threshhold(W, fdr=0.1)
        S_max01 = np.where(W>tau_max01)[0]
        S_max01_list.append(S_max01)
        tau_max02 = data_dependent_threshhold(W, fdr=0.2)
        S_max02 = np.where(W>tau_max02)[0]
        S_max02_list.append(S_max02)
    return Z_max_list, W_list, S_max01_list, S_max02_list

Z_max_list1, W_list1, S_max01_list1, S_max02_list1 = selectVarNew(YkList1, y1, exp1)
with open('/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/fittedS1.pkl', 'wb') as f:
    pickle.dump([Z_max_list1, W_list1, S_max01_list1, S_max02_list1], f)  
   

Z_max_list2, W_list2, S_max01_list2, S_max02_list2= selectVarNew(YkList2, y2, exp2)
with open('/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/fittedS2.pkl', 'wb') as f:
    pickle.dump([Z_max_list2, W_list2, S_max01_list2, S_max02_list2], f)



betaMtx11, betaMtx21 = obtainBetaMtx(exp1.values, S_max01_list1, S_max02_list1, y1)
with open('/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/betaMtx_1.pkl', 'wb') as f:
    pickle.dump([betaMtx11, betaMtx21], f) 


betaMtx12, betaMtx22 = obtainBetaMtx(exp2.values, S_max01_list2, S_max02_list2, y2)
with open('/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/betaMtx2.pkl', 'wb') as f:
    pickle.dump([betaMtx12, betaMtx22], f) 




r = GeneMtx1.shape[1]
summaryBeta11, resultSummary11, nonZeroRows11 = obtainSummaryBeta(r, betaMtx11)
summaryBeta21, resultSummary21, nonZeroRows21 = obtainSummaryBeta(r, betaMtx21)    
gene_info_path = "/Users/tingtingzhao/Documents/Research/data/GSE92742_Broad_LINCS_gene_info.txt"
outputPath1 = "/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/Cluster1"
sorted_dict011, sorted_dict021 = sortVarFreq(S_max01_list1, S_max02_list1, nonZeroRows11, nonZeroRows21)   
gene01All, gene02All, df_gene1, df_gene2, gene01Top80, gene02Top80 = writeSelectGenes(gene_info_path, GeneMtx1, sorted_dict011, sorted_dict021, 
                     summaryBeta11, summaryBeta21, combined_df1, outputPath1)
    

r = GeneMtx1.shape[1]
summaryBeta12, resultSummary12, nonZeroRows12 = obtainSummaryBeta(r, betaMtx12)
summaryBeta22, resultSummary22, nonZeroRows22 = obtainSummaryBeta(r, betaMtx22)    
outputPath2 = "/Users/tingtingzhao/Documents/Teaching/HornorsThesis/HonorThesisKatie/Cluster2"
sorted_dict012, sorted_dict022 = sortVarFreq(S_max01_list2, S_max02_list2, nonZeroRows12, nonZeroRows22)   
gene012All, gene022All, df_gene12, df_gene22, gene012Top80, gene022Top80 = writeSelectGenes(gene_info_path, GeneMtx1, sorted_dict012, sorted_dict022, 
                     summaryBeta12, summaryBeta22, combined_df1, outputPath2)


















