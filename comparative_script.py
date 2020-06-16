#!/usr/bin/env python
#Author: Yunlong Ma
#Email: glb-biotech@zju.edu.cn
#Usage: Code for comparative analysis which is a part of an analysis of integrative genomics analysis of Endometriosis (EM)


import numpy


#Part I Load data

#Read data 
f1 = open("Geneset1.txt","r+")
f2 = open("Geneset2.txt","r+")
f3 = open("Geneset13.txt","r+")
f4 = open("Geneset4.txt","r+")
f5 = open("Geneset5.txt","r+")


#Part II Extract genes according to different P values

#Establish a function for extracting subgroup genes at three different P value thresholds
#P values: 0.05, 0.01, and 0.001 as 3 comparative points from Sherlock analysis
def Subgroup_extract(x):
    Geneset_1 =[]
    Geneset_2 =[]
    Geneset_3 =[]
    for line in x:
        if line[0:3] !="Gene":
            dd = line.strip().split()
            Geneset_1.append(dd[0])
            if dd[1]!="P" and float(dd[1]) < 0.01:
                Geneset_2.append(dd[0])
            if dd[1]!="P" and float(dd[1]) < 0.001:
                Geneset_3.append(dd[0])
    
    return(Geneset_1, Geneset_2, Geneset_3)

         
#Generate subgroup of genes for Genesets #1, #2, and #3
result_1 = Subgroup_extract(f1)
result_2 = Subgroup_extract(f2)
result_3 = Subgroup_extract(f3)


#f3 MAGMA analysis of GWAS on EM (MAGMA-based P value <0.05)
Geneset4 =[]
for i in f4:
    if i[0:3]!="Gene":
        dd3 = i.strip()
        Geneset4.append(dd3)
        

#f4 MAGMA analysis of GWAS on null trait (MAGMA-based P value <0.05)
Geneset5=[]
for j in f5:
    if j[0:3]!="Gene":
        dd4 = j.strip()
        Geneset5.append(dd4)


#Part III Calculate overlapped gene rates at 3 thresholds of P values 


#Establish a novel function for calculate the overlapped gene rates
def fun_rate(x,y,z):
    sum1 = len(x)
    Combined = []
    for num in x:
        if num in y:
            Combined.append(num)
    sum2 = len(Combined)
    rate = sum2/sum1
    
    Combined_null = []
    for num2 in x:
        if num2 in z:
            Combined_null.append(num2)
    sum_null_1 = len(Combined_null)
    if sum_null_1 == 0:
        sum_null_1 = 0.0001
    rate2 = sum_null_1/sum1
    return(rate,rate2)
 
         
#For Sherlock analysis of Zeller blood-based eQTL vs. MAGMA  
#At the threshold of P value = 0.05
Rate1 = fun_rate(result_1[0],Geneset4,Geneset5)
#At the threshold of P value = 0.01
Rate2 = fun_rate(result_1[1],Geneset4,Geneset5)
#At the threshold of P value = 0.001
Rate3 = fun_rate(result_1[2],Geneset4,Geneset5)
EM_group = [Rate1[0],Rate2[0],Rate3[0]]
Null_group =  [Rate1[1],Rate2[1],Rate3[1]]

#For Sherlock analysis of Dixon blood-based eQTL vs. MAGMA
#At the threshold of P value = 0.05
Rate4 = fun_rate(result_2[0],Geneset4,Geneset5)
#At the threshold of P value = 0.01
Rate5 = fun_rate(result_2[1],Geneset4,Geneset5)
#At the threshold of P value = 0.001
Rate6 = fun_rate(result_2[2],Geneset4,Geneset5)
EM_group_2 = [Rate4[0],Rate5[0],Rate6[0]]
Null_group_2 =  [Rate4[1],Rate5[1],Rate6[1]]


#For Sherlock analysis of GTEx blood-based eQTL vs. MAGMA
#At the threshold of P value = 0.05
Rate7 = fun_rate(result_3[0], Geneset4,Geneset5)
#At the threshold of P value = 0.01
Rate8 = fun_rate(result_3[1],Geneset4,Geneset5)
#At the threshold of P value = 0.001
Rate9 = fun_rate(result_3[2],Geneset4,Geneset5)
EM_group_3 = [Rate7[0],Rate8[0],Rate9[0]]
Null_group_3 =  [Rate7[1],Rate8[1],Rate9[1]]



if __name__ == "__main__":
    print("This a comparative analysis")

#End


