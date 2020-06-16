#Rscript
#Author: Yunlong Ma
#E-mail: glb-biotech@zju.edu.cn

#100,000 times computer-based permutation analysis
#This script is designed for the integrative genomics analysis of endometriosis (EM)
#These significant gene sets were identified from both Sherlock-based Bayesian analysis and MAGMA gene-level analysis


#Set the work directory
setwd("C:\\Users\\Administrator\\Desktop\\05-Computer-based permutation analysis")
set.seed(12345)

#Part I Read data on significant genes and background genes

#1
#Read significant genes of Geneset #1
Sig_1 <- read.table("geneset1.txt", header=T)
Sig_geneset1 <- Sig_1$Gene_name

#2
#Read significant genes from Geneset #2
Sig_2 <- read.table("geneset2.txt", header=T)
Sig_geneset2 <- Sig_2$Gene_name

#Read background genes of Dixon blood eQTL data
Backgroud_2<- read.table("geneset2_all.txt", header=T)
Backgroud_geneset2 <- Backgroud_2$Gene_name

#3
#Read significant genes from Geneset #3
Sig_3 <- read.table("geneset3.txt", header=T)
Sig_geneset3 <- Sig_3$Gene_name

#Read background genes of GTEx blood eQTL data
Backgroud_3<- read.table("geneset3_all.txt", header=T)
Backgroud_geneset3 <- Backgroud_3$Gene_name

#4
#Read significant genes from Geneset #4
Sig_4 <- read.table("geneset4.txt", header=T)
Sig_geneset4 <- Sig_4$Gene_name

#Read background genes of MAGMA analysis on EM GWAS summary dataset
Background_4<- read.table("geneset4_all.txt", header=T)
Background_geneset4 <- Background_4$Gene_name



#Calculate the numebr of genes in each gene set 
len_Sig_geneset1 <- length(Sig_geneset1)
len_Sig_geneset2 <- length(Sig_geneset2)
len_Backgroud_geneset2 <- length(Backgroud_geneset2)
len_Sig_geneset3 <- length(Sig_geneset3)
len_Background_geneset3 <- length(Background_geneset3)
len_Sig_geneset4 <- length(Sig_geneset4)
len_Background_geneset4 <- length(Background_geneset4)


#Part II establish a function for permutation analysis

#Permutation Function
Permut_analysis <- function(x,y,z){
  
  random_selected_genes <- sample(x,y)
  
  temp <- match(random_selected_genes,z)
  
  random_overlaped_genes <- na.omit(temp)
  
  num<-length(random_overlaped_genes)
  
  return(num)
  
}


#100000 times permutation analysis for each gene set 
results_1 <- replicate(100000,Permut_analysis(Backgroud_geneset2,len_Sig_geneset2,Sig_geneset1))
results_2 <- replicate(100000,Permut_analysis(Background_geneset3,len_Sig_geneset3,Sig_geneset1))
results_3 <- replicate(100000,Permut_analysis(Background_geneset4,len_Sig_geneset4,Sig_geneset1))
results_4 <- replicate(100000,Permut_analysis(Backgroud_geneset3,len_Sig_geneset3,Sig_geneset2))
results_5 <- replicate(100000,Permut_analysis(Background_geneset4,len_Sig_geneset4,Sig_geneset2))
results_6 <- replicate(100000,Permut_analysis(Background_geneset4,len_Sig_geneset4,Sig_geneset3))


#Ploting function
Fig_random <- function(x,y,z){
  
  hist(x, col="red",xlab="Counts of overlapped genes",main=NULL)
  temp1 <- match(y,z)
  Observed_overlaped_genes <- na.omit(temp1)
  Observed_gene_num <- length(Observed_overlaped_genes)
  abline(v=Observed_gene_num,col="darkblue",lty="longdash")
  P_value=length(x[x>Observed_gene_num])/length(x)
  x1= Observed_gene_num
  freq <- table(x)
  y1 = max(freq)
  text(x1,y1,P_value)
  
}


#Visulization for the results of permutation analysis
Fig_random(results_1,Sig_geneset2,Sig_geneset1)
Fig_random(results_2,Sig_geneset3,Sig_geneset1)
Fig_random(results_3,Sig_geneset4,Sig_geneset1)
Fig_random(results_4,Sig_geneset3,Sig_geneset2)
Fig_random(results_5,Sig_geneset4,Sig_geneset2)
Fig_random(results_6,Sig_geneset4,Sig_geneset3)




#End

