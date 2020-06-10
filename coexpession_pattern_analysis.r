#Rscript
#Author: Yunlong Ma
#Co-expression patterns by using the Pearson correlation analysis
#Usage: it is designed to uncover the co-expression patterns among genes associated with Endometriosis (EM)

#Set up the working directory depended on user own need.
setwd("C:\\Users\\Administrator\\Desktop\\corplot")

#Part I install package

#Installing corrplot package of R platform 
if(!require("corrplot"))install.packages("corrplot")
if(!require("corrr"))install.packages("corrr")
if(!require("dplyr"))install.packages("dplyr")
if(!require("ggplot2"))install.packages("ggplot2")
if(!require("reshape2"))install.packages("reshape2")


#Load the corrplot package
library(corrplot)
library(corrr)
library(dplyr)
library(ggplot2)
library(reshape2)

set.seed(123456)


#Part II read relevant RNA expression data of co-expression analysis according to disease status


#Established a function for reading and re-organization expression data
cor_processing <- function(x){
  pearson<- read.delim(x,header = T)
  pearson_new2<-pearson[,c(-1,-2)] 
  mat_pearson<-as.matrix(pearson_new2)
  row.names(mat_pearson)<- pearson[,2]
  data_for_x <-t(mat_pearson)
  correlation_x <- cor(data_for_x )
  return(correlation_x)
}


#Part III Visualization

#calculating correlations among genes based on expression data
cor_statistics_EM <- cor_processing("EM.txt")
cor_statistics_Con <- cor_processing("Control.txt")
cor_statistics_Con_fertile <- cor_processing("Control_fertile.txt")
cor_statistics_Con_infertile <- cor_processing("Control_infertile.txt")
cor_statistics_EM_fertile <- cor_processing("Endometriosis_fertile.txt")
cor_statistics_EM_infertile <- cor_processing("Endometriosis_infertile.txt")


#make corrplot for each group
corrplot(cor_statistics_EM, method = "circle",tl.col = "black",tl.cex=0.9,number.cex = 1.8)
corrplot(cor_statistics_Con, method = "circle",tl.col = "black",tl.cex=0.9,number.cex = 1.8)
corrplot(cor_statistics_Con_fertile, method = "circle",tl.col = "black",tl.cex=0.9,number.cex = 1.8)
corrplot(cor_statistics_Con_infertile, method = "circle",tl.col = "black",tl.cex=0.9,number.cex = 1.8)
corrplot(cor_statistics_EM_fertile, method = "circle",tl.col = "black",tl.cex=0.9,number.cex = 1.8)
corrplot(cor_statistics_EM_infertile, method = "circle",tl.col = "black",tl.cex=0.9,number.cex = 1.8)


#Part IV Use corrr package to visualize gene-gene co-expression patterns in a network


#Established a function for reading and re-organization expression data
data_processing <- function(x){
  tt<- read.delim(x,header = T)
  dd<-tt[,c(-1,-2)] 
  mat_cor<-as.matrix(dd)
  row.names(mat_cor)<- tt[,2]
  data_for_x <-t(mat_cor)
  return(data_for_x)
}

#data_processing
data_for_EM <- data_processing("EM.txt")
data_for_Con <- data_processing("Control.txt")
data_for_Con_fertile <- data_processing("Control_fertile.txt")
data_for_Con_infertile <- data_processing("Control_infertile.txt")
data_for_EM_fertile <- data_processing("Endometriosis_fertile.txt")
data_for_EM_infertile <- data_processing("Endometriosis_infertile.txt")


#Visualization for co-expression network analysis
#Set the correlation efficient >=0.3 to show
data_for_EM %>% correlate() %>%  rearrange() %>% network_plot(min_cor = 0.3) 
data_for_Con %>% correlate() %>%  rearrange() %>% network_plot(min_cor = 0.3)
data_for_Con_fertile %>% correlate() %>%  rearrange() %>% network_plot(min_cor = 0.3)
data_for_Con_infertile %>% correlate() %>%  rearrange() %>% network_plot(min_cor = 0.3)
data_for_EM_fertile %>% correlate() %>%  rearrange() %>% network_plot(min_cor = 0.3)
data_for_EM_infertile %>% correlate() %>%  rearrange() %>% network_plot(min_cor = 0.3)


#Part V Paired Student's t-test

#Establish a function for generating correlation patterns
cor_pattern <- function(x){
  
  cor_pattern_data <- x %>% correlate() %>%  rearrange() %>% stretch()
  
  return(cor_pattern_data)
}


# Create correlation patterns
temp_EM <- cor_pattern(data_for_EM)
temp_Con <- cor_pattern(data_for_Con)
temp_Con_fertile <- cor_pattern(data_for_Con_fertile)
temp_Con_infertile <- cor_pattern(data_for_Con_infertile)
temp_EM_fertile <- cor_pattern(data_for_EM_fertile)
temp_EM_infertile <- cor_pattern(data_for_EM_infertile)


#define a function for match the correlation data betwen CoA and control
Cor_reshape <- function(x,y){
  len = length(x[,1])
  x[,4]<-y[,3]
  for (i in seq(1,len)){
    for(j in seq(1,len)){
      if (x[i,1]==y[j,1] & x[i,2]==y[j,2]){
        x[i,4] <- y[j,3]
      }
    }
  }
  colnames(x) <- c("Gene1","Gene2","Case","Control")
  final_data <- na.omit(x)
  return(final_data)
}

#use the Cor_reshape for data organization
data <- Cor_reshape(temp_EM,temp_Con)
data1 <- Cor_reshape(temp_EM_fertile,temp_Con_fertile)
data2 <- Cor_reshape(temp_EM_infertile,temp_Con_infertile)


#Establish a function for Paired t-test and density plot
data_Analyzing <- function(x){
  tt <- t.test(x$Case, x$Control, paired=TRUE)
  Pvalue <- tt[3]
  return(Pvalue)
}

P<-data_Analyzing(data)
P1<-data_Analyzing(data1)
P2<-data_Analyzing(data2)


#Plotting the density figure for data
md <- melt(data,id=c("Gene1","Gene2"))
p<-ggplot(md, aes(x = value))
p + geom_density(color = "black", fill = "gray")
p + geom_density(aes(color =variable)) 
density_plot<- p + geom_density(aes(fill = variable), alpha=0.6)

#Plotting the density figure for data1
md <- melt(data1,id=c("Gene1","Gene2"))
p<-ggplot(md, aes(x = value))
p + geom_density(color = "black", fill = "gray")
p + geom_density(aes(color =variable)) 
density_plot1<- p + geom_density(aes(fill = variable), alpha=0.6)

#Plotting the density figure for data2
md <- melt(data2,id=c("Gene1","Gene2"))
p<-ggplot(md, aes(x = value))
p + geom_density(color = "black", fill = "gray")
p + geom_density(aes(color =variable)) 
density_plot2<- p + geom_density(aes(fill = variable), alpha=0.6)



#End

