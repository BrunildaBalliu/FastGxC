#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% April 2nd 2020
#%%%%%%%%%%%%%%% Script to generate simulated SNP and expression data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Directories for input and output files 
args=commandArgs(TRUE)
work_dir = args[1]
data_dir=paste0(work_dir,'/data/') # directory with phenotype and genotype data
if(!dir.exists(data_dir)) dir.create(data_dir)

# Required libraries
library(mvtnorm)
library(reshape2)
library(mppa)
library(MatrixEQTL)
library(ggplot2)
library(data.table)
library(dplyr)

# Parameters
N = 300 # number of individuals
n_genes = 10 # number of genes
n_snps_per_gene=5
n_contexts = 50 # number of contexts (e.g. tissues or cell types)
mus = c(rep(5, n_contexts-1),0) # gene expression mean in each context for genotype = 0
maf = 0.1 # genotype minor allele frequency
hsq = c(rep(0.2, n_contexts-1),0.4) # expression heritability, i.e. proportion of gene expression variance explained by genetics, in each context

# Error variance-covariance matrix
w_corr = 0.2 # error covariance between contexts
v_e = 1 # error variance in each context
sigma = matrix(w_corr,nrow=n_contexts,ncol=n_contexts) # Error variance-covariance matrix
diag(sigma) = v_e

# Simulate genotypes 
genos = sapply(X = 1:(n_snps_per_gene*n_genes), FUN = function(X) rbinom(N, 2, maf))
colnames(genos) = paste0("snp",1:(n_snps_per_gene*n_genes))
rownames(genos) = paste0("ind",1:N)
write.table(x = t(data.table(genos, keep.rownames = T) %>% rename(snpid=rn)), file = paste0(data_dir,"SNPs.txt"), quote = F, sep = '\t', row.names = T, col.names = F)

# Save SNP locations
snp_loc=data.frame(snpid=colnames(genos), chr="chr1",pos=c(t(matrix(seq(1,n_genes*1e06, by = 1e06)) %*% 1:n_snps_per_gene)))
write.table(x = snp_loc, file = paste0(data_dir,"SNPs_loc.txt"), quote = F, sep = '\t', row.names = F, col.names = T)

# Use only one snp per gene to generate expression
genos_with_effect = genos[,seq(from = 1, to = (n_snps_per_gene*n_genes), by = 5)]

# Generate expression matrix
exp_mat=expand.grid(iid=paste0("ind",1:N),context=paste0("context",1:n_contexts))

for(i in 1:n_genes){
  
  # Genotypic effect in each context for assumed heritability 
  betas=sqrt((hsq*v_e)/((1-hsq)*var(genos_with_effect[,i]))) 
  
  # expression of gene per context without noise
  Y = matrix(0,nrow=N,ncol=n_contexts, dimnames = list(paste0("ind",1:N), paste0("context",1:n_contexts))) 
  for (j in 1:n_contexts)  Y[,j] = mus[j] + genos_with_effect[,i]*betas[j]
  
  # get noise per individual
  for(j in 1:N) Y[j,] = Y[j,] + rmvnorm(1, rep(0,n_contexts), sigma)
  
  data_mat=melt(data = data.table(Y,keep.rownames = T), id.vars = "rn") 
  colnames(data_mat) =  c("iid", "context", paste0("gene",i))
  exp_mat = merge(x = exp_mat, y = data_mat)
  
}

rownames(exp_mat) = paste(exp_mat$iid,exp_mat$context, sep = "_")

write.table(x = exp_mat[,-c(1:2)], file = paste0(data_dir,"exp_mat.txt"), quote = F, sep = '\t', row.names = T, col.names = T)
