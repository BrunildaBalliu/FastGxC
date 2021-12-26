#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Script to generate simulated SNP and expression data
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% December 1st 2021
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
n_genes = 100 # number of genes
n_snps_per_gene= 100 # number of cis-SNPs per gene
n_contexts = 50 # number of contexts (e.g. tissues or cell types)
mus = c(rep(0, n_contexts-1),0) # gene expression mean in each context for genotype = 0
maf = 0.1 # genotype minor allele frequency

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

# Save SNP location file
# Data frame with 3 initial columns (name, chrom, and position) that match standard SNP map file, followed by 1 column for each context with a 0/1 indicator of whether the given SNP passed QC in that tissue. 
snp_loc=data.frame(snpid=colnames(genos), chr="chr1",pos=rep(x = seq(1,n_genes*1e9, by = 1e9), each=n_snps_per_gene), 
                   matrix(data = 1, nrow = ncol(genos), ncol = n_contexts, dimnames = list(NULL, paste0("context",1:n_contexts))))
write.table(x = snp_loc, file = paste0(data_dir,"snpsloc.txt"), quote = F, sep = '\t', row.names = F, col.names = T)

# Save gene location file
# data frame with 4 initial columns (name, chrom, and start and end position) that match standard gene map file, followed by 1 column for each context with a 0/1 indicator of whether the given gene passed QC in that context 
gene_loc=data.frame(geneid=paste0("gene",1:n_genes), 
                    chr="chr1",
                    s1=seq(1,n_genes*1e9, by = 1e9), 
                    s2=seq(1,n_genes*1e9, by = 1e9)+ 1000,
                    matrix(data = 1, nrow = n_genes, ncol = n_contexts, dimnames = list(NULL, paste0("context",1:n_contexts))))
write.table(x = gene_loc, file = paste0(data_dir,"geneloc.txt"), quote = F, sep = '\t', row.names = F, col.names = T)

# Use only one snp per gene to generate expression
genos_with_effect = genos[,seq(from = 1, to = (n_snps_per_gene*n_genes), by = n_snps_per_gene)]

# Generate expression matrix
exp_mat=expand.grid(iid=paste0("ind",1:N),context=paste0("context",1:n_contexts))


which_context=rep_len(x = 1:n_contexts, length.out = n_genes)

for(i in 1:n_genes){
  
  # Genotypic effect in each context for assumed heritability 
  hsq=rep(NA, n_contexts)  # expression heritability, i.e. proportion of gene expression variance explained by genetics, in each context
  hsq[which_context[i]]= 0.4
  hsq[-which_context[i]]= rep(0.2, n_contexts-1)
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

write.table(x = exp_mat[,-c(1:2)], file = paste0(data_dir,"expression.txt"), quote = F, sep = '\t', row.names = T, col.names = T)
