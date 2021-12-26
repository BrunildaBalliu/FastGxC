#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Script to run TreeQTL by content
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% December 1st 2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Arguments and parameters
args=commandArgs(TRUE) 
work_dir = args[1];
data_dir=paste0(work_dir,'/data/') # directory with phenotype and genotype data
out_dir=paste0(work_dir,'/results/') # MatrixEQTL output
source(paste0(work_dir, "/functions.R"))

# Libraries 
library(TreeQTL)

# use a single thread
print(paste0("data.table getDTthreads(): ",getDTthreads()))
setDTthreads(1)
print(paste0("after setting as 1 thread; data.table getDTthreads(): ",getDTthreads()))

# Display all warnings as they occur
options(warn=1)

# FDR thresholds
level1=0.05
level2=0.05
level3=0.05

# Gene and SNP positions
snpspos = read.table(file = snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(file = gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# Use treeQTL to perform hierarchical FDR and get sp_eGenes, i.e. genes with at least one context-specific eQTL
sp_eGenes=get_eGenes_multi_tissue(genes_by_tissue = genepos[,-c(2:4)], 
                               snps_by_tissue = snpspos[,-c(2:3)],
                               gene_map = genepos[,1:4], 
                               snp_map = snpspos[,1:3], 
                               nearby = TRUE, dist = 1e06,
                               m_eqtl_out_dir=out_dir, 
                               tissue_names=sort(colnames(genepos)[-c(1:4)]),
                               level1 = level1, level2 = level2, level3 = level3)

fwrite(x = eGenes,file = paste0("sp_eGenes.txt"), sep = '\t', row.names = F, col.names = T)


n_tests_per_SNP <- get_n_tests_per_SNP(snp_map = snpspos[,1:3], 
                                       gene_map= genepos[,1:4], 
                                       nearby = TRUE, dist = 1e06)

sp_eSNPs=get_eSNPs_multi_tissue(genes_by_tissue = genepos[,-c(2:4)], 
                       snps_by_tissue = snpspos[,-c(2:3)], 
                       n_tests_per_SNP = n_tests_per_SNP, 
                       m_eqtl_out_dir=out_dir, 
                       tissue_names=sort(colnames(genepos)[-c(1:4)]),
                       level1 = level1, level2 = level2, level3 = level3)

