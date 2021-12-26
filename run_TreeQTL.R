#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Script to run TreeQTL by content
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% December 1st 2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Arguments and parameters
args=commandArgs(TRUE) 
work_dir = args[1];
data_dir=paste0(work_dir,'data/') # directory with phenotype and genotype data
out_dir=paste0(work_dir,'results_',args[2],'/') # MatrixEQTL output
source(paste0(work_dir, "functions.R"))

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

# Distance for local gene-SNP pairs
cisDist = 1e6;


# Gene and SNP positions
snps_location_file_name = paste0(data_dir, "snpsloc.txt"); 
gene_location_file_name = paste0(data_dir, "geneloc.txt"); 

snpspos = read.table(file = snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(file = gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# Use treeQTL to perform hierarchical FDR and get specific_eGenes, i.e. genes with at least one context-specific eQTL, and shared_eGenes, i.e. genes with at least one context-shared eQTL
if(args[2] == "specific"){ 
  
  specific_eGenes=get_eGenes_multi_tissue(genes_by_tissue = genepos[,-c(2:4)], 
                                          snps_by_tissue = snpspos[,-c(2:3)],
                                          gene_map = genepos[,1:4], 
                                          snp_map = snpspos[,1:3], 
                                          nearby = TRUE, 
                                          dist = cisDist,
                                          m_eqtl_out_dir=out_dir, 
                                          tissue_names=sort(colnames(genepos)[-c(1:4)]),
                                          level1 = level1, level2 = level2, level3 = level3)
  write.table(x = specific_eGenes, file = paste0(work_dir,"specific_eGenes.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
  
  
}

if(args[2] == "shared"){ 
  
  n_tests_per_gene = get_n_tests_per_gene(snp_map = snpspos[,1:3], gene_map = genepos[,1:4], 
                                          nearby = TRUE, dist = cisDist)

  
  shared_eGenes = get_eGenes(n_tests_per_gene = n_tests_per_gene, 
                             m_eqtl_out = paste0(out_dir,"context_shared_eQTLs.txt"), 
                             method = "BH",
                             level1 = level1, level2 = level2, 
                             slice_size = 1e+05,
                             silent = FALSE, 
                             gene_pvals = NA)
  
  write.table(x = shared_eGenes, file = paste0(work_dir,"shared_eGenes.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
  
  
  eAssociations = get_eAssociations(eDiscoveries = shared_eGenes, n_tests = n_tests_per_gene, 
                    m_eqtl_out = paste0(out_dir,"context_shared_eQTLs.txt"), 
                    out_file = paste0(work_dir,"eAssoc_2level_by_gene_context_shared.txt"), 
                    by_snp = F, slice_size = 1e+05,
                    silent = FALSE)
  
}
