#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% April 21st 2020
#%%%%%%%%%%%%%%% Script to run TreeQTL by content
#%%%%%%%%%%%%%%% For Hoffman
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# qrsh -l h_data=80G,h_rt=12:00:00,highp
# module load R/3.6.1

## Arguments and parameters
args=commandArgs(TRUE) 
exp_scale=as.numeric(args[1])

# FDR thresholds
level1=0.01
level2=0.01
level3=0.01

## Libraries 
library(MatrixEQTL)
library(TreeQTL)
library(dplyr)

# use a single thread
print(paste0("data.table getDTthreads(): ",getDTthreads()))
setDTthreads(1)
print(paste0("after setting as 1 thread; data.table getDTthreads(): ",getDTthreads()))

# Display all warnings as they occur
options(warn=1)

## Functions
get_eGenes_multi_tissue_mod = function (m_eqtl_out_dir, treeQTL_dir, tissue_names, level1 = 0.05, level2 = 0.05, level3 = 0.05, exp_suffix) {
  pattern=paste0(exp_suffix,".all_pairs.txt")
  
  print(paste("Step 0.1: Computing summary statistics for each tissue"))
  m_eqtl_outfiles <- list.files(m_eqtl_out_dir, pattern = pattern, full.names = TRUE)
  if(length(m_eqtl_outfiles)!=49) stop(sprintf("Expecting 49 MatrixEQTL files but got %i.", length(m_eqtl_outfiles)))
  
  n_SNPs_per_gene_outfiles <- list.files(treeQTL_dir, pattern = "n_SNPs_per_gene", full.names = TRUE)
  n_SNPs_per_gene_outfiles=n_SNPs_per_gene_outfiles[!grepl(pattern = "AverageTissue", x = n_SNPs_per_gene_outfiles)]
  if(length(n_SNPs_per_gene_outfiles)!=49) stop(sprintf("Expecting 49 files with nr of SNPs per gene but got %i.", length(n_SNPs_per_gene_outfiles)))

  n_tissue <- length(tissue_names)
  for (i in 1:n_tissue) {
    cur_tissue_name <- tissue_names[i]
    
    print(paste("Computing summary statistics for tissue ", cur_tissue_name, sep = ""))
    n_SNPs_per_gene_this_tissue <- data.frame(fread(input = n_SNPs_per_gene_outfiles[i], header = F), stringsAsFactors = F,check.names = F)
    colnames(n_SNPs_per_gene_this_tissue)=c("family","n_tests")
    n_SNPs_per_gene_this_tissue <- n_SNPs_per_gene_this_tissue[n_SNPs_per_gene_this_tissue$n_tests > 0, ]
    
    gene_simes_cur_tissue <- get_eGenes(n_tests_per_gene = n_SNPs_per_gene_this_tissue, m_eqtl_out = m_eqtl_outfiles[i], method = "BH", level1 = 1, level2 = 1, silent = TRUE)
    gene_simes_cur_tissue <- merge(gene_simes_cur_tissue, n_SNPs_per_gene_this_tissue, by = "family", all = TRUE)
    gene_simes_cur_tissue$fam_p[which(is.na(gene_simes_cur_tissue$fam_p))] <- 1
    
    if (i == 1) {
      eGene_pvals <- gene_simes_cur_tissue[, c("family", "fam_p")]
      n_SNPs_per_gene_xT <- n_SNPs_per_gene_this_tissue
    } else {
      eGene_pvals <- merge(eGene_pvals, gene_simes_cur_tissue[, c("family", "fam_p")], by = "family", all = TRUE)
      n_SNPs_per_gene_xT <- merge(n_SNPs_per_gene_xT, n_SNPs_per_gene_this_tissue, by = "family", all = TRUE)
    }
    names(eGene_pvals)[i + 1] <- cur_tissue_name
    names(n_SNPs_per_gene_xT)[i + 1] <- cur_tissue_name
  }
  names(eGene_pvals)[1] <- "gene"
  remove(cur_tissue_name, n_SNPs_per_gene_this_tissue, gene_simes_cur_tissue)
  
  print("Step 0.2: Computing summary statistics across tissues")
  col_ind_pvals <- 2:(n_tissue + 1)
  eGene_pvals$simes_p <- apply(eGene_pvals[, col_ind_pvals], 1, TreeQTL:::get_simes_p)
  
  print("Step 1: Selecting eGenes across tissues")
  eGene_xT_qvals <- qvalue(eGene_pvals$simes_p, lambda = 0)$qvalue
  R_G <- sum(eGene_xT_qvals <= level1)
  print(paste("Number of cross-tissue eGenes = ", R_G))
  
  print("Step 2: Selecting tissues in which eGenes are active")
  q2_adj <- R_G * level2/nrow(eGene_pvals)
  ind_sel_simes <- which(eGene_xT_qvals <= level1)
  sel_eGenes_simes <- eGene_pvals[ind_sel_simes, ]
  rej_simes <- t(1 * apply(sel_eGenes_simes[, c(col_ind_pvals)], 1, TreeQTL:::qsel_by_fam, q2_adj))
  
  print("Step 3: Selecting SNPs associated to each gene in each tissue")
  sel_eGenes_simes$n_sel_tissues <- rowSums(rej_simes)
  sel_eGenes_simes$n_tested_tissues <- rowSums(!is.na(sel_eGenes_simes[, col_ind_pvals]))
  
  for (i in 1:n_tissue) {
    cur_tissue_name <- tissue_names[i]
    print(paste("Selecting SNPs for tissue", cur_tissue_name))
    sel_gene_names_this_tissue <- sel_eGenes_simes$gene[which(rej_simes[, i] == 1)]
    sel_gene_info <- n_SNPs_per_gene_xT[which(n_SNPs_per_gene_xT$family %in% sel_gene_names_this_tissue), c(1, i + 1)]
    names(sel_gene_info)[2] <- "n_tests"
    sel_gene_info <- merge(sel_gene_info, sel_eGenes_simes[, c("gene", "n_sel_tissues", "n_tested_tissues")], 
                           by.x = "family", by.y = "gene", all.x = TRUE, all.y = FALSE)
    n_sel_per_gene <- TreeQTL:::get_nsel_SNPs_per_gene_tissue_pair(sel_gene_info, cur_tissue_name, m_eqtl_outfiles[i], R_G, nrow(eGene_pvals), 
                                                                   level3 = level3)
    
    print(paste("Total number of associations for tissue", cur_tissue_name, "=", sum(n_sel_per_gene$n_sel_snp)))
    out_file_name <- paste0(treeQTL_dir,"/eAssoc_by_gene.", cur_tissue_name,exp_suffix,".txt")
    print(paste("Writing output file", out_file_name))
    get_eAssociations(data.frame(family = n_sel_per_gene$family, pval = NA, n_sel = n_sel_per_gene$n_sel_snp), NULL, 
                      m_eqtl_outfiles[i], out_file_name, by_snp = FALSE, silent = TRUE)
  }
  eGene_xT_sel <- data.frame(gene = sel_eGenes_simes$gene)
  eGene_xT_sel <- cbind(eGene_xT_sel, rej_simes)
  names(eGene_xT_sel)[2:(n_tissue + 1)] <- tissue_names
  eGene_xT_sel
}

get_pvals_and_fam_p_mod = function(genes_by_tissue, snps_by_tissue, m_eqtl_out_dir, tissue_names, exp_suffix) {
  pattern=paste0(exp_suffix,".all_pairs.txt")
  m_eqtl_outfiles <- list.files(path = m_eqtl_out_dir, pattern = pattern, full.names = TRUE)
  n_tissue <- length(m_eqtl_outfiles)
  pvals_all_tissues <- data.table()
  
  for (i in 1:n_tissue) {
    print(paste("Reading output for tissue ", tissue_names[i], sep = ""))
    cur_data <- data.frame(fread(input = m_eqtl_outfiles[i], header = TRUE,  stringsAsFactors = FALSE))
    cur_data_table <- data.frame(cur_data %>% mutate(pair_names = paste(SNP,gene,sep = "*"))  %>% 
                                   select(pair_names,SNP,gene,p.value),check.names = F, stringsAsFactors = F)
    
    
    names(cur_data_table)[4] <- paste("p.value", i, sep = "_")
    cur_data_table <- data.table(cur_data_table)
    setkey(cur_data_table, pair_names)
    if (sum(duplicated(cur_data_table)) > 0) {
      cur_data_table <- unique(cur_data_table)
      print("Warning: Duplicate key in current output")
    }
    if (i == 1) {
      pvals_all_tissues <- cur_data_table
    } else {
      pvals_all_tissues <- merge(pvals_all_tissues, cur_data_table, by = c("pair_names", "SNP", "gene"), all = TRUE)
    }
  }
  
  print("Calculating p-values for each SNP-gene pair")
  genes_by_tissue <- data.table(genes_by_tissue)
  setkey(genes_by_tissue, gene)
  snps_by_tissue <- data.table(snps_by_tissue)
  setkey(snps_by_tissue, snp)
  pvals_all_tissues$n_tests_pair <- rowSums(genes_by_tissue[J(pvals_all_tissues$gene), 2:ncol(genes_by_tissue), with = FALSE] + 
                                              snps_by_tissue[J(pvals_all_tissues$SNP), 2:ncol(snps_by_tissue), with = FALSE] == 2)
  col_ind_ntests <- which(names(pvals_all_tissues) == "n_tests_pair")
  pvals_all_tissues$fam_p <- apply(pvals_all_tissues[, c(4:(n_tissue + 3), col_ind_ntests), with = FALSE], 1, TreeQTL:::get_simes_p_given_n_tests)
  pvals_all_tissues
}

get_eSNPs_multi_tissue_mod = function(genes_by_tissue, snps_by_tissue, n_tests_per_SNP, m_eqtl_out_dir, tissue_names, level1 = 0.05, level2 = 0.05, level3 = 0.05, exp_suffix) {
  names(snps_by_tissue)[1] <- "snp"
  names(genes_by_tissue)[1] <- "gene"
  
  pvals_all_tissues <- get_pvals_and_fam_p_mod(genes_by_tissue, snps_by_tissue, m_eqtl_out_dir, tissue_names,exp_suffix)
  n_tissue <- length(tissue_names)
  print("Applying error control procedure")
  xT_meqtl_out <- data.frame(SNP = as.character(pvals_all_tissues$SNP), 
                             gene = as.character(pvals_all_tissues$gene), beta = NA, 
                             `t-stat` = NA, `p-value` = pvals_all_tissues$fam_p, FDR = NA)
  xT_meqtl_out <- xT_meqtl_out[order(xT_meqtl_out$p.value), ]
  names(xT_meqtl_out) <- c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
  m_eqtl_out_filename <- tempfile(tmpdir = getwd())
  write.table(xT_meqtl_out, m_eqtl_out_filename, quote = FALSE, row.names = FALSE)
  remove(xT_meqtl_out)
  
  eSNPs <- get_eSNPs(n_tests_per_SNP, m_eqtl_out_filename, method = "BH", level1 = level1, level2 = level2)
  print(paste("Number of eSNPs = ", nrow(eSNPs)))
  
  eAssoc_filename <- tempfile(tmpdir = getwd())
  get_eAssociations(eSNPs, n_tests_per_SNP, m_eqtl_out_filename, eAssoc_filename, by_snp = TRUE)
  level2_sel <- read.table(eAssoc_filename, header = TRUE, stringsAsFactors = FALSE)
  unlink(m_eqtl_out_filename)
  unlink(eAssoc_filename)
  R <- nrow(level2_sel)
  print(paste("Number of selected SNPxgene pairs =", R))
  level2_sel$pair_names <- paste(level2_sel$SNP, "*", level2_sel$gene, sep = "")
  sel_families <- pvals_all_tissues[which(pvals_all_tissues$pair_names %in% level2_sel$pair_names), ]
  n_sel_per_SNP <- table(sel_families$SNP)
  n_sel_per_SNP <- data.frame(SNP = names(n_sel_per_SNP), n_sel = as.numeric(n_sel_per_SNP))
  sel_families <- merge(sel_families, n_sel_per_SNP, by = "SNP", all.x = TRUE)
  names(n_tests_per_SNP) <- c("SNP", "n_genes_per_SNP")
  sel_families <- merge(sel_families, n_tests_per_SNP, by = "SNP", all.x = TRUE)
  sel_families$q_adj <- sel_families$n_sel * nrow(eSNPs) * level3/nrow(n_tests_per_SNP)/sel_families$n_genes_per_SNP
  sel_families$n_sel <- NULL
  sel_families$n_genes_per_SNP <- NULL
  max_recorded <- max(pvals_all_tissues[, 4:(n_tissue + 3), with = FALSE], na.rm = TRUE)
  if (min(sel_families$q_adj) >= max_recorded) {
    warning("Matrix eQTL output threshold may be too small for given levels")
  }
  col_ind_pvals <- which(sapply(names(sel_families), grep, fixed = TRUE, pattern = "p.value") == 1)
  col_ind_ntests <- which(names(sel_families) == "n_tests_pair")
  col_ind_qadj <- which(names(sel_families) == "q_adj")
  rej_by_fam_q <- apply(sel_families[, c(col_ind_pvals, col_ind_ntests, col_ind_qadj), with = FALSE], 1, TreeQTL:::bh_by_fam_q_adj)
  sel_pairs <- data.frame(matrix(t(rej_by_fam_q), nrow = ncol(rej_by_fam_q)))
  return_value <- data.frame(family = sel_families$pair_names, fam_p = sel_families$fam_p, as.matrix(sel_families[, col_ind_pvals, with = FALSE]) * sel_pairs)
  names(return_value) <- c("family", "fam_p", tissue_names)
  return_value$family <- as.character(return_value$family)
  return_value[is.na(return_value)] <- 0
  return_value[, 2] <- as.character(sapply(as.character(return_value$family), TreeQTL:::get_gene_name))
  return_value[, 1] <- as.character(sapply(as.character(return_value$family), TreeQTL:::get_snp_name))
  names(return_value)[1:2] <- c("SNP", "gene")
  return_value
}

## Location of data files.
project.dir = '/u/project/zaitlenlab/bballiu/FastGxE';

## Assign TISSUE_NAME and scale
if(exp_scale == 1) exp_suffix = ".v8.EUR.normalized_expression"
if(exp_scale == 2) exp_suffix = ".v8.EUR.normalized_and_residualized_expression"
if(exp_scale == 3) exp_suffix = ".v8.EUR.normalized_expression_heterogeneous"
if(exp_scale == 4) exp_suffix = ".v8.EUR.normalized_and_residualized_expression_heterogeneous"
if(exp_scale == 5) exp_suffix = ".v8.EUR.normalized_expression_homogeneous"
if(exp_scale == 6) exp_suffix = ".v8.EUR.normalized_and_residualized_expression_homogeneous"

#### Perform hierarchical selection for eGenes across multiple tissues by first selecting eGenes which are subject to genetic regulation in any tissue, then the tissues in which they are regulated, and finally the specific SNPs associated to each selected gene x tissue pair

#m_eqtl_out_dir: path to output files with eQTL output per tissue. 
# These are assumed to follow the output format generated by Matrix eQTL.
m_eqtl_out_dir=paste0(project.dir, "/results/eQTL_mapping/MatrixEQTL_FDRthreshold");
# m_eqtl_out_dir=paste0(project.dir, "/results/eQTL_mapping/MatrixEQTL");
treeQTL_dir=paste0(project.dir, "/results/eQTL_mapping/TreeQTL");

#tissue_names: vector of names for each tissue in alphanumeric order
pattern=paste0(exp_suffix,".all_pairs.txt")
tissue_names=sort(gsub(pattern = paste(exp_suffix,".all_pairs.txt",sep = "|"), replacement = "", x = list.files(m_eqtl_out_dir, pattern = pattern, full.names = FALSE)))

# Get eGenes. 
if(exp_scale %in% 1:4){
  
  if(length(tissue_names)!=49) stop(sprintf("Expecting 49 tissues but got %i.", length(tissue_names)))
  
  eGenes = get_eGenes_multi_tissue_mod(m_eqtl_out_dir = m_eqtl_out_dir, 
                                       treeQTL_dir=treeQTL_dir, 
                                       tissue_names=tissue_names, 
                                       level1 = level1, level2 = level1, level3 = level1, 
                                       exp_suffix=exp_suffix)
  fwrite(x = eGenes,file = paste0(treeQTL_dir,"/eGenes",exp_suffix,".txt"), sep = '\t', row.names = F, col.names = T)
  
  
}

if(exp_scale %in% 5:6){
  
  m_eqtl_outfile <- list.files(m_eqtl_out_dir, pattern = pattern, full.names = TRUE)
  n_SNPs_per_gene_outfile <- list.files(treeQTL_dir, pattern = "n_SNPs_per_gene_AverageTissue.txt", full.names = TRUE)
  
  n_tests_per_gene <- data.frame(fread(input = n_SNPs_per_gene_outfile, header = F), stringsAsFactors = F,check.names = F)
  colnames(n_tests_per_gene)=c("family","n_tests")
  n_tests_per_gene <- n_tests_per_gene[n_tests_per_gene$n_tests > 0, ]
  
  eGenes=get_eGenes(n_tests_per_gene = n_tests_per_gene, m_eqtl_out=m_eqtl_outfile, method = "BY", level1 = level1, level2 = level2) 
  fwrite(x = eGenes,file = paste0(treeQTL_dir,"/eGenes",exp_suffix,".txt"), sep = '\t', row.names = F, col.names = T)
  
  get_eAssociations(eDiscoveries=eGenes, n_tests =n_tests_per_gene, m_eqtl_out = m_eqtl_outfile, 
                    out_file = paste0(treeQTL_dir,"/eAssoc_by_gene.AverageTissue",exp_suffix,".txt"), by_snp = FALSE)
}

# Get eSNPs
if(0){
if(exp_scale %in% 1:4){
  
  snps_by_tissue_file=paste0(project.dir, "/data/GTEx_v8/misc/GTEx_v8_SNPs_by_Tissue.txt")
  if(file.exists(snps_by_tissue_file)){
    print("Reading snps-by-tissue file")
    snps_by_tissue=fread(input = snps_by_tissue_file, header = T, sep = ',')  
  } else {
    print("Making snps-by-tissue file")
    MAF=data.frame(fread(input = paste0(project.dir,"/data/GTEx_v8/misc/GTEx_v8_SNPs_by_Tissue_MAF.txt"), header = T))
    snps_by_tissue = MAF
    for(i in 2:ncol(snps_by_tissue)){
      snps_by_tissue[,i] = (rowSums(MAF[,-1]>=0.05)>=2)&(MAF[,i]>=0.05)
    }
    fwrite(x = snps_by_tissue, file = paste0(project.dir, "/data/GTEx_v8/misc/GTEx_v8_SNPs_by_Tissue.txt"), col.names = T, row.names = F)
  }
  
  genes_by_tissue = data.frame(fread(input = "/u/project/zaitlenlab/bballiu/FastGxE/data/GTEx_v8/misc/GTEx_v8_Gene_by_Tissue_Expression.txt", header = T), 
                               check.names = F,stringsAsFactors = F)
  genes_by_tissue[rowSums(genes_by_tissue[,-1])<2,-1]=0
  
  n_tests_per_SNP=data.frame(fread(input = paste0(project.dir, "/results/eQTL_mapping/TreeQTL/n_genes_per_SNPs_Heterogeneous.txt"),
                                   header = F),check.names = F, stringsAsFactors = F)
  colnames(n_tests_per_SNP)=c("family", "n_tests")

  eSNPs = get_eSNPs_multi_tissue_mod(genes_by_tissue = genes_by_tissue, snps_by_tissue = snps_by_tissue, 
                                     n_tests_per_SNP = n_tests_per_SNP, m_eqtl_out_dir = m_eqtl_out_dir, 
                                     tissue_names=tissue_names, level1 = level1, level2 = level2, level3 = level3, 
                                     exp_suffix=exp_suffix)   
  
  fwrite(x = eSNPs,file = paste0(treeQTL_dir,"/eSNPs",exp_suffix,".txt"), sep = '\t', row.names = F, col.names = T)
  
}

if(exp_scale %in% 5:6){
  
  n_tests_per_SNP=data.frame(fread(input = paste0(project.dir, "/results/eQTL_mapping/TreeQTL/n_genes_per_SNPs_Homogeneous.txt"),
                                   header = F),check.names = F, stringsAsFactors = F)
  colnames(n_tests_per_SNP)=c("family", "n_tests")
  
  m_eqtl_out <- list.files(m_eqtl_out_dir, pattern = pattern, full.names = TRUE)
  
  eSNPs = get_eSNPs(n_tests_per_SNP = n_tests_per_SNP, m_eqtl_out = m_eqtl_out, level1 = level1, level2 = level2)
  fwrite(x = eSNPs,file = paste0(treeQTL_dir,"/eSNPs",exp_suffix,".txt"), sep = '\t', row.names = F, col.names = T)
  
}
}


print("run_TreeQTL.R finished!")