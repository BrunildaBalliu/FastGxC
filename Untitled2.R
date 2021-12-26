genes_by_tissue = genepos[,-c(2:4)]
snps_by_tissue = snpspos[,-c(2:3)]
gene_map = genepos[,1:4]
snp_map = snpspos[,1:3] 
nearby = TRUE 
dist = 1e06
m_eqtl_out_dir=out_dir 
tissue_names=sort(colnames(genepos)[-c(1:4)])
level1 = 0.05 
level2 = 0.05 
level3 = 0.05


names(snp_map)[1] <- "snp"
names(gene_map)[1] <- "gene"
names(snps_by_tissue)[1] <- "snp"
names(genes_by_tissue)[1] <- "gene"
print(paste("Step 0.1: Computing summary statistics for each tissue"))
m_eqtl_outfiles <- list.files(m_eqtl_out_dir, full.names = TRUE)
n_tissue <- length(tissue_names)
for (i in 1:n_tissue) {
  cur_tissue_name <- tissue_names[i]
  print(paste("Computing summary statistics for tissue",cur_tissue_name, sep = ""))
  snps_this_tissue <- snps_by_tissue$snp[which(snps_by_tissue[, i + 1] == 1)]
  genes_this_tissue <- genes_by_tissue$gene[which(genes_by_tissue[, i + 1] == 1)]
  n_SNPs_per_gene_this_tissue <- get_n_tests_per_gene(snp_map[which(snp_map$snp %in% snps_this_tissue), ], 
                                                      gene_map[which(gene_map$gene %in% genes_this_tissue), ], 
                                                      nearby = nearby, dist = dist)
  n_SNPs_per_gene_this_tissue <- n_SNPs_per_gene_this_tissue[n_SNPs_per_gene_this_tissue$n_tests > 0, ]
  gene_simes_cur_tissue <- get_eGenes(n_SNPs_per_gene_this_tissue, m_eqtl_outfiles[i], method = "BH", level1 = 1, level2 = 1, silent = TRUE)
  gene_simes_cur_tissue <- merge(gene_simes_cur_tissue, n_SNPs_per_gene_this_tissue, by = "family", all = TRUE)
  gene_simes_cur_tissue$fam_p[which(is.na(gene_simes_cur_tissue$fam_p))] <- 1
  if (i == 1) {
    eGene_pvals <- gene_simes_cur_tissue[, c("family", 
                                             "fam_p")]
    n_SNPs_per_gene_xT <- n_SNPs_per_gene_this_tissue
  } else {
    eGene_pvals <- merge(eGene_pvals, gene_simes_cur_tissue[, 
                                                            c("family", "fam_p")], by = "family", all = TRUE)
    n_SNPs_per_gene_xT <- merge(n_SNPs_per_gene_xT, n_SNPs_per_gene_this_tissue, 
                                by = "family", all = TRUE)
  }
  names(eGene_pvals)[i + 1] <- cur_tissue_name
  names(n_SNPs_per_gene_xT)[i + 1] <- cur_tissue_name
}


names(eGene_pvals)[1] <- "gene"
remove(cur_tissue_name, n_SNPs_per_gene_this_tissue, gene_simes_cur_tissue, 
       snps_this_tissue, genes_this_tissue)
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
m_eqtl_out_files <- list.files(m_eqtl_out_dir)

for (i in 1:n_tissue) {
  cur_tissue_name <- tissue_names[i]
  print(paste("Selecting SNPs for tissue", cur_tissue_name))
  sel_gene_names_this_tissue <- sel_eGenes_simes$gene[which(rej_simes[, i] == 1)]
  sel_gene_info <- n_SNPs_per_gene_xT[which(n_SNPs_per_gene_xT$family %in% sel_gene_names_this_tissue), c(1, i + 1)]
  names(sel_gene_info)[2] <- "n_tests"
  sel_gene_info <- merge(sel_gene_info, sel_eGenes_simes[, 
                                                         c("gene", "n_sel_tissues", "n_tested_tissues")], 
                         by.x = "family", by.y = "gene", all.x = TRUE, all.y = FALSE)
  n_sel_per_gene <- TreeQTL:::get_nsel_SNPs_per_gene_tissue_pair(sel_gene_info, cur_tissue_name, m_eqtl_outfiles[i], R_G, nrow(eGene_pvals), level3 = level3)
  print(paste("Total number of associations for tissue", 
              cur_tissue_name, "=", sum(n_sel_per_gene$n_sel_snp)))
  out_file_name <- paste("eAssoc_3level_by_gene_tissue_", 
                         cur_tissue_name, ".txt", sep = "")
  print(paste("Writing output file", out_file_name))
  get_eAssociations(data.frame(family = n_sel_per_gene$family,  pval = NA, n_sel = n_sel_per_gene$n_sel_snp), NULL, m_eqtl_outfiles[i], out_file_name, by_snp = FALSE, 
                    silent = TRUE)
}

eGene_xT_sel <- data.frame(gene = sel_eGenes_simes$gene)
eGene_xT_sel <- cbind(eGene_xT_sel, rej_simes)
names(eGene_xT_sel)[2:(n_tissue + 1)] <- tissue_names
eGene_xT_sel
}