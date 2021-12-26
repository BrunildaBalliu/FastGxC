library(TreeQTL)



snpspos = read.table(file = snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(file = gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

eGenes=get_eGenes_multi_tissue(genes_by_tissue = genepos[,-c(2:4)], 
                               snps_by_tissue = snpspos[,-c(2:3)],
                               gene_map = genepos[,1:4], 
                               snp_map = snpspos[,1:3], 
                               nearby = TRUE, dist = 1e06,
                               m_eqtl_out_dir=out_dir, 
                               tissue_names=sort(colnames(genepos)[-c(1:4)]),
                               level1 = 0.05, level2 = 0.05, level3 = 0.05)


