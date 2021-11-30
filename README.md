# FastGxC
Computationally efficient and statistically powerful software for detecting context-specific eQTL effects in multi-context genomic studies. 

Preprint available on BioRxiv at https://www.biorxiv.org/content/10.1101/2021.06.17.448889v1 

Extended data with FastGxC results on GTEx and CLUES cohorts can be found at https://zenodo.org/record/5015123#.YNJ1WpNKjOR

Scripts are still under construction but please email us (bballiu@ucla.edu and andrew.lu.chun@gmail.com) with comments / questions. 

# Running FastGxC

FastGxC works in three steps. 

*Step 0*: If you want to run a toy example, you can generate simulated data by running 

  Rscript generate_simulated_data.R "/Users/bballiu/Documents/GitHub/FastGxC"

This script will make a _data_ folder in _your_project_directory_ (if one does not already exists) and generate and save the following files 
(1) SNPs.txt: snp genotype data for 50 SNPs and 300 individuals (MatrixEQTL input format), 
(2) snpsloc.txt: location information of these 50 SNPs (MatrixEQTL input format), 
(3) expression.txt: gene expression data for 300 individuals across 10 genes and 50 contexts
(2) geneloc.txt: location information of these 10 genes (MatrixEQTL input format), 

*Step 1*: For each individual, decompose the phenotype of interest (e.g. gene expression) across C contexts (e.g. tissues or cell-types) into one context-shared and C context-specific components by running. 
  
  Rscript decompose_expression.R "/Users/bballiu/Documents/GitHub/FastGxC" "expression.txt"

This script will take as an imput a file with gene expression data for all individuals, genes, and contexts (see _expression.txt_ for right format) and outputs one file with context-shared expression (context_shared_expression.txt) and C files with expression specific to each context (CONTEXT_NAME_specific_expression.txt). 

*Step 2*: FastGxC estimates genetic effects on the context-shared component and each of the context-specific components separately using simple linear models. Here we use the R package MatrixEQTL but any other software that can perform quick linear regression can be used. 

  Rscript run_MatrixEQTL_by_context.R "/Users/bballiu/Documents/GitHub/FastGxC" SNPs.txt snpsloc.txt context_shared_expression.txt geneloc.txt "context_shared_eQTLs.txt"

  for i in $(seq 1 5); do
    Rscript run_MatrixEQTL_by_context.R "/Users/bballiu/Documents/GitHub/FastGxC" SNPs.txt snpsloc.txt context$i\_specific\_expression.txt geneloc.txt context$i\_specific\_eQTLs.txt
  done



*Step 3*: FastGxC performs multiple testing adjustment across all contexts, genes, and genetic variants tested using the hierarchical FDR procedures implemented in the R package TreeQTL. This is done using the script _run_TreeQTL.R_ This step can be replaced with other methods to adjust for multiple testing, e.g. mashR https://github.com/stephenslab/mashr. This can lead to a considerable increase in power! 
