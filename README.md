# FastGxC
Computationally efficient and statistically powerful software for detecting context-specific eQTL effects in multi-context genomic studies. Preprint available on BioRxiv at https://www.biorxiv.org/content/10.1101/2021.06.17.448889v1 

Scripts are still under construction but please email us (bballiu@ucla.edu and andrew.lu.chun@gmail.com) with comments / questions. 
Extended data with FastGxC results on GTEx and CLUES cohorts can be found at https://zenodo.org/record/5015123#.YNJ1WpNKjOR

# Running FastGxC

FastGxC works in three steps. 

Step 1: For each individual, FastGxC decompose the phenotype of interest (e.g. gene expression) across C contexts (e.g. tissues or cell-types) into one context-shared and C context-specific components. This is done using the script _decompose_expression.R_.

Step 2: FastGxC estimates genetic effects on the context-shared component and each of the context-specific components separately using simple linear models as implemented in thr R package MatrixEQTL. This is done using the script _run_MatrixEQTL_by_context.R_.

Step 3: FastGxC performs multiple testing adjustment across all contexts, genes, and genetic variants tested using the hierarchical FDR procedures implemented in the R package TreeQTL. This is done using the script _run_TreeQTL.R_ This step can be replaced with other methods to adjust for multiple testing, e.g. mashR https://github.com/stephenslab/mashr. This can lead to a considerable increase in power! 
