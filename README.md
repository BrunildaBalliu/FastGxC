# FastGxC
Computationally efficient and statistically powerful software for detecting context-specific eQTL effects in multi-context genomic studies. 

Preprint available on [BioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.17.448889v1) 

Extended data with FastGxC results on GTEx and CLUES cohorts can be found [here](https://zenodo.org/record/5015123#.YNJ1WpNKjOR)

Scripts are still under construction but please email me (bballiu at ucla dot edu) with comments / questions. 

# Simulate toy data
Please download the required R packages inside `generate_simulated_data.R` before running the toy example. 

If you want to run a toy example, you can generate simulated data by running 
```
  project_directory=your_project_directory
  
  Rscript generate_simulated_data.R $project_directory
```

This script will make a _data_ folder in your _project_directory_ (if one does not already exists) and generate and save the following files 
(1) SNPs.txt: snp genotype data for 10,000 SNPs and 300 individuals (MatrixEQTL input format), 
(2) snpsloc.txt: location information of these 10,000 SNPs (MatrixEQTL input format), 
(3) expression.txt: gene expression data for 300 individuals across 100 genes and 50 contexts
(2) geneloc.txt: location information of these 10 genes (MatrixEQTL input format), 

# Running FastGxC
Please download the required R packages inside `decompose_expression.R` and `run_MatrixEQTL.R` before running FastGxC. 

FastGxC works in two steps. In the first step, expression is decomposed into shared and context-specific components. In the second step, eQTLs are separately mapped on these components.

*Step 1 - Decomposition:* For each individual, decompose the phenotype of interest (e.g. gene expression) across C contexts (e.g. tissues or cell-types) into one context-shared and C context-specific components by running. 
  
  ```
  project_directory=your_project_directory
  exp_file_name="expression.txt"
  Rscript decompose_expression.R $project_directory $exp_file_name
  ```

This script will take as an imput a file with gene expression data for all individuals, genes, and contexts (see _expression.txt_ for right format) and outputs one file with context-shared expression (context_shared_expression.txt) and C files with expression specific to each context (CONTEXT_NAME_specific_expression.txt). The files will be saved in the  _data_ folder in your _project_directory_. 

*Step 2 - eQTL mapping:* FastGxC estimates genetic effects on the context-shared component and each of the C context-specific components separately using simple linear models. Note: Here we use the R package [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) but any other software that can perform quick linear regression can be used (e.g. [FastQTL](http://fastqtl.sourceforge.net/) or [tensorqtl](https://github.com/broadinstitute/tensorqtl)). 

Map context-specific eQTLs
```  
project_directory=your_project_directory

for i in $(seq 1 50); do
    Rscript run_MatrixEQTL.R $project_directory SNPs.txt snpsloc.txt context$i\_specific\_expression.txt geneloc.txt  context$i\_specific\_eQTLs.txt specific
 done
```
Map context-shared eQTLs
```
  Rscript run_MatrixEQTL.R $project_directory SNPs.txt snpsloc.txt context_shared_expression.txt geneloc.txt  "context_shared_eQTLs.txt" shared
```

This script take as input data needed to run MatrixEQTL and outputs eQTL summary statistics in the MatrixEQTL format. In the end, you should have one file with summary statistics for shared eQTL and C files with summary statistics for each context C. 

# Multiple testing adjustment

Please download the required R packages inside `run_TreeQTL.R` before running TreeQTL. 

To adjust for multiple testing across all contexts, genes, and genetic variants tested you can use the hierarchical FDR procedures implemented in the R package [TreeQTL](http://bioinformatics.org/treeqtl/). 

TreeQTL requires that you run MatrixEQTL to do eQTL mapping (see step 2 above). If you used another eQTL mapping softwares, please make sure the output is in the format required by TreeQTL. You can also replace TreeQTL with other methods, e.g. [mashR](https://github.com/stephenslab/mashr), which can also lead to a considerable increase in power. 

Map specific-eGenes, i.e., genes with at least one context-specific eQT  
```  
project_directory=your_project_directory
Rscript run_TreeQTL.R $project_directory specific
```

Map shared-eGenes, i.e., genes with at least one context-shared eQT  
```  
project_directory=your_project_directory
Rscript run_TreeQTL.R $project_directory shared
```

This script take as input data needed to run TreeQTL and outputs shared and specific eGenes (two files) and eAssociation (C+1 files) summary statistics in the TreeQTL format. 


