#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Script to run Matrix EQTL by content
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% December 1st 2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%% Input and output files
#%%%%%%%%%%%%%%%%%%%%%%%%

args=commandArgs(TRUE)
work_dir = args[1]
data_dir=paste0(work_dir,'/data/') # directory with phenotype and genotype data
out_dir=paste0(work_dir,'/results/') # MatrixEQTL output
if(!dir.exists(out_dir)) dir.create(out_dir)

# Genotype file name
SNP_file_name = paste0(data_dir, args[2]); 
snps_location_file_name = paste0(data_dir, args[3]); 

# Gene expression file name
expression_file_name = paste0(data_dir, args[4]); 
gene_location_file_name = paste0(data_dir, args[5]); 

# Output file name
output_file_name_cis = paste0(out_dir, args[6]); 
output_file_name_tra = tempfile();

# work_dir = getwd()
# data_dir=paste0(work_dir,'/data/') # directory with phenotype and genotype data
# out_dir=paste0(work_dir,'/results/') # MatrixEQTL output
# SNP_file_name = paste0(data_dir, "SNPs.txt");
# snps_location_file_name = paste0(data_dir, "snpsloc.txt");
# expression_file_name = paste0(data_dir, "context_shared_expression.txt");
# gene_location_file_name = paste0(data_dir, "geneloc.txt");
# output_file_name_cis = paste0(out_dir, "context_shared_eQTLs.txt");
# output_file_name_tra = tempfile();

#%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%% Libraries and functions
#%%%%%%%%%%%%%%%%%%%%%%%%
library(data.table)
library(MatrixEQTL)

setDTthreads(1)
print(paste0("data.table getDTthreads(): ",getDTthreads()))

sprintf("Running analysis for %s", args[4])

#%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%% MatrixEQTL parameters
#%%%%%%%%%%%%%%%%%%%%%%%%

# Linear model to use
useModel = modelLINEAR; 

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1; 
pvOutputThreshold_tra = 0;

# Error covariance matrix
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 1e6;

#%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%% Read files
#%%%%%%%%%%%%%%%%%%%%%%%%

# genes_by_tissue=data.frame(fread(input = paste0(work_dir, "/data/GTEx_v8/misc/GTEx_v8_Gene_by_Tissue_Expression.txt")), row.names = 1)

## Raw gene expression data with gene position
expression_mat=as.matrix(data.frame(fread(input = expression_file_name, header = T),row.names = 1, check.names = F))
genepos = read.table(file = gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)[,1:4];

## Genotype data with snp position
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

snpspos = read.table(file = snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)[,1:3];

## Make sure order of individuals is the same in gene expression and genotype matrices  
snps$ColumnSubsample(which(colnames(snps) %in% colnames(expression_mat))) # Match SNP and expression individuals
expression_mat=expression_mat[,colnames(snps)]
gene = SlicedData$new();
gene$CreateFromMatrix(expression_mat)

#%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%% Run the analysis
#%%%%%%%%%%%%%%%%%%%%%%%%

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = SlicedData$new(),
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

## Results:
cat('Analysis finished in: ', me$time.in.sec, ' seconds', '\n')

