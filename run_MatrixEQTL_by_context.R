#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Andrew Lu & Brunilda Balliu 
#%%%%%%%%%%%%%%% April 2nd 2020
#%%%%%%%%%%%%%%% Script to run Matrix EQTL by content
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Arguments
args=commandArgs(TRUE)
exp_scale=as.numeric(args[1])
i=as.numeric(args[2])

## Libraries and functions
library(data.table)
library(MatrixEQTL)
library(TreeQTL)

setDTthreads(1)
print(paste0("data.table getDTthreads(): ",getDTthreads()))

## Location of data files.
project.dir = getwd(); # change to desired directory

## Assign TISSUE_NAME and scale
if(exp_scale == 1) exp_suffix = ".v8.EUR.normalized_expression"
if(exp_scale == 2) exp_suffix = ".v8.EUR.normalized_and_residualized_expression"
if(exp_scale == 3) exp_suffix = ".v8.EUR.normalized_expression_heterogeneous"
if(exp_scale == 4) exp_suffix = ".v8.EUR.normalized_and_residualized_expression_heterogeneous"
if(exp_scale == 5) exp_suffix = ".v8.EUR.normalized_expression_homogeneous"
if(exp_scale == 6) exp_suffix = ".v8.EUR.normalized_and_residualized_expression_homogeneous"

## Get tissue name
all_files=list.files(path = paste0(project.dir, "/data/GTEx_v8/MatrixEQTL_input/"), pattern = paste0(exp_suffix,".txt"), all.files = FALSE,full.names = F)

if((exp_scale %in% 1:4) & length(all_files)!=49) stop(sprintf("Expecting gene expression files for 49 tissues but got %i.", length(all_files)))

TISSUE_NAME = gsub(x = list.files(path = paste0(project.dir, "/data/GTEx_v8/MatrixEQTL_input/"),
                                  pattern = paste0(exp_suffix,".txt"),
                                  all.files = FALSE,full.names = F),
                   pattern = paste(exp_suffix,".txt",sep = "|"),
                   replacement = "")[i]

sprintf("Running analysis for %s tissue on %s data", TISSUE_NAME, gsub(pattern = "_", replacement = " ", gsub(pattern = ".v8.EUR.",replacement = "",x = exp_suffix)))

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste0(project.dir, "/data/GTEx_v8/MatrixEQTL_input/GTEx_v8_WGS_838Indiv_Analysis_Freeze_EUR_SNPs_1prcMAF.txt");
snps_location_file_name = paste0(project.dir, "/data/GTEx_v8/MatrixEQTL_input/GTEx_v8_WGS_838Indiv_Analysis_Freeze_EUR_SNPs_1prcMAF_snpsloc.txt");

# Gene expression file name
expression_file_name = paste0(project.dir, "/data/GTEx_v8/MatrixEQTL_input/",TISSUE_NAME,exp_suffix,".txt");
gene_location_file_name = paste0(project.dir, "/data/GTEx_v8/MatrixEQTL_input/GTEx_v8_geneloc.txt");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character() #paste(project.dir, "/data/GTEx_v8/expression_covariates/",TISSUE_NAME,"v8.EUR.covariates.txt", sep="");

# Output file name
output_file_name_cis = paste0(project.dir, "/results/eQTL_mapping/MatrixEQTL/",TISSUE_NAME,exp_suffix,".all_pairs.txt");
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1; #if((exp_scale == 1)|(exp_scale == 3)|(exp_scale == 5)) pvOutputThreshold_cis = 0.05 else
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load gene expression data
genes_by_tissue=data.frame(fread(input = paste0(project.dir, "/data/GTEx_v8/misc/GTEx_v8_Gene_by_Tissue_Expression.txt")), row.names = 1)
expression_mat=as.matrix(data.frame(fread(input = expression_file_name, header = T),row.names = 1, check.names = F))
expression_mat=expression_mat[!apply(is.na(expression_mat), 1, all),] # Filter genes with NA across all samples
expression_mat=expression_mat[,!apply(is.na(expression_mat), 2, all)] # Filter samples with NA across all genes
expression_mat=expression_mat[intersect(rownames(expression_mat),names(which(rowSums(genes_by_tissue)!=1))), ] # Filter genes expressed only in single tissue
gene = SlicedData$new();
gene$CreateFromMatrix(expression_mat)

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);
#rownames(snps) = gsub(pattern = "_A|_T|_C|_G",replacement = "", x = snps$GetAllRowNames())

# Match SNP and expression indivuduals
snps$ColumnSubsample(which(colnames(snps) %in% colnames(gene)))

## Run the analysis
snpspos = fread(input = snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# Keep SNP with MAF>5% in each tissue 
if(0){
  passSNPs=data.frame(fread(input = paste0(project.dir,"/data/GTEx_v8/misc/GTEx_v8_SNPs_MAFgeq5_each_tissue.txt"), header = T))
  snps$RowReorder(rownames(snps) %in% passSNPs$SNP);
  snpspos=snpspos[snpspos$snp %in% passSNPs$SNP,]

}

# Keep SNP with MAF>5% in the tissue of interest
if(1){
  MAF=data.frame(fread(input = paste0(project.dir,"/data/GTEx_v8/misc/GTEx_v8_SNPs_by_Tissue_MAF.txt"), header = T,check.names = F),check.names = F)
  passSNPs=MAF[(MAF[,TISSUE_NAME]>=0.05),"SNP"]
  snps$RowReorder(rownames(snps) %in% passSNPs);
  snpspos=snpspos[snpspos$snp %in% passSNPs,]
}

genepos = fread(input = gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos=genepos[genepos$geneid %in% rownames(expression_mat),] # keep positions only for tested genes

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
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
cat('Analyse gedaan in: ', me$time.in.sec, ' seconden', '\n')
