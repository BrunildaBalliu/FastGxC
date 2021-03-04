#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% April 2nd 2020
#%%%%%%%%%%%%%%% Script to decompose expression of GTEx samples and run PCA
#%%%%%%%%%%%%%%% For Hoffman
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# qrsh -l h_data=16G,h_rt=12:00:00,highp
# module load R/3.6.1

# Change 
# normalized_expression
# to 
# normalized_and_residualized
# to merge the files residualized for covariates

#%%%%%%%%%%%%%%% R libraries 
library(data.table)
library(reshape2)
library(pcaMethods)
library(magrittr)

#%%%%%%%%%%%%%%% Function to decompose expression into homogeneous and heterogeneous tissue expression
decompose=function(X,design){
  X = as.matrix(X)
  rep.measures = factor(design)
  if (any(summary(as.factor(rep.measures)) == 1)) 
    stop("A multilevel analysis can not be performed when at least one some sample is not repeated.")
  indiv.names = rownames(X)
  rownames(X) = as.character(rep.measures)
  X.mean.indiv = matrix(apply(X, 2, tapply, rep.measures, 
                              mean, na.rm = TRUE), nrow = length(unique(rep.measures)), 
                        ncol = dim(X)[2], dimnames = list(levels(as.factor(rep.measures)), 
                                                          colnames(X)))
  Xb = X.mean.indiv[as.character(rep.measures), ]
  Xw = X - Xb
  dimnames(Xw) = list(indiv.names, colnames(X))
  
  return(list(Xw=Xw,Xb=X.mean.indiv))
}

#%%%%%%%%%%%%%%% Location of data files.
data.dir = '/u/home/b/bballiu/FastGxE/data/GTEx_v8/MatrixEQTL_input/'
cov.dir = '/u/home/b/bballiu/FastGxE/data/GTEx_v8/expression_covariates/'
res.dir = '/u/home/b/bballiu/FastGxE/results/pca/'

#%%%%%%%%%%%%%%% Read expression matrix, genes in columns, samples in rows.
exp_all=data.frame(fread(input = paste0(data.dir,"all_tissues.v8.EUR.normalized_and_residualized_expression_merged.txt"), header = T, sep='\t'),  check.names = F, stringsAsFactors = F, row.names = 1) 

#%%%%%%%%%%%%%%% Print number of genes and samples
sprintf("There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(exp_all), ncol(exp_all),max(colSums(is.na(exp_all))),max(rowSums(is.na(exp_all))))

#%%%%%%%%%%%%%%% Sample and tissue names
sample_tissue_names=matrix(unlist(strsplit(rownames(exp_all), split = " - ")), ncol = 2,byrow = T)
sample_names=sample_tissue_names[,1]
tissue_names=sample_tissue_names[,2]
tissues=unique(tissue_names)

#%%%%%%%%%%%%%%% Decompose expression into homogeneous and heterogeneous tissue expression
print("Decomposing data")
dec_exp_all=decompose(X = exp_all, design = sample_names)
bexp_all=dec_exp_all$Xb
wexp_all=dec_exp_all$Xw
bexp_all[is.nan(bexp_all)]=NA
wexp_all[is.nan(wexp_all)]=NA

# Make sure the output from these two is the same
# Change the sample name
# exp_all[grepl(x = rownames(exp_all), pattern = sample_names[1]),colnames(exp_all)[1]]
# wexp_all[grepl(x = rownames(wexp_all), pattern = sample_names[1]),colnames(exp_all)[1]] +
# bexp_all[grepl(x = rownames(bexp_all), pattern = sample_names[1]),colnames(exp_all)[1]]

sprintf("Between individual matrix: There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(bexp_all), ncol(bexp_all),max(colSums(is.na(bexp_all))),max(rowSums(is.na(bexp_all))))

sprintf("Within individual matrix: There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(wexp_all), ncol(wexp_all),max(colSums(is.na(wexp_all))),max(rowSums(is.na(wexp_all))))

#%%%%%%%%%%%%%%% Save decomposed expression files 
print("Finished decomposition, saving files")

print("Saving between-individuals variation matrix")
fwrite(x = data.table(t(bexp_all),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]},  
       file = paste0(data.dir,"AverageTissue.v8.EUR.normalized_and_residualized_expression_homogeneous.txt"), quote = F, row.names = F, 
       col.names = T, append = F, sep = '\t')

print("Saving within-individuals variation matrix for tissue: ")
for(i in 1:length(tissues)){
  print(tissues[i])
  wexp_t = wexp_all[grep(pattern = tissues[i], rownames(wexp_all)),]
  rownames(wexp_t)=gsub(pattern = paste0(" - ",tissues[i]), replacement = "", x = rownames(wexp_t))
  fwrite(x = data.table(t(wexp_t),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]}, 
         file = paste0(data.dir,tissues[i],".v8.EUR.normalized_and_residualized_expression_heterogeneous.txt"),quote = F, row.names = F, 
         col.names = T, append = F, sep = '\t')
}

#%%%%%%%%%%%%%%% PCA of decomposed data
if(0){
print("Running PCA on full and decomposed data")

#%%%%% Full matrix 
# Remove all samples with NA values across all genes
exp_all_noNA = exp_all[rowSums(is.na(exp_all))<ncol(exp_all)/2,]
# Remove genes with NA values across all samples
exp_all_noNA = exp_all_noNA[,colSums(is.na(exp_all_noNA))<nrow(exp_all_noNA)]
# Keep 10K genes with least amount of NAs
exp_all_noNA = exp_all_noNA[,names(sort(colSums(is.na(exp_all_noNA))))[1:10000]]
# Run PCA and extract PCs
ppca_res=pca(exp_all_noNA, method="ppca", center = T,scale. = T, nPcs=40)
PCs=scores(ppca_res)
fwrite(x = data.table(PCs,keep.rownames = T) %>% {setnames(., old = "rn", new = "sampleID")[]},  
       file = paste0(res.dir,"PCs.v8.EUR.normalized_and_residualized_expression_full.txt"), quote = F, row.names = F, 
       col.names = T, append = F, sep = '\t')
print("Finished extracting and saving PCs from full expression matrix")


#%%%%% Between matrix 
# Remove all samples with NA values across all genes
bexp_all_noNA = bexp_all[rowSums(is.na(bexp_all))<ncol(bexp_all)/2,]
# Remove genes with NA values across all samples
bexp_all_noNA = bexp_all_noNA[,colSums(is.na(bexp_all_noNA))<nrow(bexp_all_noNA)]
# Keep 10K genes with least amount of NAs
bexp_all_noNA = bexp_all_noNA[,names(sort(colSums(is.na(bexp_all_noNA))))[1:10000]]
# Run PCA and extract PCs
b_ppca_res=pca(bexp_all_noNA, method="ppca", center = T,scale. = T, nPcs=40)
bPCs=scores(b_ppca_res)
fwrite(x = data.table(bPCs,keep.rownames = T) %>% {setnames(., old = "rn", new = "sampleID")[]},  
       file = paste0(res.dir,"PCs.v8.EUR.normalized_and_residualized_expression_homogeneous.txt"), quote = F, row.names = F, 
       col.names = T, append = F, sep = '\t')
print("Finished extracting and saving PCs from between-individuals variation matrix")

#%%%%% Within matrix 
# Remove all samples with NA values across all genes
wexp_all_noNA = wexp_all[rowSums(is.na(wexp_all))<ncol(wexp_all)/2,] 
# Remove all genes with NA values across all samples
wexp_all_noNA = wexp_all_noNA[,colSums(is.na(wexp_all_noNA))<nrow(wexp_all_noNA)]
# Keep 10K genes with least amount of NAs
wexp_all_noNA = wexp_all_noNA[,names(sort(colSums(is.na(wexp_all_noNA))))[1:10000]]
# Run PCA and extract PCs
w_ppca_res=pca(wexp_all_noNA, method="ppca", center = T,scale. = T, nPcs=40)
wPCs=scores(w_ppca_res)
fwrite(x = data.table(wPCs,keep.rownames = T) %>% {setnames(., old = "rn", new = "sampleID")[]},  
       file = paste0(res.dir,"PCs.v8.EUR.normalized_and_residualized_expression_heterogeneous.txt"), quote = F, row.names = F, 
       col.names = T,append = F, sep = '\t')
print("Finished extracting and saving PCs from within-individuals variation matrix")
}
