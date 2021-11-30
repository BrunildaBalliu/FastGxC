#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% April 2nd 2020
#%%%%%%%%%%%%%%% Script to decompose expression into shared and context specific components
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%  Directories for input and output files 
args=commandArgs(TRUE)
work_dir = args[1]
exp_mat_filename = args[2]
data_dir=paste0(work_dir,'/data/') # directory with phenotype and genotype data
if(!dir.exists(data_dir)) dir.create(data_dir)

#%%%%%%%%%%%%%%% R libraries 
library(data.table)
library(reshape2)
library(magrittr)

#%%%%%%%%%%%%%%% Function to decompose expression into context-shared and context-specific components
decompose=function(X,design){
  X = as.matrix(X)
  rep.measures = factor(design)
  if (any(summary(as.factor(rep.measures)) == 1)) 
    stop("A multilevel analysis can not be performed when at least one some sample is not repeated.")
  indiv.names = rownames(X)
  rownames(X) = as.character(rep.measures)
  X.mean.indiv = matrix(apply(X, 2, tapply, rep.measures, mean, na.rm = TRUE), 
                        nrow = length(unique(rep.measures)), 
                        ncol = dim(X)[2], 
                        dimnames = list(levels(as.factor(rep.measures)), colnames(X)))
  Xb = X.mean.indiv[as.character(rep.measures), ]
  Xw = X - Xb
  dimnames(Xw) = list(indiv.names, colnames(X))
  return(list(Xw=Xw,Xb=X.mean.indiv))
}

#%%%%%%%%%%%%%%% Read expression matrix, genes in columns, samples in rows.
exp_mat=read.table(file = paste0(data_dir,exp_mat_filename), sep = '\t')

#%%%%%%%%%%%%%%% Print number of genes and samples
sprintf("There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(exp_mat), ncol(exp_mat),max(colSums(is.na(exp_mat))),max(rowSums(is.na(exp_mat))))

#%%%%%%%%%%%%%%% Sample and context names
design=sapply(1:nrow(exp_mat), function(i) unlist(strsplit(rownames(exp_mat)[i], split = "_"))[1])
context_names=sapply(1:nrow(exp_mat), function(i) unlist(strsplit(rownames(exp_mat)[i], split = "_"))[2])
contexts=unique(context_names)

#%%%%%%%%%%%%%%% Decompose expression into homogeneous and heterogeneous context expression
print("Decomposing data")
dec_exp_all=decompose(X = exp_mat, design = design)
bexp_all=dec_exp_all$Xb
wexp_all=dec_exp_all$Xw
bexp_all[is.nan(bexp_all)]=NA
wexp_all[is.nan(wexp_all)]=NA

sprintf("Between individual matrix: There are %s individuals and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(bexp_all), ncol(bexp_all),max(colSums(is.na(bexp_all))),max(rowSums(is.na(bexp_all))))

sprintf("Within individual matrix: There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(wexp_all), ncol(wexp_all),max(colSums(is.na(wexp_all))),max(rowSums(is.na(wexp_all))))

#%%%%%%%%%%%%%%% Save decomposed expression files 
print("Finished decomposition, saving files")

print("Saving between-individuals variation matrix")
fwrite(x = data.table(t(bexp_all),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]}, 
       file = paste0(data_dir,"context_shared_expression.txt"), quote = F, row.names = F, 
       col.names = T, append = F, sep = '\t')

print("Saving within-individuals variation matrix for context: ")
for(i in 1:length(contexts)){
  print(contexts[i])
  wexp_t = wexp_all[grep(pattern = contexts[i], rownames(wexp_all)),]
  rownames(wexp_t)=gsub(pattern = paste0("_",contexts[i]), replacement = "", x = rownames(wexp_t))
  fwrite(x = data.table(t(wexp_t),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]}, 
         file = paste0(data_dir,contexts[i],"_specific_expression.txt"),quote = F, row.names = F, 
         col.names = T, append = F, sep = '\t')
}

