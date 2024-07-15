pacman::p_load(SpatialExperiment, batchelor, BiocParallel, readr, purrr, dplyr, 
               dittoSeq, viridis, ggplot2, ggrepel, scater, corrplot, cowplot, 
               scuttle, caret, RcppAlgos)

spe<-readRDS('spe_raw.rds')


# Scaling
assay(spe, 'exprs')<-t(scale(t(assay(spe, 'counts')), center= rowMeans(assay(spe, 'counts')), scale=rowSds(assay(spe, 'counts'))))

# Removing cels positive for more than 10 markers
spe<-spe[,colSums(assay(spe, 'exprs')>1)<=10]

# Rescaling
assay(spe, 'exprs')<-t(scale(t(assay(spe, 'counts')), center= rowMeans(assay(spe, 'counts')), scale=rowSds(assay(spe, 'counts'))))

saveRDS(spe, 'spe.rds')

# Co-expression analysis
co_exp_matrix<-matrix(0, 47, 47)
rownames(co_exp_matrix)<- colnames(co_exp_matrix) <- rownames(spe)[1:47]
for(i in 1:46) {
  i_intensities <- assay(spe, "exprs")[i,]
  for (j in (i+1):47){
    j_intensities<- assay(spe, "exprs")[j,]
    co_exp_matrix[i,j] <- co_exp_matrix[j,i] <- sum(j_intensities>1 & i_intensities>1)/min(sum(i_intensities>1), sum(j_intensities>1))
    
  }
}
diag(co_exp_matrix)<-1
dir.create('data')
write.csv(co_exp_matrix, 'data/co_exp_matrix.csv')

# Category assignment & Cell typing

cat_type_map<-list(
  'Epithelial'=c( 'Neuroendocrine', 'Goblet', 'Paneth', 'TA'), 
  'Stromal'=c('Smooth muscle', 'Endothelial', 'Nerve', 'Lymphatic', 'ICC'),
  'Immune'=c('CD8+ T', 'CD4+ T', 'B', 'Plasma', 'M1 Macrophage', 'M2 Macrophage', 'Neutrophil', 'DC', 'NK')
)
cell_categories<-list(
  'Epithelial'=c('Cytokeratin', 'SOX9', 'aDef5'), 
  'Immune'=c('CD45', 'CD38', 'CD161', 'CD16', 'CD163', 'CD68', 'CD11c')
)
category_markers<-list_c(cell_categories)


cell_types<-list(
  'Goblet'= c('MUC2'), 
  'Paneth'= c('aDef5'),
  'Neuroendocrine'=c('CHGA'), 
  'TA'=c('SOX9'), 
  
  'B'=c('CD19', 'CD21'),
  'CD4+ T'=c('CD4', 'CD3'),
  'CD8+ T'=c('CD8', 'CD3'),
  'NK' =c('CD161', 'NKG2D'),
  'Plasma'=c('CD38', 'CD138'),
  'DC'=c('HLADR', 'CD11c'),
  'Neutrophil'=c('CD15', 'CD66'), 
  'M1 Macrophage'=c('CD68', 'HLADR'), 
  'M2 Macrophage'=c('CD163', 'CD206'), 
  
  'Smooth muscle'=c('aSMA'),
  'Endothelial'=c('CD31'),
  'Lymphatic'=c('Podoplanin'),
  'Nerve'=c('Synapto', 'CD56'),
  'ICC'=c('CD117') 
)
type_markers<-unique(list_c(cell_types)) 


assign_cat<-function(exprs){
  
  if(any(exprs[cell_categories[['Epithelial']]]>1)){
    if(any(exprs[cell_categories[['Immune']]]>1)) return('unknown')
    else return('Epithelial')
  }
  if(any(exprs[cell_categories[['Immune']]]>1)) return('Immune')
  if(all(exprs[cell_categories[['Epithelial']]]< 0) & all(exprs[cell_categories[['Immune']]]< 0)) return('Stromal')
  return('unknown')
  
}  

assign_type<-function(exprs, cat){
  cat_markers<-unique(list_c(cell_types[cat_type_map[[cat]]]))
  scores<-sapply(cat_type_map[[cat]], function(t)  score<- sum(exprs[cell_types[[t]]]>1)/length(cell_types[[t]]) - (sum(exprs[setdiff(cat_markers,cell_types[[t]])]>1)/length(setdiff(cat_markers,cell_types[[t]]))))
  # Adjust for the NK cells, otherwise they are to many
  if(cat=='Immune') scores['NK']<- sum(exprs[cell_types[['NK']]]>2)/2 - (sum(exprs[setdiff(cat_markers,cell_types[['NK']])]>1)/length(setdiff(cat_markers,cell_types[['NK']])))
  if(max(scores)==1) return(names(which.max(scores)))
  if(cat=='Epithelial' & all(exprs[cat_markers] < 0))  return('Enterocyte')
  if(cat=='Stromal' & all(exprs[cat_markers] < 0))  return('Stroma')
  return('unknown')
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

spe$prior<-apply(assay(spe, 'exprs'), 2, function(x) assign_cat(x))


set.seed(68)
indices<-sample(which(spe$prior!='unknown'), 100000)
df_train<- cbind(t(as.matrix(assay(spe, 'exprs')[category_markers, indices])), spe$prior[indices])
colnames(df_train)<-c(category_markers, 'cell_type')
write.csv(df_train, 'categories_training.csv', row.names = F)

df_all<- t(as.matrix(assay(spe, 'exprs')[category_markers,]))
colnames(df_all)<-category_markers
write.csv(df_all, 'categories_all.csv', row.names = F)

system(paste0('python3 nn_categories.py ', length(category_markers)))

probs <- read.csv('probabilities.csv')
spe$cell_category<- apply(probs, 1, function(p) names(which.max(p)))

# Scaling within category
assay(spe, 'exprs_cat', withDimnames=FALSE)<-matrix(0, nrow(spe), ncol(spe))
for (cat in names(cat_type_map)){
  assay(spe, 'exprs_cat')[,spe$cell_category==cat]<-t(scale(t(assay(spe, 'counts')[,spe$cell_category==cat]), center= rowMeans(assay(spe, 'counts')[,spe$cell_category==cat]), scale=rowSds(assay(spe, 'counts')[,spe$cell_category==cat])))
}
spe$prior_type<-NA
for (cat in names(cat_type_map)){
  spe$prior_type[spe$cell_category==cat]<-apply(assay(spe, 'exprs_cat')[,spe$cell_category==cat], 2, function(x) assign_type(x, cat))
}
# Restrict stroma cells in the training dataset
set.seed(43)
spe$prior_type[sample(which(spe$prior_type=='Stroma'), sum(spe$prior_type=='Stroma')-7500)]<-'unknown'

df_train<- cbind(t(as.matrix(assay(spe, 'exprs_cat')[type_markers, spe$prior_type!='unknown'])), spe$cell_category[spe$prior_type!='unknown'], spe$prior_type[spe$prior_type!='unknown'])
colnames(df_train)<-c(type_markers, 'cell_category', 'cell_type')
write.csv(df_train, 'types_training.csv', row.names = F)

df_all<- cbind(t(as.matrix(assay(spe, 'exprs_cat')[type_markers, ])), spe$cell_category)
colnames(df_all)<- c(type_markers, 'cell_category')
write.csv(df_all, 'types_all.csv', row.names = F)

system('python3 kNN_types.py ')

spe$cell_type<-read.csv('cell_types.csv')$cell_type

# STELLAR classification
df_train<- cbind(t(as.matrix(assay(spe, 'exprs')[1:47,spe$prior_type!='unknown'])), spatialCoords(spe)[spe$prior_type!='unknown',], spe$unique_region[spe$prior_type!='unknown'], spe$prior_type[spe$prior_type!='unknown'])
colnames(df_train)<-c(rownames(spe)[1:47], 'x', 'y', 'unique_region', 'cell_type')
write.csv(df_train, 'stellar/data/types_training.csv', row.names = F)

df_all<- cbind(t(as.matrix(assay(spe, 'exprs')[1:47,])), spatialCoords(spe), spe$unique_region)
colnames(df_all)<-c(rownames(spe[1:47]), 'x', 'y', 'unique_region')
write.csv(df_all, 'stellar/data/types_all.csv', row.names = F)

setwd('./stellar')
system('python3 STELLAR_run.py --dataset types')
setwd('..')

stellar_types<-read.csv('stellar/experiments/STELLAR_run/types_STELLAR/types_results.csv', header=F)
type_ids<-read.csv('stellar/data/types.csv', row.names = 2, header=F)

spe$stellar_type <- type_ids[as.character(stellar_types[,1]), 1]
spe$stellar_category <- ifelse(spe$stellar_type %in% c(cat_type_map[['Epithelial']], 'Enterocyte'), 'Epithelial',
                               ifelse(spe$stellar_type %in% cat_type_map[['Immune']], 'Immune', 'Stromal'))  

saveRDS(spe, 'spe.rds')

# Compression & Decompression

immune <- c("CD8",    "CD3" ,   "CD4"  ,  "CD19"  , "CD21" ,  "CD138",  "CD206" , "NKG2D" , "CD7"  ,  "CD123" , "CD127",  "CD69"  , "CD45RO")
epithelial <- c( "MUC1" , "CHGA" , "MUC2" ,"CDX2"  ,"ITLN1")
stromal <-c("aSMA",  "Podoplanin" , "CD90"  )

pairs<-rbind(expand.grid(epithelial, immune), expand.grid(stromal, immune))
pairs<-t(apply(pairs, 1, function(p) if(co_exp_matrix[p[1],p[2]]<=.1) p else c(NA, NA)))
pairs<-pairs[!rowSums(is.na(pairs)),]

corr_by_pair <- min_corr_by_pair <- rep(0, nrow(pairs))
f1_by_pair <- min_f1_by_pair <-rep(0, nrow(pairs))
names(corr_by_pair) <- names(f1_by_pair) <- names(min_corr_by_pair) <- names(min_f1_by_pair) <-apply(pairs, 1, function(x) paste0(x, collapse='+'))

# Compression by pair - Basic decompression

for(i in 1:nrow(pairs)){
  p<-unlist(pairs[i,])
  if(p[1]%in% immune) {
    cat1<-'Immune'} else {
      if (p[1]%in% epithelial) cat1<-'Epithelial'
      else cat1<-'Stromal'}
  
  if(p[2]%in% immune) {
    cat2<-'Immune'} else {
      if (p[2]%in% epithelial) cat2<-'Epithelial'
      else cat2<-'Stromal'}
  # Compression
  comp_ch <-colSums(assay(spe, 'counts')[p,])
  # Decompression
  decompr_ch1<-(spe$cell_category == cat1) * comp_ch
  decompr_ch2<-(spe$cell_category == cat2) * comp_ch
  # Scaling
  decompr_ch1<-(decompr_ch1 - mean(decompr_ch1))/sd(decompr_ch1)
  decompr_ch2<-(decompr_ch2 - mean(decompr_ch2))/sd(decompr_ch2)
  
  corr_by_pair[paste0(p, collapse='+')] <- mean(c(cor(assay(spe, 'exprs')[p[1],], decompr_ch1), cor(assay(spe, 'exprs')[p[2],], decompr_ch2)))
  f1_by_pair[paste0(p, collapse='+')] <-mean(c(confusionMatrix(as.factor(decompr_ch1>1), as.factor(assay(spe, "exprs")[p[1],]>1),mode='everything', positive='TRUE')$byClass['F1'],confusionMatrix(as.factor(decompr_ch2>1), as.factor(assay(spe, "exprs")[p[2],]>1),mode='everything', positive='TRUE')$byClass['F1']))
  
}

# Compression by pair - LOESS decompression

set.seed(100)
spe_mini<-spe[,sample(seq_len(ncol(spe)), 100000)]

corr_by_pair_loess <- f1_by_pair_loess <- rep(0, nrow(pairs))
names(corr_by_pair_loess) <- names(f1_by_pair_loess) <-apply(pairs, 1, function(x) paste0(x, collapse='+'))

set.seed(436)
for(i in 1:nrow(pairs)){
  p<-unlist(pairs[i,])
  if(p[1]%in% immune) {
    cat1<-'Immune'} else {
      if (p[1]%in% epithelial) cat1<-'Epithelial'
      else cat1<-'Stromal'}
  
  if(p[2]%in% immune) {
    cat2<-'Immune'} else {
      if (p[2]%in% epithelial) cat2<-'Epithelial'
      else cat2<-'Stromal'}
  
  # Cross-validation
  
  for(k in 1:10){
    
    
    sample <- sample(c(TRUE, FALSE), ncol(spe_mini), replace=TRUE, prob=c(0.9, 0.1))
    
    spe_train <-spe_mini[,sample]
    spe_test <-spe_mini[,!sample]
    # Compression
    comp_ch_train <-colSums(assay(spe_train, 'counts')[p,])
    comp_ch_test <-colSums(assay(spe_test, 'counts')[p,])
    # Model training
    fit1 <- loess(assay(spe_train, 'counts')[p[1],spe_train$cell_category==cat1]~comp_ch_train[spe_train$cell_category==cat1], degree = 1, span = 0.1, control=loess.control(surface="direct"))
    fit2 <- loess(assay(spe_train, 'counts')[p[2],spe_train$cell_category==cat2]~comp_ch_train[spe_train$cell_category==cat2], degree = 1, span = 0.1, control=loess.control(surface="direct"))
    fit3 <- loess(assay(spe_train, 'counts')[p[1],!(spe_train$cell_category%in%c(cat1, cat2))]~comp_ch_train[!(spe_train$cell_category%in%c(cat1, cat2))], degree = 1, span = 0.1, control=loess.control(surface="direct"))
    # Decompression
    decompr_ch1<- decompr_ch2 <-rep(NA, length(comp_ch_test))
    decompr_ch1[spe_test$cell_category==cat1]<-predict(fit1, comp_ch_test[spe_test$cell_category==cat1])
    decompr_ch2[spe_test$cell_category==cat1]<-comp_ch_test[spe_test$cell_category==cat1]-decompr_ch1[spe_test$cell_category==cat1]
    decompr_ch2[spe_test$cell_category==cat2]<-predict(fit2, comp_ch_test[spe_test$cell_category==cat2])
    decompr_ch1[spe_test$cell_category==cat2]<-comp_ch_test[spe_test$cell_category==cat2]-decompr_ch2[spe_test$cell_category==cat2]
    decompr_ch1[!(spe_test$cell_category%in%c(cat1, cat2))]<-predict(fit3, comp_ch_test[!(spe_test$cell_category%in%c(cat1, cat2))])
    decompr_ch2[!(spe_test$cell_category%in%c(cat1, cat2))]<-comp_ch_test[!(spe_test$cell_category%in%c(cat1, cat2))]-decompr_ch1[!(spe_test$cell_category%in%c(cat1, cat2))]
    # Scaling
    decompr_ch1<-(decompr_ch1 - mean(decompr_ch1))/sd(decompr_ch1)
    decompr_ch2<-(decompr_ch2 - mean(decompr_ch2))/sd(decompr_ch2)
    
    corr_by_pair_loess[paste0(p, collapse='+')] <-corr_by_pair_loess[paste0(p, collapse='+')] + mean(c(cor(assay(spe_test, 'exprs')[p[1],], decompr_ch1), cor(assay(spe_test, 'exprs')[p[2],], decompr_ch2)))
    f1_by_pair_loess[paste0(p, collapse='+')] <-f1_by_pair_loess[paste0(p, collapse='+')] + mean(c(confusionMatrix(as.factor(decompr_ch1>1), as.factor(assay(spe_test, "exprs")[p[1],]>1),mode='everything', positive='TRUE')$byClass['F1'],confusionMatrix(as.factor(decompr_ch2>1), as.factor(assay(spe_test, "exprs")[p[2],]>1),mode='everything', positive='TRUE')$byClass['F1']))
  }
}

corr_by_pair_loess <- corr_by_pair_loess/10
f1_by_pair_loess <- f1_by_pair_loess/10

df<-data.frame(f1=f1_by_pair, f1_loess=f1_by_pair_loess, corr=corr_by_pair, corr_loess=corr_by_pair_loess)

write.csv(df, 'data/reconstruction_pairs.csv')

comp <- data.frame( cor = c(corr_by_pair, corr_by_pair_loess), 
                    f1 = c(f1_by_pair, f1_by_pair_loess),
                    method = c(rep('Basic', nrow(pairs)), rep('LOESS', nrow(pairs))), 
                    paired = c(1:nrow(pairs), 1:nrow(pairs))) 


write.csv(comp, 'data/reconstruction_pairs_long.csv')


# Random combinations of 8 pairs

optimal<- corr_by_pair>0.6 & f1_by_pair>0.5

pairs<-pairs[optimal,]

set.seed(49)
combinations<-comboSample(nrow(pairs), 8, n = 1e7)
combinations<- t(apply(combinations, 1, function(x) if (!any(duplicated(as.vector(pairs[x,])))) x else rep(NA, length(x) )))
combinations<-combinations[!rowSums(is.na(combinations)),]

set.seed(53)
combinations<-combinations[sample(1:nrow(combinations), 20),]

write.csv(sapply(1:nrow(combinations), function(comb) apply(pairs[combinations[comb,],], 1, function(x) paste0(x, collapse='+'))), 'data/random_combinations.csv')

# Compression by pair combination - Basic decompression

dir.create("compression_results")
for (comb in 1:nrow(combinations)) {
  
  ch_pairs<-split(t(pairs[combinations[comb,],]), rep(1:8, each = 2))
  names(ch_pairs)<-NULL
  compressed_channels<-c(setdiff(rownames(spe), list_c(ch_pairs)), ch_pairs)
  
  # Compression
  rowData<-rowData(spe)[!rowData(spe)$name%in%list_c(ch_pairs),c('channel', 'name')]
  rowData<-rbind(rowData, data.frame(channel=(length(compressed_channels)-length(ch_pairs)+1):length(compressed_channels), name= sapply(ch_pairs, function(p) paste0(p, collapse='+'))))
  rownames(rowData)<-rowData$name
  
  spe_compr<-SpatialExperiment(
    metadata=metadata(spe),
    rowData=rowData,
    colData=colData(spe),
    spatialCoords=spatialCoords(spe),
    imgData=imgData(spe)
  )
  
  assay(spe_compr, 'counts', withDimnames=FALSE) <-  t(sapply(compressed_channels, function(ch) {
    if(length(ch)==1) assay(spe, 'counts')[ch,]
    else colSums(assay(spe, 'counts')[ch,])
  }))
  
  # Decompression
  spe_decompr<-SpatialExperiment(
    metadata=metadata(spe),
    rowData=rowData(spe),
    colData=colData(spe_compr),
    spatialCoords=spatialCoords(spe_compr),
    imgData=imgData(spe_compr)
  )
  
  assay(spe_decompr, 'counts', withDimnames=FALSE)<-matrix(0, 49, ncol(spe_compr))
  assay(spe_decompr, 'counts')[setdiff(rownames(spe), list_c(ch_pairs)),] <-  assay(spe_compr, 'counts')[setdiff(rownames(spe), list_c(ch_pairs)),]
  
  for(p in 1:length(ch_pairs)){
    
    if(ch_pairs[[p]][1]%in% immune) {
      cat1<-'Immune'} else {
        if (ch_pairs[[p]][1]%in% epithelial) cat1<-'Epithelial'
        else cat1<-'Stromal'}
    
    
    if(ch_pairs[[p]][2]%in% immune) {
      cat2<-'Immune'} else {
        if (ch_pairs[[p]][2]%in% epithelial) cat2<-'Epithelial'
        else cat2<-'Stromal'}
    
    assay(spe_decompr, 'counts')[ch_pairs[[p]],]<- rbind((spe_compr$cell_category == cat1) * assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),] , (spe_compr$cell_category == cat2) * assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),] )
    
  }
  
  # Scaling
  # assay(spe_decompr, 'exprs', withDimnames=FALSE)<-t(scale(t(counts(spe_decompr)), center= rowMeans(counts(spe_decompr)), scale=rowSds(counts(spe_decompr))))
  
  assay(spe_decompr, 'exprs_cat', withDimnames=FALSE)<-matrix(0, nrow(spe_decompr), ncol(spe_decompr))
  for (cat in names(cat_type_map)){
    assay(spe_decompr, 'exprs_cat')[,spe_decompr$cell_category==cat]<-t(scale(t(assay(spe_decompr, 'counts')[,spe_decompr$cell_category==cat]), center= rowMeans(assay(spe_decompr, 'counts')[,spe_decompr$cell_category==cat]), scale=rowSds(assay(spe_decompr, 'counts')[,spe_decompr$cell_category==cat])))
  }
  
  # Cell typing
  spe_decompr$prior_type<-NA
  for (cat in names(cat_type_map)){
    spe_decompr$prior_type[spe_decompr$cell_category==cat]<-apply(assay(spe_decompr, 'exprs_cat')[,spe_decompr$cell_category==cat], 2, function(x) assign_type(x, cat))
  }
  
  set.seed(43)
  spe_decompr$prior_type[sample(which(spe_decompr$prior_type=='Stroma'), sum(spe_decompr$prior_type=='Stroma')-7500)]<-'unknown'
  #spe_decompr$prior_type[sample(which(spe_decompr$prior_type=='Enterocyte'), sum(spe_decompr$prior_type=='Enterocyte')-7500)]<-'unknown'
  
  df <- data.frame(cbind(t(assay(spe_decompr, "exprs_cat")[type_markers,spe_decompr$prior_type!='unknown']), spe_decompr$prior_type[spe_decompr$prior_type!='unknown'], spe_decompr$cell_category[spe_decompr$prior_type!='unknown']))
  colnames(df)<-c(type_markers, 'cell_type', 'cell_category')
  write.csv(df, 'types_train.csv', row.names = F)
  
  df <- cbind(t(assay(spe_decompr, "exprs_cat")[type_markers,]), spe_decompr$cell_category)
  colnames(df)<-c(type_markers, 'cell_category')
  write.csv(df, 'types_all.csv', row.names = F)
  
  system('python3 kNN_types.py')
  
  spe_decompr$cell_type<-read.csv('cell_types.csv')$cell_type
  
  saveRDS(spe_decompr, paste0("compression_results/spe_decompr_", comb, ".rds"))
  
}

typing_accuracy<-sapply(1:nrow(combinations), function(i) {
  spe_decompr<-readRDS(paste0("compression_results/spe_decompr_", i, ".rds"))
  sum(spe$cell_type==spe_decompr$cell_type)/ncol(spe)
})


spe_8<-readRDS(paste0('compression_results/spe_decompr_', which.max(typing_accuracy), '.rds'))
write.csv(apply(pairs[combinations[which.max(typing_accuracy),],], 1, function(x) paste0(x, collapse='+')), 'data/best_compression.csv')

saveRDS(spe_8, 'spe_decompr.rds')


# Compression by pair combination - LOESS decompression

typing_accuracy_loess <- rep(0, nrow(combinations))

set.seed(65)
for (comb in 1:nrow(combinations)) {
  
  ch_pairs<-split(t(pairs[combinations[comb,],]), rep(1:8, each = 2))
  names(ch_pairs)<-NULL
  compressed_channels<-c(setdiff(rownames(spe), list_c(ch_pairs)), ch_pairs)
  
  # Cross-validation
  for(k in 1:10){
    
    
    sample <- sample(c(TRUE, FALSE), ncol(spe_mini), replace=TRUE, prob=c(0.9, 0.1))
    
    spe_train <-spe_mini[,sample]
    spe_test <-spe_mini[,!sample]
    
    # Rescaling & cell typing within test dataset
    
    assay(spe_test, 'exprs_cat', withDimnames=FALSE)<-matrix(0, nrow(spe_test), ncol(spe_test))
    for (cat in names(cat_type_map)){
      assay(spe_test, 'exprs_cat')[,spe_test$cell_category==cat]<-t(scale(t(assay(spe_test, 'counts')[,spe_test$cell_category==cat]), center= rowMeans(assay(spe_test, 'counts')[,spe_test$cell_category==cat]), scale=rowSds(assay(spe_test, 'counts')[,spe_test$cell_category==cat])))
    }
    spe_test$prior_type<-NA
    for (cat in names(cat_type_map)){
      spe_test$prior_type[spe_test$cell_category==cat]<-apply(assay(spe_test, 'exprs_cat')[,spe_test$cell_category==cat], 2, function(x) assign_type(x, cat))
    }
    
    spe_test$prior_type[sample(which(spe_test$prior_type=='Stroma'), sum(spe_test$prior_type=='Stroma')-100)]<-'unknown'

    df <- data.frame(cbind(t(assay(spe_test, "exprs_cat")[type_markers,spe_test$prior_type!='unknown']), spe_test$prior_type[spe_test$prior_type!='unknown'], spe_test$cell_category[spe_test$prior_type!='unknown']))
    colnames(df)<-c(type_markers, 'cell_type', 'cell_category')
    write.csv(df, 'types_train.csv', row.names = F)
    
    df <- cbind(t(assay(spe_test, "exprs_cat")[type_markers,]), spe_test$cell_category)
    colnames(df)<-c(type_markers, 'cell_category')
    write.csv(df, 'types_all.csv', row.names = F)
    
    system('python3 kNN_types.py ')
    
    spe_test$cell_type<-read.csv('cell_types.csv')$cell_type
    
    # Compression
    rowData<-rowData(spe)[!rowData(spe)$name%in%list_c(ch_pairs),c('channel', 'name')]
    rowData<-rbind(rowData, data.frame(channel=(length(compressed_channels)-length(ch_pairs)+1):length(compressed_channels), name= sapply(ch_pairs, function(p) paste0(p, collapse='+'))))
    rownames(rowData)<-rowData$name
    
    spe_compr<-SpatialExperiment(
      metadata=metadata(spe),
      rowData=rowData,
      colData=colData(spe_test),
      spatialCoords=spatialCoords(spe),
      imgData=imgData(spe)
    )
    
    assay(spe_compr, 'counts', withDimnames=FALSE) <-  t(sapply(compressed_channels, function(ch) {
      if(length(ch)==1) assay(spe_test, 'counts')[ch,]
      else colSums(assay(spe_test, 'counts')[ch,])
    }))
    
    # Decompression
    
    spe_decompr<-SpatialExperiment(
      metadata=metadata(spe),
      rowData=rowData(spe),
      colData=colData(spe_compr),
      spatialCoords=spatialCoords(spe_compr),
      imgData=imgData(spe_compr)
    )
    
    assay(spe_decompr, 'counts', withDimnames=FALSE)<-matrix(0, 49, ncol(spe_compr))
    assay(spe_decompr, 'counts')[setdiff(rownames(spe), list_c(ch_pairs)),] <-  assay(spe_compr, 'counts')[setdiff(rownames(spe), list_c(ch_pairs)),]
    
    for(p in 1:length(ch_pairs)){
      
      if(ch_pairs[[p]][1]%in% immune) {
        cat1<-'Immune'} else {
          if (ch_pairs[[p]][1]%in% epithelial) cat1<-'Epithelial'
          else cat1<-'Stromal'}
      
      
      if(ch_pairs[[p]][2]%in% immune) {
        cat2<-'Immune'} else {
          if (ch_pairs[[p]][2]%in% epithelial) cat2<-'Epithelial'
          else cat2<-'Stromal'}
      
      fit1 <- loess(assay(spe_train, 'counts')[ch_pairs[[p]][1],spe_train$cell_category==cat1]~colSums(assay(spe_train, 'counts')[ch_pairs[[p]],])[spe_train$cell_category==cat1], degree = 1, span = .1)
      fit2 <- loess(assay(spe_train, 'counts')[ch_pairs[[p]][2],spe_train$cell_category==cat2]~colSums(assay(spe_train, 'counts')[ch_pairs[[p]],])[spe_train$cell_category==cat2], degree = 1, span = .1)
      fit3 <- loess(assay(spe_train, 'counts')[ch_pairs[[p]][1],!(spe_train$cell_category%in%c(cat1, cat2))]~colSums(assay(spe_train, 'counts')[ch_pairs[[p]],])[!(spe_train$cell_category%in%c(cat1, cat2))], degree = 1, span = .1)
      
      assay(spe_decompr, 'counts')[ch_pairs[[p]][1],spe_decompr$cell_category==cat1]<-predict(fit1, assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),spe_compr$cell_category==cat1])
      assay(spe_decompr, 'counts')[ch_pairs[[p]][2],spe_decompr$cell_category==cat1]<-assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),spe_compr$cell_category==cat1]-assay(spe_decompr, 'counts')[ch_pairs[[p]][1],spe_decompr$cell_category==cat1]
      assay(spe_decompr, 'counts')[ch_pairs[[p]][2],spe_decompr$cell_category==cat2]<-predict(fit2, assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),spe_compr$cell_category==cat2])
      assay(spe_decompr, 'counts')[ch_pairs[[p]][1],spe_decompr$cell_category==cat2]<-assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),spe_compr$cell_category==cat2]-assay(spe_decompr, 'counts')[ch_pairs[[p]][2],spe_decompr$cell_category==cat2]
      assay(spe_decompr, 'counts')[ch_pairs[[p]][1],!(spe_decompr$cell_category%in%c(cat1, cat2))]<-predict(fit3, assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),!(spe_compr$cell_category%in%c(cat1, cat2))])
      assay(spe_decompr, 'counts')[ch_pairs[[p]][2],!(spe_decompr$cell_category%in%c(cat1, cat2))]<-assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),!(spe_compr$cell_category%in%c(cat1, cat2))]-assay(spe_decompr, 'counts')[ch_pairs[[p]][1],!(spe_decompr$cell_category%in%c(cat1, cat2))]
      
    }
    
    
    assay(spe_decompr, 'counts')[is.na(assay(spe_decompr, 'counts'))]<-0
    
    #Scaling
    
    assay(spe_decompr, 'exprs_cat', withDimnames=FALSE)<-matrix(0, nrow(spe_decompr), ncol(spe_decompr))
    for (cat in names(cat_type_map)){
      assay(spe_decompr, 'exprs_cat')[,spe_decompr$cell_category==cat]<-t(scale(t(assay(spe_decompr, 'counts')[,spe_decompr$cell_category==cat]), center= rowMeans(assay(spe_decompr, 'counts')[,spe_decompr$cell_category==cat]), scale=rowSds(assay(spe_decompr, 'counts')[,spe_decompr$cell_category==cat])))
    }
    
    # Cell typing
    
    spe_decompr$prior_type<-NA
    for (cat in names(cat_type_map)){
      spe_decompr$prior_type[spe_decompr$cell_category==cat]<-apply(assay(spe_decompr, 'exprs_cat')[,spe_decompr$cell_category==cat], 2, function(x) assign_type(x, cat))
    }
    
    set.seed(43)
    spe_decompr$prior_type[sample(which(spe_decompr$prior_type=='Stroma'), sum(spe_decompr$prior_type=='Stroma')-100)]<-'unknown'

    df <- data.frame(cbind(t(assay(spe_decompr, "exprs_cat")[type_markers,spe_decompr$prior_type!='unknown']), spe_decompr$prior_type[spe_decompr$prior_type!='unknown'], spe_decompr$cell_category[spe_decompr$prior_type!='unknown']))
    colnames(df)<-c(type_markers, 'cell_type', 'cell_category')
    write.csv(df, 'types_train.csv', row.names = F)
    
    df <- cbind(t(assay(spe_decompr, "exprs_cat")[type_markers,]), spe_decompr$cell_category)
    colnames(df)<-c(type_markers, 'cell_category')
    write.csv(df, 'types_all.csv', row.names = F)
    
    system('python3 kNN_types.py')
    
    spe_decompr$cell_type<-read.csv('cell_types.csv')$cell_type
    
    typing_accuracy_loess[comb] <- typing_accuracy_loess[comb] + sum(spe_test$cell_type==spe_decompr$cell_type)/ncol(spe_test)
  }
}
typing_accuracy_loess <- typing_accuracy_loess/10


comp <- data.frame( accuracy = c(typing_accuracy, typing_accuracy_loess), 
                    method = c(rep('Basic', nrow(combinations)), rep('LOESS', nrow(combinations))), 
                    paired = c(1:nrow(combinations), 1:nrow(combinations))) 

write.csv(comp, 'data/typing_accuracy.csv')
