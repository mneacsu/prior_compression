pacman::p_load(imcRtools, CATALYST, cytomapper, BiocParallel, 
               plyr, tiff, corrplot, rlist, caret, ggrepel, purrr, ROCR,
               dittoSeq, viridis, scuttle, tidyverse, SpatialExperiment, scater,
               tibble, cvms, cowplot, dplyr, SCORPIUS, pliman)

# Spillover compensation

spillover_matrix<- read.csv("spillover_matrix.csv", row.names=1)

isotopes <- c("In113","In115","Pr141","Nd142","Nd143","Nd144","Nd145","Nd146","Sm147","Nd148",
              "Sm149","Nd150","Sm152","Eu153","Sm154","Gd156","Gd158","Tb159","Gd160","Dy161","Dy162","Dy163",
              "Dy164","Ho165","Er166","Er167","Er168","Tm169","Er170","Yb171","Yb172","Yb173","Yb174","Lu175","Yb176","Ir191","Ir193")
targets <-c("H3", "SMA", "INS", "CD38", "CD44", "PCSK2", "CD99", "CD68" , "MPO", "SLC2A1",
            "CD20", "AMY2A",  "CD3e", "PPY", "PIN", "GCG", "PDX1", "SST",
            "SYP", "KRT19" ,"CD45","FOXP3","CD45RA" ,"CD8a",
            "CA9", "IAPP", "KI-67", "NKX6-1", "pH3", "CD4",
            "CD31", "CDH", "PTPRN", "pRB", "cCASP3-cPARP1", "Ir191", "Ir193")

images <- list.files('img')

rownames(spillover_matrix) <- colnames(spillover_matrix) <- paste0(colnames(spillover_matrix), 'Di')
adapted_sm <- adaptSpillmat(as.matrix(spillover_matrix), paste0(isotopes, 'Di'),
                            isotope_list = CATALYST::isotope_list)

dir.create("steinbock/img")
for(i in 1:85){
  print(i)
  images_batch<-loadImages(paste0('img/',images[((i-1)*10+1):(i*10)][!is.na(images[((i-1)*10+1):(i*10)])]))
  channelNames(images_batch) <- paste0(isotopes, 'Di')
  images_comp <- compImage(images_batch, adapted_sm, BPPARAM = MulticoreParam())
  lapply(names(images_comp), function(x){
    writeImage(as.array(images_comp[[x]])/(2**32-1),
               paste0("steinbock/img/", x, ".tiff"),
               bits.per.sample = 32)
  })
}


system("mv img img_raw")

# Segmentation

system("./segmentation.sh")


# Cell matching

cells<-readRDS('paper_cells.rds')
spe<-  read_steinbock("steinbock", image_file = "metadata.csv", panel_file = "panel.csv", extract_imagemetadata_from = c("case", "slide", "part",	"group",	"stage"))

spe$match<-NA
for(im in unique(spe$sample_id)) {
  dist_matrix<-proxy::dist(spatialCoords(spe[,spe$sample_id==im])%>%as.data.frame, cells[cells$ImageNumber==im, c('Location_Center_X', 'Location_Center_Y')]%>%as.data.frame(), method='euclidean')
  spe$match[spe$sample_id==im]<-apply(dist_matrix, 1, function(x) {
    # Match cells from the two segmentations if the distance between centroids < 5 pixels 
    if(min(x)<5) which(cells$ImageNumber==im)[which.min(x)]
    else NA
  })
}
# Remove ambiguous matches
spe$match[duplicated(spe$match)]<-NA 

# cells$CellCat <- plyr::mapvalues(cells$CellCat, from = c('other'), to = c("stromal"))
# cells$CellType <- plyr::mapvalues(cells$CellType, from = c( 'macrophage',    'naiveTc',  'otherimmune', 'stromal'),
#                                   to = c("monocyte", "Tc", "other immune", "smooth muscle"))
spe$CellCat<-NA
spe$CellType<-NA
spe$CellCat[!is.na(spe$match)]<-cells$CellCat[spe$match[!is.na(spe$match)]]
spe$CellType[!is.na(spe$match)]<-cells$CellType[spe$match[!is.na(spe$match)]]

#  Scaling

assay(spe, 'asinh', withDimnames=FALSE)<- asinh(counts(spe)/3)
assay(spe, 'exprs', withDimnames=FALSE)<- t(scale(t(assay(spe, 'asinh')), center= rowMins(assay(spe, 'asinh')), scale=rowMaxs(assay(spe, 'asinh'))-rowMins(assay(spe, 'asinh'))))


# Co-expression analysis


thresholds_otsu <- sapply(rownames(spe), function(t)  pliman::otsu(assay(spe, "exprs")[t,]))

co_exp_matrix<-matrix(0, 37, 37)
rownames(co_exp_matrix)<- colnames(co_exp_matrix) <- rownames(spe)
for(i in 1:36) {
  i_intensities <- assay(spe, "exprs")[i,]
  for (j in (i+1):37){
    j_intensities<- assay(spe, "exprs")[j,]
    co_exp_matrix[i,j] <- co_exp_matrix[j,i] <- sum(j_intensities>thresholds_otsu[j] & i_intensities>thresholds_otsu[i])/min(sum(i_intensities>thresholds_otsu[i]), sum(j_intensities>thresholds_otsu[j]))
    
  }
}
diag(co_exp_matrix)<-1

dir.create('data')
write.csv(co_exp_matrix, 'data/co_exp_matrix.csv')


#Signal-to-noise ratio

sig_noise <- t(sapply(rownames(spe), function(x){
  # Signal-to-noise ratio
  snr <- (mean(assay(spe, 'exprs')[x,assay(spe, 'exprs')[x,]> thresholds_otsu[x]]) - mean(assay(spe, 'exprs')[x,assay(spe, 'exprs')[x,] <= thresholds_otsu[x]]))/sd(assay(spe, 'exprs')[x,assay(spe, 'exprs')[x,] <= thresholds_otsu[x]])
  # Signal intensity
  s <- mean(assay(spe, 'exprs')[x,assay(spe, 'exprs')[x,]> thresholds_otsu[x]])
  
  return(c(snr = snr, s = s))
}))

write.csv(sig_noise, 'data/snr.csv')


#Cell category assignment

cell_categories <- list(
  "islet"=c("SYP"),
  "exocrine"=c("AMY2A","KRT19"),
  "immune"=c("CD45", "MPO", "CD68"),
  "stromal"=c("SMA", "CD31")
)
category_markers<-list_c(cell_categories)

prior_frequencies<-c(
  "SYP"=.1,
  "AMY2A"=.55, 
  "KRT19"=.2,
  "CD45"=.02, 
  "MPO"=.01,
  "CD68"=.01,
  "SMA"=.02, 
  "CD31"=.02
)

thresholds_prior <- sapply(category_markers, function(t)
  quantile(assay(spe, "exprs")[t,], 1-prior_frequencies[t])
)
names(thresholds_prior)<-category_markers

assign_cat<-function(exprs){
  
  scores<-sapply(names(cell_categories), function(t)
    any(exprs[cell_categories[[t]]]>thresholds_prior[cell_categories[[t]]]) - any(exprs[setdiff(category_markers, cell_categories[[t]])]>thresholds_prior[setdiff(category_markers, cell_categories[[t]])])
  )
  if(max(scores)==1) return(names(which.max(scores)))
  else return('unknown')
}   

spe$prior<-apply(assay(spe, 'exprs'), 2, function(x) assign_cat(x))

set.seed(34)
df_train<-cbind(t(assay(spe, "exprs")[category_markers, spe$prior!='unknown']), spe$prior[spe$prior!='unknown'])
colnames(df_train)<-c(category_markers, 'cell_category')
write.csv(df_train[sample(1:nrow(df_train), 100000),], 'df_train.csv', row.names = F)

df_all<- t(assay(spe, "exprs")[category_markers,])
df_all<-data.frame(df_all)
write.csv(df_all, 'df_all.csv', row.names = F)

system(paste0('python3 nn_categories.py ', length(category_markers)))

probs <- read.csv('probabilities.csv')
spe$cell_category <- apply(probs, 1, function(p)
  if(max(p)>0.9) colnames(probs)[which.max(p)]
  else 'unknown'
)

# Cell typing

cat_type_map<-list(
  "islet"=c("alpha", "beta", "gamma", "delta"),
  "exocrine"=c("acinar", "ductal"),
  "immune"=c("Tc", "Th", "B", "monocyte", "neutrophil"),
  "stromal"=c("endothelial", "smooth muscle")
)

cell_types<-list(
  "alpha"=c("GCG", "PCSK2"),
  "beta"=c( "INS", "PIN", "NKX6-1", "IAPP", "PDX1"),
  "delta"=c("SST", "PDX1"),
  "gamma"=c("PPY"),
  "acinar"=c("AMY2A"),
  "ductal"=c("KRT19", "PDX1"),
  "Tc"=c("CD8a", "CD3e"),
  "Th"=c("CD4", "CD3e"),
  "B"=c("CD20"),
  "monocyte"=c("CD68"),
  "neutrophil"=c("MPO"),
  "endothelial"=c("CD31"),
  "smooth muscle"=c("SMA")
)

type_markers<-unique(list_c(cell_types))

assign_type<-function(exprs, cat){
  cat_markers<-unique(list_c(cell_types[cat_type_map[[cat]]]))
  scores<-sapply(cat_type_map[[cat]], function(t) sum(exprs[cell_types[[t]]]>thresholds_otsu[cell_types[[t]]])/length(cell_types[[t]]) - (sum(exprs[setdiff(cat_markers,cell_types[[t]])]>thresholds_otsu[setdiff(cat_markers,cell_types[[t]])])/length(setdiff(cat_markers,cell_types[[t]]))))
  if(max(scores)==1) return(names(which.max(scores)))
  if(cat=='immune' && all(scores==0))  return('other immune')
  return('unknown')
}


spe$prior_type<-NA
for (cat in names(cell_categories)){
  thresholds_otsu <- sapply(unique(list_c(cell_types[cat_type_map[[cat]]])), function(t) pliman::otsu(assay(spe, "exprs")[t,spe$cell_category==cat]))
  spe$prior_type[spe$cell_category==cat]<-apply(assay(spe, 'exprs')[,spe$cell_category==cat], 2, function(x) assign_type(x, cat))
}

spe$prior_type[spe$cell_category=='unknown']<- 'unknown'

set.seed(100)
df <- data.frame(cbind(t(assay(spe, "exprs")[type_markers,spe$prior_type!='unknown']), spe$prior_type[spe$prior_type!='unknown'], spe$cell_category[spe$prior_type!='unknown']))
colnames(df)<-c(type_markers, 'cell_type', 'cell_category')
write.csv(df[sample(1:nrow(df), 100000),], 'types_train.csv', row.names = F)

df <- cbind(t(assay(spe, "exprs")[type_markers,]), spe$cell_category)
colnames(df)<-c(type_markers, 'cell_category')
write.csv(df, 'types_all.csv', row.names = F)

system('python3 kNN_types.py')

spe$cell_type<-'unknown'
spe$cell_type[spe$cell_category!='unknown']<-read.csv('cell_types.csv')$cell_type[spe$cell_category!='unknown']

confusionMatrix(as.factor(spe$cell_category[!is.na(spe$match)]), as.factor(spe$CellCat[!is.na(spe$match)]))
confusionMatrix(as.factor(spe$cell_type[!is.na(spe$match)]), as.factor(spe$CellType[!is.na(spe$match)]))

saveRDS(spe, "spe.rds")

# Compression & Decompression

immune<-c("CD8a", "CD3e", "CD20", "CD4", "CD45RA", "CD38", "FOXP3")
islet<-c("GCG", "PCSK2",  "INS", "PIN", "NKX6-1", "IAPP", "PPY",  
         "SST", "SLC2A1", "PTPRN", "CD99")
exocrine<-c("CD44")

pairs<-rbind(expand.grid(immune, islet), expand.grid(exocrine, islet))
# Filter by co-expression rate
pairs<-t(apply(pairs, 1, function(p) if(co_exp_matrix[p[1],p[2]] < 0.15) p else c(NA, NA)))
pairs<-pairs[!rowSums(is.na(pairs)),]

# Simulation by pair - Basic decompression
corr_by_pair <- f1_by_pair <- rep(0, nrow(pairs))
names(corr_by_pair) <- names(f1_by_pair) <- apply(pairs, 1, function(x) paste0(x, collapse='+'))

thresholds_otsu <- sapply(rownames(spe), function(t)  pliman::otsu(assay(spe, "exprs")[t,]))
signal_intensities<-sapply(rownames(spe), function(x) mean(assay(spe, 'exprs')[x,assay(spe, 'exprs')[x,]> thresholds_otsu[x]]))

for(i in 1:nrow(pairs)){
  p<-unlist(pairs[i,])
  if(p[1]%in% immune) {
    cat1<-'immune'} else {
      if (p[1]%in% exocrine) cat1<-'exocrine'
      else cat1<-'islet'}
  
  if(p[2]%in% immune) {
    cat2<-'immune'} else {
      if (p[2]%in% exocrine) cat2<-'exocrine'
      else cat2<-'islet'}
  
  # Compression
  comp_ch <-colSums(counts(spe)[p,])
  # Decompression
  decompr_ch1<-(spe$cell_category == cat1) * comp_ch
  decompr_ch2<-(spe$cell_category == cat2) * comp_ch
  # Scaling
  decompr_ch1<-asinh(decompr_ch1/3)
  decompr_ch2<-asinh(decompr_ch2/3)
  decompr_ch1<-decompr_ch1/max(decompr_ch1)
  decompr_ch2<-decompr_ch2/max(decompr_ch2)
  
  t1<-pliman::otsu(decompr_ch1)
  t2<-pliman::otsu(decompr_ch2)

  corr_by_pair[paste0(p, collapse='+')] <- mean(c(cor(assay(spe, 'exprs')[p[1],], decompr_ch1), cor(assay(spe, 'exprs')[p[2],], decompr_ch2)))
  f1_by_pair[paste0(p, collapse='+')] <-mean(c(confusionMatrix(as.factor(decompr_ch1>t1), as.factor(assay(spe, "exprs")[p[1],]>thresholds_otsu[p[1]]),mode='everything', positive='TRUE')$byClass['F1'],
                                               confusionMatrix(as.factor(decompr_ch2>t2), as.factor(assay(spe, "exprs")[p[2],]>thresholds_otsu[p[2]]),mode='everything', positive='TRUE')$byClass['F1']))

}

# Simulation by pair - LOESS decompression

corr_by_pair_loess <- f1_by_pair_loess <- rep(0, nrow(pairs))
names(corr_by_pair_loess) <- names(f1_by_pair_loess) <-apply(pairs, 1, function(x) paste0(x, collapse='+'))

set.seed(55)
spe_mini<-spe[,sample(seq_len(ncol(spe)), 100000)]

set.seed(436)
for(i in 1:nrow(pairs)){
  p<-unlist(pairs[i,])
  if(p[1]%in% immune) {
    cat1<-'immune'} else {
      if (p[1]%in% exocrine) cat1<-'exocrine'
      else cat1<-'islet'}
  
  if(p[2]%in% immune) {
    cat2<-'immune'} else {
      if (p[2]%in% exocrine) cat2<-'exocrine'
      else cat2<-'islet'}
  
  
  # 10 x 10-fold cross-validation
  for(k in 1:10){
    
    
    sample <- sample(c(TRUE, FALSE), ncol(spe_mini), replace=TRUE, prob=c(0.9, 0.1))
    
    spe_train <-spe_mini[,sample]
    spe_test <-spe_mini[,!sample]
    
    # Compression
    
    comp_ch_train <-colSums(assay(spe_train, 'counts')[p,])
    comp_ch_test <-colSums(assay(spe_test, 'counts')[p,])
    
    thresholds_otsu<-sapply(rownames(spe), function(t) pliman::otsu(assay(spe_test, 'exprs')[t,]))
    
    # Model training
    
    fit1 <- loess(assay(spe_train, 'counts')[p[1],spe_train$cell_category==cat1]~comp_ch_train[spe_train$cell_category==cat1], degree = 1, span = 0.1, control=loess.control(surface="direct"))
    fit2 <- loess(assay(spe_train, 'counts')[p[2],spe_train$cell_category==cat2]~comp_ch_train[spe_train$cell_category==cat2], degree = 1, span = 0.1, control=loess.control(surface="direct"))
    fit3 <- loess(assay(spe_train, 'counts')[names(which.max(signal_intensities[p])),!(spe_train$cell_category%in%c(cat1, cat2))]~comp_ch_train[!(spe_train$cell_category%in%c(cat1, cat2))], degree = 1, span = 0.1, control=loess.control(surface="direct"))
    
    # Decompression
    
    decompr_ch1<- decompr_ch2 <-rep(NA, length(comp_ch_test))
    decompr_ch1[spe_test$cell_category==cat1]<-predict(fit1, comp_ch_test[spe_test$cell_category==cat1])
    decompr_ch2[spe_test$cell_category==cat1]<-comp_ch_test[spe_test$cell_category==cat1]-decompr_ch1[spe_test$cell_category==cat1]
    decompr_ch2[spe_test$cell_category==cat2]<-predict(fit2, comp_ch_test[spe_test$cell_category==cat2])
    decompr_ch1[spe_test$cell_category==cat2]<-comp_ch_test[spe_test$cell_category==cat2]-decompr_ch2[spe_test$cell_category==cat2]
    
    if(signal_intensities[p[1]]>signal_intensities[p[2]]){
      decompr_ch1[!(spe_test$cell_category%in%c(cat1, cat2))]<-predict(fit3, comp_ch_test[!(spe_test$cell_category%in%c(cat1, cat2))])
      decompr_ch2[!(spe_test$cell_category%in%c(cat1, cat2))]<-comp_ch_test[!(spe_test$cell_category%in%c(cat1, cat2))]-decompr_ch1[!(spe_test$cell_category%in%c(cat1, cat2))]
    }else{
      decompr_ch2[!(spe_test$cell_category%in%c(cat1, cat2))]<-predict(fit3, comp_ch_test[!(spe_test$cell_category%in%c(cat1, cat2))])
      decompr_ch1[!(spe_test$cell_category%in%c(cat1, cat2))]<-comp_ch_test[!(spe_test$cell_category%in%c(cat1, cat2))]-decompr_ch2[!(spe_test$cell_category%in%c(cat1, cat2))]
    }
    
    decompr_ch1[is.na(decompr_ch1)]<-0
    decompr_ch2[is.na(decompr_ch2)]<-0
    
    # Scaling
    
    decompr_ch1<-asinh(decompr_ch1/3)
    decompr_ch2<-asinh(decompr_ch2/3)
    decompr_ch1<-decompr_ch1/max(decompr_ch1)
    decompr_ch2<-decompr_ch2/max(decompr_ch2)
    
    t1<-pliman::otsu(decompr_ch1)
    t2<-pliman::otsu(decompr_ch2)

    corr_by_pair_loess[paste0(p, collapse='+')] <-corr_by_pair_loess[paste0(p, collapse='+')] + mean(c(cor(assay(spe_test, 'exprs')[p[1],], decompr_ch1), cor(assay(spe_test, 'exprs')[p[2],], decompr_ch2)))
    f1_by_pair_loess[paste0(p, collapse='+')] <-f1_by_pair_loess[paste0(p, collapse='+')] + mean(c(confusionMatrix(as.factor(decompr_ch1>t1), as.factor(assay(spe_test, "exprs")[p[1],]>thresholds_otsu[p[1]]),mode='everything', positive='TRUE')$byClass['F1'],confusionMatrix(as.factor(decompr_ch2>t2), as.factor(assay(spe_test, "exprs")[p[2],]>thresholds_otsu[p[2]]),mode='everything', positive='TRUE')$byClass['F1']))
  }
}

corr_by_pair_loess <- corr_by_pair_loess/10
f1_by_pair_loess <- f1_by_pair_loess/10

comp <- data.frame( cor = c(corr_by_pair, corr_by_pair_loess), 
                    f1 = c(f1_by_pair, f1_by_pair_loess),
                    method = c(rep('Basic', nrow(pairs)), rep('LOESS', nrow(pairs))), 
                    paired = c(1:nrow(pairs), 1:nrow(pairs))) 

write.csv(comp, 'data/reconstruction_pairs_long.csv')

# Statistics by pair

noise_intensities<-sapply(rownames(spe), function(x) mean(assay(spe, 'exprs')[x,assay(spe, 'exprs')[x,]<= thresholds_otsu[x]]))
snr <- signal_intensities/noise_intensities
frequencies<-sapply(rownames(spe), function(x) sum(assay(spe, 'exprs')[x,]> thresholds_otsu[x]))/ncol(spe)

mean_snr<-apply(pairs, 1, function(p) mean(c(snr[p[1]], snr[p[2]])))
diff_sig_intensities<-apply(pairs, 1, function(p) abs(signal_intensities[p[1]]-signal_intensities[p[2]]))
diff_noise_intensities<-apply(pairs, 1, function(p) abs(noise_intensities[p[1]]-noise_intensities[p[2]]))
co_exp_rate<-apply(pairs, 1, function(p) co_exp_matrix[p[1],p[2]])
mean_frequency<-apply(pairs, 1, function(p) mean(c(frequencies[p[1]],frequencies[p[2]])))

df<-data.frame(f1=f1_by_pair, f1_loess=f1_by_pair_loess, corr=corr_by_pair, corr_loess=corr_by_pair_loess,
                freq=mean_frequency, snr=mean_snr, co_exp_rate=co_exp_rate, 
               diff_noise=diff_noise_intensities, diff_sig=diff_sig_intensities)

write.csv(df, 'data/reconstruction_pairs.csv')


# Filtering by correlation coefficients and F1 scores

optimal<- corr_by_pair>0.5 & f1_by_pair>0.5
pairs<-pairs[optimal,]

# Random combinations of 6 pairs

combinations<-combn(1:nrow(pairs), 6)
combinations<- t(apply(combinations, 2, function(x) if (!any(duplicated(as.vector(pairs[x,])))) x else rep(NA, length(x) )))
combinations<-combinations[!rowSums(is.na(combinations)),]

set.seed(55)
combinations<-combinations[sample(1:nrow(combinations), 20),]

write.csv(sapply(1:nrow(combinations), function(comb) apply(pairs[combinations[comb,],], 1, function(x) paste0(x, collapse='+'))), 'data/random_combinations.csv')

# Simulation by pair combination - Basic decompression

dir.create("compression_results")
for (comb in 1:nrow(combinations)) {
  
  ch_pairs<-split(t(pairs[combinations[comb,],]), rep(1:6, each = 2))
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
  
  assay(spe_decompr, 'counts', withDimnames=FALSE)<-matrix(0, 37, ncol(spe_compr))
  assay(spe_decompr, 'counts')[setdiff(rownames(spe), list_c(ch_pairs)),] <-  assay(spe_compr, 'counts')[setdiff(rownames(spe), list_c(ch_pairs)),]
  
  for(p in 1:length(ch_pairs)){
    
    if(ch_pairs[[p]][1]%in% immune) {
      cat1<-'immune'} else {
        if (ch_pairs[[p]][1]%in% exocrine) cat1<-'exocrine'
        else cat1<-'islet'}
    
    
    if(ch_pairs[[p]][2]%in% immune) {
      cat2<-'immune'} else {
        if (ch_pairs[[p]][2]%in% exocrine) cat2<-'exocrine'
        else cat2<-'islet'}
    
    assay(spe_decompr, 'counts')[ch_pairs[[p]],]<- rbind((spe_compr$cell_category == cat1) * assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),] , (spe_compr$cell_category == cat2) * assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),] )
    
  }
  
  # Scaling
  
  assay(spe_decompr, "asinh") <- asinh(counts(spe_decompr)/3)
  assay(spe_decompr, 'exprs', withDimnames=FALSE)<-t(scale(t(assay(spe_decompr, 'asinh')), center= rowMins(assay(spe_decompr, 'asinh')), scale=rowMaxs(assay(spe_decompr, 'asinh'))-rowMins(assay(spe_decompr, 'asinh'))))
  
  # Cell typing
  
  spe_decompr$prior_type<-NA
  for (cat in names(cell_categories)){
    thresholds_otsu <- sapply(unique(list_c(cell_types[cat_type_map[[cat]]])), function(t){
      
      pliman::otsu(assay(spe_decompr, "exprs")[t,spe_decompr$cell_category==cat] )
    })
    spe_decompr$prior_type[spe_decompr$cell_category==cat]<-apply(assay(spe_decompr, 'exprs')[,spe_decompr$cell_category==cat], 2, function(x) assign_type(x, cat))
    
  }
  
  spe_decompr$prior_type[spe_decompr$cell_category=='unknown']<- 'unknown'
  
  df <- data.frame(cbind(t(assay(spe_decompr, "exprs")[type_markers,spe_decompr$prior_type!='unknown']), spe_decompr$prior_type[spe_decompr$prior_type!='unknown'], spe_decompr$cell_category[spe_decompr$prior_type!='unknown']))
  colnames(df)<-c(type_markers, 'cell_type', 'cell_category')
  write.csv(df, 'types_train.csv', row.names = F)
  
  df <- cbind(t(assay(spe_decompr, "exprs")[type_markers,]), spe_decompr$cell_category)
  colnames(df)<-c(type_markers, 'cell_category')
  write.csv(df, 'types_all.csv', row.names = F)
  
  system('python3 kNN_types.py')
  
  spe_decompr$cell_type<-'unknown'
  spe_decompr$cell_type[spe_decompr$cell_category!='unknown']<-read.csv('cell_types.csv')$cell_type[spe_decompr$cell_category!='unknown']
  
  saveRDS(spe_decompr, paste0("compression_results/spe_decompr_", comb, ".rds"))
  
}

typing_accuracy<-sapply(1:nrow(combinations), function(i) {
  spe_decompr<-readRDS(paste0("compression_results/spe_decompr_", i, ".rds"))
  sum(spe$cell_type==spe_decompr$cell_type)/ncol(spe)
})

# Simulation by pair combination - LOESS decompression

typing_accuracy_loess <- rep(0, nrow(combinations))

set.seed(65)
for (comb in 1:nrow(combinations)) {
  
  ch_pairs<-split(t(pairs[combinations[comb,],]), rep(1:6, each = 2))
  names(ch_pairs)<-NULL
  compressed_channels<-c(setdiff(rownames(spe), list_c(ch_pairs)), ch_pairs)
  
  # 10 x 10-fold cross-validation
  
  for(k in 1:10){
    
    
    sample <- sample(c(TRUE, FALSE), ncol(spe_mini), replace=TRUE, prob=c(0.9, 0.1))
    
    spe_train <-spe_mini[,sample]
    spe_test <-spe_mini[,!sample]
    
    # Re-scaling and re-typing within test dataset to ensure comparability
    
    assay(spe_test, 'exprs', withDimnames=FALSE)<-t(scale(t(assay(spe_test, 'asinh')), center= rowMins(assay(spe_test, 'asinh')), scale=rowMaxs(assay(spe_test, 'asinh'))-rowMins(assay(spe_test, 'asinh'))))
    
    spe_test$prior_type<-NA
    for (cat in names(cell_categories)){
      thresholds_otsu <- sapply(unique(list_c(cell_types[cat_type_map[[cat]]])), function(t) pliman::otsu(assay(spe_test, "exprs")[t,spe_test$cell_category==cat]))
      spe_test$prior_type[spe_test$cell_category==cat]<-apply(assay(spe_test, 'exprs')[,spe_test$cell_category==cat], 2, function(x) assign_type(x, cat))
    }
    
    spe_test$prior_type[spe_test$cell_category=='unknown']<- 'unknown'
    
    df <- data.frame(cbind(t(assay(spe_test, "exprs")[type_markers,spe_test$prior_type!='unknown']), spe_test$prior_type[spe_test$prior_type!='unknown'], spe_test$cell_category[spe_test$prior_type!='unknown']))
    colnames(df)<-c(type_markers, 'cell_type', 'cell_category')
    write.csv(df, 'types_train.csv', row.names = F)
    
    df <- cbind(t(assay(spe_test, "exprs")[type_markers,]), spe_test$cell_category)
    colnames(df)<-c(type_markers, 'cell_category')
    write.csv(df, 'types_all.csv', row.names = F)
    
    system('python3 kNN_types.py')
    
    spe_test$cell_type<-'unknown'
    spe_test$cell_type[spe_test$cell_category!='unknown']<-read.csv('cell_types.csv')$cell_type[spe_test$cell_category!='unknown']
    
    # Compression
    
    rowData<-rowData(spe)[!rowData(spe)$name%in%list_c(ch_pairs),c('channel', 'name')]
    rowData<-rbind(rowData, data.frame(channel=(length(compressed_channels)-length(ch_pairs)+1):length(compressed_channels), name= sapply(ch_pairs, function(p) paste0(p, collapse='+'))))
    rownames(rowData)<-rowData$name
    
    spe_compr<-SpatialExperiment(
      metadata=metadata(spe),
      rowData=rowData,
      colData=colData(spe_test),
      spatialCoords=spatialCoords(spe_test),
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
      imgData=imgData(spe)
    )
    
    assay(spe_decompr, 'counts', withDimnames=FALSE)<-matrix(0, 37, ncol(spe_compr))
    assay(spe_decompr, 'counts')[setdiff(rownames(spe), list_c(ch_pairs)),] <-  assay(spe_compr, 'counts')[setdiff(rownames(spe), list_c(ch_pairs)),]
    
    
    for(p in 1:length(ch_pairs)){
      
      if(ch_pairs[[p]][1]%in% immune) {
        cat1<-'immune'} else {
          if (ch_pairs[[p]][1]%in% exocrine) cat1<-'exocrine'
          else cat1<-'islet'}
      
      
      if(ch_pairs[[p]][2]%in% immune) {
        cat2<-'immune'} else {
          if (ch_pairs[[p]][2]%in% exocrine) cat2<-'exocrine'
          else cat2<-'islet'}
      
      fit1 <- loess(assay(spe_train, 'counts')[ch_pairs[[p]][1],spe_train$cell_category==cat1]~colSums(assay(spe_train, 'counts')[ch_pairs[[p]],])[spe_train$cell_category==cat1], degree = 1, span = 0.1)
      fit2 <- loess(assay(spe_train, 'counts')[ch_pairs[[p]][2],spe_train$cell_category==cat2]~colSums(assay(spe_train, 'counts')[ch_pairs[[p]],])[spe_train$cell_category==cat2], degree = 1, span = 0.1)
      fit3 <- loess(assay(spe_train, 'counts')[names(which.max(signal_intensities[ch_pairs[[p]]])),!(spe_train$cell_category%in%c(cat1, cat2))]~colSums(assay(spe_train, 'counts')[ch_pairs[[p]],])[!(spe_train$cell_category%in%c(cat1, cat2))], degree = 1, span = 0.1)
      
      assay(spe_decompr, 'counts')[ch_pairs[[p]][1],spe_decompr$cell_category==cat1]<-predict(fit1, assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),spe_compr$cell_category==cat1])
      assay(spe_decompr, 'counts')[ch_pairs[[p]][2],spe_decompr$cell_category==cat1]<-assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),spe_compr$cell_category==cat1]-assay(spe_decompr, 'counts')[ch_pairs[[p]][1],spe_decompr$cell_category==cat1]
      assay(spe_decompr, 'counts')[ch_pairs[[p]][2],spe_decompr$cell_category==cat2]<-predict(fit2, assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),spe_compr$cell_category==cat2])
      assay(spe_decompr, 'counts')[ch_pairs[[p]][1],spe_decompr$cell_category==cat2]<-assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),spe_compr$cell_category==cat2]-assay(spe_decompr, 'counts')[ch_pairs[[p]][2],spe_decompr$cell_category==cat2]
      if(signal_intensities[ch_pairs[[p]][1]]>signal_intensities[ch_pairs[[p]][2]]){
        assay(spe_decompr, 'counts')[ch_pairs[[p]][1],!(spe_decompr$cell_category%in%c(cat1, cat2))]<-predict(fit3, assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),!(spe_compr$cell_category%in%c(cat1, cat2))])
        assay(spe_decompr, 'counts')[ch_pairs[[p]][2],!(spe_decompr$cell_category%in%c(cat1, cat2))]<-assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),!(spe_compr$cell_category%in%c(cat1, cat2))]-assay(spe_decompr, 'counts')[ch_pairs[[p]][1],!(spe_decompr$cell_category%in%c(cat1, cat2))]
      }else{
        assay(spe_decompr, 'counts')[ch_pairs[[p]][2],!(spe_decompr$cell_category%in%c(cat1, cat2))]<-predict(fit3, assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),!(spe_compr$cell_category%in%c(cat1, cat2))])
        assay(spe_decompr, 'counts')[ch_pairs[[p]][1],!(spe_decompr$cell_category%in%c(cat1, cat2))]<-assay(spe_compr, 'counts')[paste0(ch_pairs[[p]], collapse='+'),!(spe_compr$cell_category%in%c(cat1, cat2))]-assay(spe_decompr, 'counts')[ch_pairs[[p]][2],!(spe_decompr$cell_category%in%c(cat1, cat2))]
      }
    }
    
    assay(spe_decompr, 'counts')[is.na(assay(spe_decompr, 'counts'))]<-0
    
    #Scaling
    
    assay(spe_decompr, "asinh") <- asinh(counts(spe_decompr)/3)
    assay(spe_decompr, 'exprs', withDimnames=FALSE)<-t(scale(t(assay(spe_decompr, 'asinh')), center= rowMins(assay(spe_decompr, 'asinh')), scale=rowMaxs(assay(spe_decompr, 'asinh'))-rowMins(assay(spe_decompr, 'asinh'))))
    
    # Cell typing
    
    spe_decompr$prior_type<-NA
    for (cat in names(cell_categories)){
      thresholds_otsu <- sapply(unique(list_c(cell_types[cat_type_map[[cat]]])), function(t) pliman::otsu(assay(spe_decompr, "exprs")[t,spe_decompr$cell_category==cat]))
      spe_decompr$prior_type[spe_decompr$cell_category==cat]<-apply(assay(spe_decompr, 'exprs')[,spe_decompr$cell_category==cat], 2, function(x) assign_type(x, cat))
    }
    
    spe_decompr$prior_type[spe_decompr$cell_category=='unknown']<- 'unknown'
    
    df <- data.frame(cbind(t(assay(spe_decompr, "exprs")[type_markers,spe_decompr$prior_type!='unknown']), spe_decompr$prior_type[spe_decompr$prior_type!='unknown'], spe_decompr$cell_category[spe_decompr$prior_type!='unknown']))
    colnames(df)<-c(type_markers, 'cell_type', 'cell_category')
    write.csv(df, 'types_train.csv', row.names = F)
    
    df <- cbind(t(assay(spe_decompr, "exprs")[type_markers,]), spe_decompr$cell_category)
    colnames(df)<-c(type_markers, 'cell_category')
    write.csv(df, 'types_all.csv', row.names = F)
    
    system('python3 kNN_types.py')
    
    spe_decompr$cell_type<-'unknown'
    spe_decompr$cell_type[spe_decompr$cell_category!='unknown']<-read.csv('cell_types.csv')$cell_type[spe_decompr$cell_category!='unknown']
    
    typing_accuracy_loess[comb] <- typing_accuracy_loess[comb] + sum(spe_test$cell_type==spe_decompr$cell_type)/ncol(spe_test)
  } 
}

typing_accuracy_loess <- typing_accuracy_loess/10


comp <- data.frame( accuracy = c(typing_accuracy, typing_accuracy_loess), 
                    method = c(rep('Basic', nrow(combinations)), rep('LOESS', nrow(combinations))), 
                    paired = c(1:nrow(combinations), 1:nrow(combinations))) 


write.csv(comp, 'data/typing_accuracy.csv')

spe_6<-readRDS(paste0('compression_results/spe_decompr_', which.max(typing_accuracy), '.rds'))
saveRDS(spe_6, 'spe_decompr.rds')
write.csv(apply(pairs[combinations[which.max(typing_accuracy),],], 1, function(x) paste0(x, collapse='+')), 'best_compression.csv')