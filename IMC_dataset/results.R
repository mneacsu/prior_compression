pacman::p_load(SpatialExperiment, cowplot, tidyverse, caret, cvms, ggimage, ggrepel, ROCR, rstatix)

# Co-expression analysis
co_exp_matrix<-read.csv('data/co_exp_matrix.csv', row.names=1, check.names=F)

functional<-c("SYP", "CD99", "PTPRN", "SLC2A1", "GCG", "PCSK2",  "INS", "PIN", "NKX6-1", "IAPP",  
              "SST", "PPY", "PDX1",  "AMY2A", "KRT19", "CD44", "CD45", "CD45RA", "CD68", "MPO",
              "CD20", "CD3e", "CD4", "CD8a", "FOXP3", "CD38", "CD31", "SMA")
co_exp_matrix<-co_exp_matrix[functional, functional]

pdf('plots/co_oc_mat.pdf')
corrplot(as.matrix(co_exp_matrix), method="shade", addCoef.col = 1, tl.cex = 0.8, number.cex = 0.5, type="upper",
         is.corr=F, col=colorRampPalette(c("blue","white","red"))(100))
dev.off()

prior<-read.csv("prior.csv", check.names = F, row.names = 1)[functional, functional]

pred <- prediction(c(co_exp_matrix[upper.tri(co_exp_matrix)]),c(prior[upper.tri(prior)]))
perf <- performance(pred,"tpr","fpr")
(auc<-performance(pred,"auc")@y.values[[1]])

snr_tbl<-read.csv('data/snr.csv', row.names=1, check.names=F)[functional,]
pdf('plots/snr.pdf')
ggplot(snr_tbl, aes(log2(s), log2(snr), label = rownames(snr_tbl))) +
  geom_point() +
  theme_minimal(base_size = 15) + ylab("Signal-to-noise ratio [log2]") +
  xlab("Signal intensity [log2]")+
  geom_label_repel(force=7, max.overlaps = 20)
dev.off()

# Cell typing
spe<-readRDS('spe.rds')

category_markers <- c('SYP', 'AMY2A', 'KRT19', 'CD45', 'CD68', 'MPO', 'SMA', 'CD31')
cell_types <- c("alpha" ,        "beta"    ,      "delta"      ,   "gamma"    ,  "acinar"    ,    "ductal"   ,     
                "Tc"     ,       "Th"  ,    "B"    ,    "monocyte"  ,    "neutrophil"   , "other immune",
                "endothelial"  , "smooth muscle", "unknown")

set.seed(55)
spe_mini<-spe[,sample(seq_len(ncol(spe)), 100000)]
spe_mini<-runUMAP(spe_mini, exprs_values = "exprs", subset_row = category_markers)
colnames(spe_mini) <- paste0('cell', seq_len(ncol(spe_mini)))
saveRDS(spe_mini, 'spe_mini.rds')


png('plots/UMAP_categories.png', height=1000, width = 2200)
p<-dittoDimPlot(spe_mini, 
                var = "cell_category", 
                reduction.use = "UMAP", 
                size = 0.2,
                labels.size =10, 
                legend.size = 20,
                do.label = TRUE,
                legend.show=TRUE) +
  ggtitle("Neural Network")+
  theme(plot.title = element_text(size = 40), legend.title = element_text('Cell Type'), legend.text = element_text(size = 30), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30))
q<-dittoDimPlot(spe_mini, 
                var = "CellCat", 
                reduction.use = "UMAP", 
                size = 0.2,
                labels.size =10, 
                do.label = TRUE,
                legend.show=FALSE) +
  ggtitle("Ilastik (Damond et al., 2009)")+
  theme(plot.title = element_text(size = 40), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30))
plot_grid(q,p, rel_widths = c(1, 1.2), scale = 0.9)
dev.off()

celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$cell_type, 
                                      statistics = "mean",
                                      use.assay.type = "exprs")

category_colors <- setNames(brewer.pal(length(unique(spe$cell_category)), name = "Set1"), 
                       unique(spe$cell_category))

pdf('plots/cell_types.pdf')
dittoHeatmap(celltype_mean,
             genes=c("SYP", "GCG", "PCSK2",  "INS", "PIN", "NKX6-1", "IAPP",  
                     "SST", "PPY", "PDX1",  "AMY2A", "KRT19", "CD45", "CD68", "MPO",
                     "CD20", "CD3e", "CD4", "CD8a", "CD31", "SMA"),
             assay = "exprs", 
             cluster_cols = F, 
             cluster_rows= F,
             cells.use = cell_types,
             scaled.to.max = TRUE,
             heatmap.colors.max.scaled = plasma(100),
             annot.by = c('cell_type', 'cell_category'),
             annotation_colors = list(cell_category= setNames(c("#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2"), c('islet', 'exocrine','immune','stromal','unknown')),
                                      cell_type=setNames(c("#6A3D9A", "#FF7F00", "skyblue2",  "palegreen2","#CAB2D6", "steelblue4","#FDBF6F",
                                                           "darkorange4", "gray70", "gold","darkturquoise", "beige", "royalblue","pink"), cell_types ))
             )
dev.off()
 

# Compression
spe_6 <-readRDS('spe_decompr.rds')

conf_mat6<-confusionMatrix(as.factor(spe_6$cell_type), as.factor(spe$cell_type), mode="prec_recall")
tbl6_immune<-as.tibble(conf_mat6$table[c("Tc", "Th", "B", "monocyte", "neutrophil", "other immune"), c("Tc", "Th", "B", "monocyte", "neutrophil", "other immune")])
a<-plot_confusion_matrix(tbl6_immune, 
                         target_col = "Reference", 
                         prediction_col = "Prediction",
                         class_order = c('other immune', 'monocyte', 'neutrophil', 'Th', 'Tc', 'B'),
                         rm_zero_text =F, rm_zero_percentages = F, add_zero_shading = F,
                         add_arrows = F, 
                         font_counts = font(size=5), font_row_percentages = font(size=3), font_col_percentages = font(size=3),
                         counts_col = "n", add_normalized = FALSE,  rotate_y_text = F)+ 
  theme(axis.title.x.top=element_text(size=20), axis.title.y=element_text(size=20), axis.text.x = element_text(angle = 45, size = 15, hjust=-.05), axis.text.y = element_text(size = 15))+
  labs(x='Uncompressed', y='Compressed')


tbl6_islet<-as.tibble(conf_mat6$table[c('gamma', 'delta', 'beta', 'alpha'), c('gamma', 'delta', 'beta', 'alpha')])
b<-plot_confusion_matrix(tbl6_islet, 
                         target_col = "Reference", 
                         prediction_col = "Prediction",
                         class_order = c('gamma', 'delta', 'beta', 'alpha'),
                         add_arrows = F, 
                         font_counts = font(size=5), font_row_percentages = font(size=3), font_col_percentages = font(size=3),
                         counts_col = "n", add_normalized = FALSE, rotate_y_text = F)+ 
  theme(axis.title.x.top=element_text(size=20), axis.title.y=element_text(size=20), axis.text.x = element_text(angle = 45, size = 15, hjust=-.05), axis.text.y = element_text(size = 15))+
  labs(x='Uncompressed', y='Compressed')

png('plots/confusion_matrix.png', width=1200, height = 600)
plot_grid(b, a, align='hv', scale=.9)
dev.off()

df <- read.csv('data/reconstruction_pairs.csv', row.names=1)
df$co_exp_rate_bin <- ifelse(df$co_exp_rate <=.05, '<= 0.05', '> 0.05')
t.test(df$f1~df$co_exp_rate_bin)
t.test(df$corr~df$co_exp_rate_bin)

comp<-read.csv('data/reconstruction_pairs_long.csv', row.names = 1)
comp_type<-read.csv('data/typing_accuracy.csv', row.names = 1)

png('plots/boxplots.png', width = 750, height = 250)
g1<-ggplot(comp[!(comp$paired%in%comp$paired[is.na(comp$f1)]),], aes(method, f1, fill=method)) + 
  geom_boxplot()+ 
  geom_line(aes(group = paired), color='darkgray', alpha=0.5)+ 
  geom_point(aes(fill=method,group=paired)) +
  xlab('Decompression method')+
  ylab('F1 score')+
  theme(legend.position = 'None')+
  stat_compare_means(method = "t.test", label.x = 0.75, paired=T)
g2<-ggplot(comp, aes(method, cor, fill=method)) + 
  geom_boxplot()+ 
  geom_line(aes(group = paired), color='darkgray', alpha=0.5)+ 
  geom_point(aes(fill=method,group=paired))+
  xlab('Decompression method')+
  ylab('Correlation coefficient')+
  theme(legend.position = 'None')+
  stat_compare_means(method = "t.test", label.x = 0.75, paired=T)
g3<-ggplot(comp_type, aes(method, accuracy, fill=method)) + 
  geom_boxplot()+ 
  geom_line(aes(group = paired), color='darkgray', alpha=0.5)+ 
  geom_point(aes(fill=method,group=paired))+
  xlab('Decompression method')+
  ylab('Cell typing accuracy')+
  theme(legend.position = 'None')+
  stat_compare_means(method = "t.test", label.x = 0.75, paired=T)
plot_grid(g1, g2, g3, nrow=1, labels='AUTO', scale=0.95)
dev.off()

# Downstream effects
image_mean <- aggregateAcrossCells(as(spe[,spe$cell_category=='islet'], "SingleCellExperiment"),  
                                   ids = spe[,spe$cell_category=='islet']$sample_id, 
                                   statistics = "mean",
                                   use.assay.type = "exprs")
image_mean <- image_mean[,image_mean$part!='Head'] 

image_means<-as.data.frame(t(assay(image_mean, 'exprs')))

set.seed(30)
space <- reduce_dimensionality(as.matrix(image_means), "pearson")

traj <- infer_trajectory(space)
traj$time<-1-traj$time

image_mean$pseudotime<-traj$time

saveRDS(image_mean, 'data/image_means_uncompressed.rds' )

pdf('plots/markers_islets_uncompressed.pdf')
dittoHeatmap(image_mean,
             genes=c("GCG", "PCSK2", "INS", "PIN", "NKX6-1", "IAPP", "PDX1", "SST", "PPY"),
             assay = "exprs", 
             cluster_cols = F, 
             order.by = "pseudotime",
             scaled.to.max = T,
             heatmap.colors.max.scaled = plasma(100),
             annot.by = c("pseudotime", "stage"),
             annotation_colors = list(stage= setNames(c("skyblue2", "#6A3D9A", "#FF7F00"), c('Non-diabetic', 'Onset', 'Long-duration'))))

dev.off()



image_mean_decompr <- aggregateAcrossCells(as(spe_6[,spe_6$cell_category=='islet'], "SingleCellExperiment"),  
                                   ids = spe_6[,spe_6$cell_category=='islet']$sample_id, 
                                   statistics = "mean",
                                   use.assay.type = "exprs")
image_mean_decompr <- image_mean_decompr[,image_mean_decompr$part!='Head'] 

image_means<-as.data.frame(t(assay(image_mean_decompr, 'exprs')))

set.seed(30)
space <- reduce_dimensionality(as.matrix(image_means), "pearson")

traj <- infer_trajectory(space)
traj$time<-1-traj$time

image_mean_decompr$pseudotime<-traj$time

saveRDS(image_mean_decompr, 'data/image_means_compressed.rds' )

df<- cbind(t(assay(image_mean_decompr, 'exprs')[c('PIN', 'INS', 'IAPP', 'NKX6-1'),]), as.character(image_mean_decompr$stage))
df<-as.data.frame(df)
colnames(df)[4]<-'NKX6.1'
colnames(df)[5]<-'stage'

df$stage<-factor(df$stage, levels=c("Non-diabetic", "Onset", "Long-duration"))
df$PIN<-as.numeric(df$PIN)
df$INS<-as.numeric(df$PIN)
df$IAPP<-as.numeric(df$PIN)
df$NKX6.1<-as.numeric(df$NKX6.1)

g1<-ggplot(df, aes(y=PIN, x=stage))+
  geom_boxplot(aes(fill=stage))+
  xlab("Stage")+
  ylab("PIN")+
  scale_fill_discrete(name="Stage")+
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list( c("Non-diabetic", "Onset"), c("Onset", "Long-duration"), c("Non-diabetic", "Long-duration")), label.x = 1.5)+
  theme(legend.position='None')
g2<-ggplot(df, aes(y=INS, x=stage))+
  geom_boxplot(aes(fill=stage))+
  xlab("Stage")+
  ylab("INS")+
  scale_fill_discrete(name="Stage")+
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list( c("Non-diabetic", "Onset"), c("Onset", "Long-duration"), c("Non-diabetic", "Long-duration")), label.x = 1.5)+
  theme(legend.position='None')
g3<-ggplot(df, aes(y=IAPP, x=stage))+
  geom_boxplot(aes(fill=stage))+
  xlab("Stage")+
  ylab("IAPP")+
  scale_fill_discrete(name="Stage")+
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list( c("Non-diabetic", "Onset"), c("Onset", "Long-duration"), c("Non-diabetic", "Long-duration")), label.x = 1.5)+
  theme(legend.position='None')
g4<-ggplot(df, aes(y=NKX6.1, x=stage))+
  geom_boxplot(aes(fill=stage))+
  xlab("Stage")+
  ylab("NKX6-1")+
  scale_fill_discrete(name="Stage")+
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list( c("Non-diabetic", "Onset"), c("Onset", "Long-duration"), c("Non-diabetic", "Long-duration")), label.x = 1.5)+
  theme(legend.position='None')

pdf('plots/boxplots_markers.pdf')
plot_grid(g1, g2, g3, g4, nrow=2)
dev.off()

pdf('plots/markers_islets_compressed.pdf')
dittoHeatmap(image_mean_decompr,
             genes=unique(list_c(cell_types[cat_type_map[['islet']]])),
             assay = "exprs", 
             cluster_cols = F, 
             order.by = "pseudotime",
             scaled.to.max = T,
             heatmap.colors.max.scaled = plasma(100),
             annot.by = c("pseudotime", "stage"),
             annotation_colors = list(stage= setNames(c("skyblue2", "#6A3D9A", "#FF7F00"), c('Non-diabetic', 'Onset', 'Long-duration'))))
dev.off()

cor(image_mean$pseudotime, image_mean_decompr$pseudotime)

png('plots/pseudotime.png', height=500, width=625)
ggplot(data.frame(pt_uncomp = image_mean$pseudotime, pt_comp = image_mean_decompr$pseudotime, stage = image_mean$stage) ,
       aes(x = pt_uncomp, y= pt_comp, color  = stage))+
  geom_point()+
  scale_color_manual(name='Stage', values=c("skyblue2", "#6A3D9A", "#FF7F00"))+
  xlab('Pseudotime (without compression)')+
  ylab('Pseudotime (with compression)')+
  theme(text=element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()


spe$stage<-factor(spe$stage, levels=c('Non-diabetic', 'Onset', 'Long-duration'))
spe_6$stage<-factor(spe_6$stage, levels=c('Non-diabetic', 'Onset', 'Long-duration'))

df <- as.data.frame(colData(spe[,!is.na(spe$match)&spe$CellCat=='islet'])) %>%
  select(case,stage, CellType) %>%
  arrange(stage)

p1<-ggplot(df) +
  geom_bar(aes(x=case, fill=factor(CellType, levels = c('gamma', 'delta', 'alpha', 'beta'))), position="fill")+
  facet_wrap(~stage, scales = "free_x", nrow = 1)+
  scale_fill_manual(name = "Cell Type", values = c("#CAB2D6", "darkseagreen1", "steelblue4", "bisque")) +
  theme(panel.margin = unit(0, "lines"), 
        strip.background = element_blank(),
        legend.position = "none")+
  xlab('Case')+
  ylab('Cell fraction')

df <- as.data.frame(colData(spe[,spe$cell_category=='islet'])) %>%
  select(case,stage, cell_type) %>%
  arrange(stage)

p2<-ggplot(df) +
  geom_bar(aes(x=case, fill=factor(cell_type, levels = c('gamma', 'delta', 'alpha', 'beta'))), position="fill")+
  facet_wrap(~stage, scales = "free_x", nrow = 1)+
  scale_fill_manual(name = "Cell Type", values = c("#CAB2D6", "darkseagreen1", "steelblue4", "bisque")) +
  theme(panel.margin = unit(0, "lines"), 
        strip.background = element_blank(),
        legend.position = "none")+
  xlab('Case')+
  ylab('Cell fraction')

df <- as.data.frame(colData(spe_6[,spe_6$cell_category=='islet'])) %>%
  select(case,stage, cell_type) %>%
  arrange(stage)

p3<-ggplot(df) +
  geom_bar(aes(x=case, fill=factor(cell_type, levels = c('gamma', 'delta', 'alpha', 'beta'))), position="fill")+
  facet_wrap(~stage, scales = "free_x", nrow = 1)+
  scale_fill_manual(name = "Cell Type", values = c("#CAB2D6", "darkseagreen1", "steelblue4", "bisque")) +
  theme(panel.margin = unit(0, "lines"), 
        strip.background = element_blank(),
        legend.position = "none")+
  xlab('Case')+
  ylab('Cell fraction')


df <- as.data.frame(colData(spe[,!is.na(spe$match)&spe$CellCat=='immune'])) %>%
  select(case,stage, CellType) %>%
  arrange(stage)

q1 <- ggplot(df) +
  geom_bar(aes(y=factor(CellType, levels=c('other immune','neutrophil','monocyte','B','Th' ,'Tc')), fill=stage), position="fill")+
  ylab('Cell Type')+
  xlab('Percentage')+
  scale_fill_manual(name = "Stage", values=c("skyblue2", "#6A3D9A", "#FF7F00")) +
  theme(legend.position = "none")

df <- as.data.frame(colData(spe[,spe$cell_category=='immune'])) %>%
  select(case,stage, cell_type) %>%
  arrange(stage)

q2 <- ggplot(df) +
  geom_bar(aes(y=factor(cell_type, levels=c('other immune','neutrophil','monocyte','B','Th' ,'Tc')), fill=stage), position="fill")+
  ylab('Cell Type')+
  xlab('Percentage')+
  scale_fill_manual(name = "Stage", values=c("skyblue2", "#6A3D9A", "#FF7F00")) +
  theme(legend.position = "none")

df <- as.data.frame(colData(spe_6[,spe$cell_category=='immune'])) %>%
  select(case,stage, cell_type) %>%
  arrange(stage)

q3 <- ggplot(df) +
  geom_bar(aes(y=factor(cell_type, levels=c('other immune','neutrophil','monocyte','B','Th' ,'Tc')), fill=stage), position="fill")+
  ylab('Cell Type')+
  xlab('Percentage')+
  scale_fill_manual(name = "Stage", values=c(c("skyblue2", "#6A3D9A", "#FF7F00"))) +
  theme(legend.position = "none")


row1 <- plot_grid(p1, q1, rel_widths = c(3,2), scale = 0.9)

title1 <- ggdraw() + 
  draw_label(
    "Damond et al., 2009",
    fontface = 'bold', x=.01, hjust=0
  )

row2 <- plot_grid(p2, q2, rel_widths = c(3,2), scale = 0.9)

title2 <- ggdraw() + 
  draw_label(
    "Uncompressed",
    fontface = 'bold', x=.01, hjust=0
  )

row3 <- plot_grid(p3, q3, rel_widths = c(3,2), scale = 0.9)

title3 <- ggdraw() + 
  draw_label(
    "Compressed",
    fontface = 'bold', x=.01, hjust=0
  )

legend1 <- get_legend(
  p1 +theme(legend.position = 'bottom')
)

legend2 <- get_legend(
  q1+theme(legend.position = 'bottom')
)
legend <-plot_grid(legend1, legend2, rel_widths = c(3,2), scale = 0.9)

grid<-plot_grid(title1, row1, title2, row2, title3, row3, ncol=1, rel_heights = c(.15,1,.15,1,.15,1))

png('plots/downstream_effects.png', width=650, height=500)
plot_grid(grid, legend, ncol=1, rel_heights = c(4,.5))
dev.off()

tbl_2009<-table(spe$CellType[spe$CellCat=='islet'], spe$stage[spe$CellCat=='islet'])
tbl<-table(spe$cell_type[spe$cell_category=='islet'], spe$stage[spe$cell_category=='islet'])
tbl_decompr<-table(spe_6$cell_type[spe$cell_category=='islet'], spe$stage[spe$cell_category=='islet'])

tbl_2009<-scale(tbl_2009, center=F, scale=colSums(tbl_2009))
tbl<-scale(tbl, center=F, scale=colSums(tbl))
tbl_decompr<-scale(tbl_decompr, center=F, scale=colSums(tbl_decompr))

tbl['beta', 'Non-diabetic'] - tbl['beta', 'Onset']
tbl_decompr['beta', 'Non-diabetic'] - tbl_decompr['beta', 'Onset']

t.test(tbl_2009['alpha',], tbl['alpha',], paired=T)
t.test(tbl_2009['beta',], tbl['beta',], paired=T)
t.test(tbl_2009['delta',], tbl['delta',], paired=T)
t.test(tbl_2009['gamma',], tbl['gamma',], paired=T)

t.test(tbl_decompr['alpha',], tbl['alpha',], paired=T)
t.test(tbl_decompr['beta',], tbl['beta',], paired=T)
t.test(tbl_decompr['delta',], tbl['delta',], paired=T)
t.test(tbl_decompr['gamma',], tbl['gamma',], paired=T)

tbl_2009<-table(spe$CellType[spe$CellCat=='immune'], spe$stage[spe$CellCat=='immune'])
tbl<-table(spe$cell_type[spe$cell_category=='immune'], spe$stage[spe$cell_category=='immune'])
tbl_decompr<-table(spe_6$cell_type[spe$cell_category=='immune'], spe$stage[spe$cell_category=='immune'])

(tbl_2009<-t(scale(t(tbl_2009), center=F, scale=rowSums(tbl_2009))))
(tbl<-t(scale(t(tbl), center=F, scale=rowSums(tbl))))
(tbl_decompr<-t(scale(t(tbl_decompr), center=F, scale=rowSums(tbl_decompr))))
