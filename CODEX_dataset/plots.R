pacman::p_load(SpatialExperiment, cowplot, tidyverse, caret, cvms, ggimage, ggrepel,
               corrplot, randomcoloR, dittoSeq, scater, viridis, ggpubr,
               Rphenograph ,igraph, rstatix, scran, bluster)

spe<-readRDS('spe.rds')

dir.create('plots')

# Co-expression analysis
co_exp_matrix<-read.csv('data/co_exp_matrix.csv', row.names=1, check.names=F)

functional<-c( 'CD45', "CD45RO", 'CD38', 'CD138', 'CD161', 'NKG2D' , 'CD16', 'CD163', 'CD206', 'CD68',
               'CD11c', "CD7"  ,"CD3" , "CD8",  "CD4"  ,  "CD19"  , "CD21" ,  "CD123" , "CD127",  "CD69"  , 
               'Cytokeratin', 'SOX9', 'aDef5', "MUC1" , "MUC2" , "CHGA" , "CDX2"  ,"ITLN1", 
               "aSMA", "Podoplanin" , "CD90")
co_exp_matrix<-co_exp_matrix[functional, functional]

pdf('plots/co_oc_mat.pdf')
corrplot(as.matrix(co_exp_matrix), method="shade", addCoef.col = 1, tl.cex = 0.8, number.cex = 0.5, type="upper",
         is.corr=F, col=colorRampPalette(c("blue","white","red"))(100))
dev.off()

co_exp_matrix_imc<-read.csv('../IMC_dataset/data/co_exp_matrix.csv', row.names=1, check.names=F)

islet<-c("SYP", "CD99", "PTPRN", "SLC2A1", "GCG", "PCSK2",  "INS", "PIN", "NKX6-1", "IAPP",  "SST", "PPY")
exocrine<-c("AMY2A", "KRT19", "CD44")
immune<-c("CD45", "CD45RA", "CD68", "MPO", "CD20", "CD3e", "CD4", "CD8a", "FOXP3", "CD38")
stromal<-c("CD31", "SMA")

inter_cat_pairs<-rbind(expand.grid(immune, islet), expand.grid(exocrine, islet),expand.grid(stromal, islet), expand.grid(exocrine, immune), expand.grid(stromal, exocrine), expand.grid(stromal, immune))

imc_rates<-apply(inter_cat_pairs, 1, function(p) co_exp_matrix_imc[p[1],p[2]])

immune<-c( 'CD45', "CD45RO", 'CD38', 'CD138', 'CD161', 'NKG2D' , 'CD16', 'CD163', 'CD206', 'CD68',
           'CD11c', "CD7"  ,"CD3" , "CD8",  "CD4"  ,  "CD19"  , "CD21" ,  "CD123" , "CD127",  "CD69" )
epithelial<-c('Cytokeratin', 'SOX9', 'aDef5', "MUC1" , "MUC2" , "CHGA" , "CDX2"  ,"ITLN1")
stromal<-c("aSMA", "Podoplanin" , "CD90")

inter_cat_pairs<-rbind(expand.grid(immune, epithelial), expand.grid(epithelial, stromal),expand.grid(stromal, immune))

codex_rates<-apply(inter_cat_pairs, 1, function(p) co_exp_matrix[p[1],p[2]])

t.test(imc_rates, codex_rates)

# Cell typing

celltype_mean <- aggregateAcrossCells(as(spe, "SingleCellExperiment"),  
                                      ids = spe$cell_type, 
                                      statistics = "mean",
                                      use.assay.type = "counts")

cell_types<-c('Enterocyte', 'Goblet', 'Paneth', 'Neuroendocrine', 'TA', 
  'B', 'CD4+ T', 'CD8+ T', 'NK', 'Plasma', 'DC', 'Neutrophil', 'M1 Macrophage', 'M2 Macrophage',
  'Smooth muscle', 'Endothelial', 'Lymphatic', 'Nerve', 'ICC', 'Stroma')

celltype_mean$cell_type <- factor(celltype_mean$cell_type, levels=cell_types)

type_markers <- c("Cytokeratin", "MUC2" ,      "aDef5"   ,   "CHGA"    ,   "SOX9"     ,  "CD45",  "CD19"     ,  "CD21"  ,     "CD4"   ,     "CD3" ,      
                   "CD8" ,       "CD161" ,     "NKG2D"  ,    "CD38"   ,    "CD138"   ,   "HLADR"   ,   "CD11c" ,     "CD15" ,     
                   "CD66" ,      "CD68"  ,     "CD163"  ,    "CD206"   ,   "aSMA"   ,    "CD31" ,      "Podoplanin" ,"Synapto" ,   "CD56",      
                   "CD117"  )


pdf('plots/cell_types.pdf')
dittoHeatmap(celltype_mean,
             genes=type_markers,
             assay = "counts", 
             cluster_cols = F, 
             cluster_rows= F,
             scaled.to.max = TRUE,
             heatmap.colors.max.scaled = plasma(100),
             annot.by = c('cell_type', 'cell_category'),
             annotation_colors = list(cell_category= setNames(c("pink", "gray50", "gold"), c('Epithelial', 'Immune','Stromal'))
                                      ))
dev.off()

set.seed(976)
spe_mini<-spe[,sample(seq_len(ncol(spe)), 100000)]

set.seed(38)
spe_mini<-runUMAP(spe_mini, exprs_values = "exprs", subset_row = type_markers)
colnames(spe_mini) <- paste0('cell', seq_len(ncol(spe_mini)))
saveRDS(spe_mini, 'spe_mini.rds')

# colors<-sample(distinctColorPalette(20))
colors<-c("forestgreen",  "#6875D9",  "lightsalmon",  "grey", "tomato", "#E158D3" , "royalblue",
  "#848877", "#FDBF6F", "sienna", "#CAB2D6",  "darkblue", "skyblue", 
  "#74457B", "#BAE6B0",  "pink", "gold", "#A968E1", "turquoise", "#E9D588")
png('plots/UMAP_types.png', height=1000, width = 2500)
p<-dittoDimPlot(spe_mini, 
                var = "cell_type", 
                reduction.use = "UMAP", 
                size = 0.2,
                labels.size =10, 
                legend.size = 20,
                do.label = F,
                legend.show=TRUE,
                color.panel = colors) +
  ggtitle("Neural Network + kNN")+
  theme(plot.title = element_text(size = 40), legend.title = element_text('Cell Type'), legend.text = element_text(size = 30), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30))+
  guides(color = guide_legend(override.aes = list(size=10), ncol=2) )
q<-dittoDimPlot(spe_mini, 
                var = "stellar_type", 
                reduction.use = "UMAP", 
                size = 0.2,
                labels.size =10, 
                do.label = F,
                legend.show=FALSE,
                color.panel = colors) +
  ggtitle("STELLAR")+
  theme(plot.title = element_text(size = 40), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30))
plot_grid(q,p, rel_widths = c(3, 5), scale = 0.9)
dev.off()

# Compression

spe_8<-readRDS('spe_decompr.rds')

sum(spe$cell_type==spe_8$cell_type)/ncol(spe)
sum(spe$cell_type[spe$cell_category=='Epithelial']==spe_8$cell_type[spe$cell_category=='Epithelial'])/sum(spe$cell_category=='Epithelial')
sum(spe$cell_type[spe$cell_category=='Immune']==spe_8$cell_type[spe$cell_category=='Immune'])/sum(spe$cell_category=='Immune')
sum(spe$cell_type[spe$cell_category=='Stromal']==spe_8$cell_type[spe$cell_category=='Stromal'])/sum(spe$cell_category=='Stromal')

conf_mat<-confusionMatrix(as.factor(spe$cell_type), as.factor(spe_8$cell_type), mode='prec_recall')
cat_type_map<-list(
  'Epithelial'=c( 'Neuroendocrine', 'Goblet', 'Paneth', 'TA'), 
  'Stromal'=c('Smooth muscle', 'Endothelial', 'Nerve', 'Lymphatic', 'ICC'),
  'Immune'=c('CD8+ T', 'CD4+ T', 'B', 'Plasma', 'M1 Macrophage', 'M2 Macrophage', 'Neutrophil', 'DC', 'NK')
)

tbl_immune<-as.tibble(conf_mat$table[cat_type_map[['Immune']], cat_type_map[['Immune']]])
a<-plot_confusion_matrix(tbl_immune, 
                         target_col = "Reference", 
                         prediction_col = "Prediction",
                         class_order = c( 'M2 Macrophage', 'M1 Macrophage', 'Neutrophil', 'Plasma', 'DC',  'NK', 'CD8+ T', 'CD4+ T',  'B'),
                         rm_zero_text =F, rm_zero_percentages = F, add_zero_shading = F,
                         add_arrows = F, 
                         # font_counts = font(size=5), font_row_percentages = font(size=3), font_col_percentages = font(size=3),
                         counts_col = "n", add_normalized = FALSE,  rotate_y_text = F)+ 
  theme(axis.title.x.top=element_text(size=20), axis.title.y=element_text(size=20), axis.text.x = element_text(angle = 45, size = 15, hjust=-.05), axis.text.y = element_text(size = 15))+
  labs(x='Uncompressed', y='Compressed')


tbl_epithelial<-as.tibble(conf_mat$table[c(cat_type_map[['Epithelial']],'Enterocyte'), c(cat_type_map[['Epithelial']],'Enterocyte')])
b<-plot_confusion_matrix(tbl_epithelial, 
                         target_col = "Reference", 
                         prediction_col = "Prediction",
                         class_order = c('Neuroendocrine', 'Paneth', 'Goblet', 'TA', 'Enterocyte'),
                         rm_zero_text =F, rm_zero_percentages = F, add_zero_shading = F,
                         add_arrows = F, 
                         # font_counts = font(size=5), font_row_percentages = font(size=3), font_col_percentages = font(size=3),
                         counts_col = "n", add_normalized = FALSE, rotate_y_text = F)+ 
  theme(axis.title.x.top=element_text(size=20), axis.title.y=element_text(size=20), axis.text.x = element_text(angle = 45, size = 15, hjust=-.05), axis.text.y = element_text(size = 15))+
  labs(x='Uncompressed', y='Compressed')

tbl_stromal<-as.tibble(conf_mat$table[c(cat_type_map[['Stromal']],'Stroma'), c(cat_type_map[['Stromal']],'Stroma')])
c<-plot_confusion_matrix(tbl_stromal, 
                         target_col = "Reference", 
                         prediction_col = "Prediction",
                         class_order = c('Stroma', 'Smooth muscle', 'Endothelial', 'Lymphatic', 'Nerve', 'ICC'),
                         rm_zero_text =F, rm_zero_percentages = F, add_zero_shading = F,
                         add_arrows = F, 
                         # font_counts = font(size=5), font_row_percentages = font(size=3), font_col_percentages = font(size=3),
                         counts_col = "n", add_normalized = FALSE, rotate_y_text = F)+ 
  theme(axis.title.x.top=element_text(size=20), axis.title.y=element_text(size=20), axis.text.x = element_text(angle = 45, size = 15, hjust=-.05), axis.text.y = element_text(size = 15))+
  labs(x='Uncompressed', y='Compressed')

png('plots/confusion_matrix.png', width=1800, height = 600)
plot_grid(b, a, c, align='hv', nrow=1, scale=.9)
dev.off()

comp<-read.csv('data/reconstruction_pairs_long.csv', row.names = 1)
comp_type<-read.csv('data/typing_accuracy.csv', row.names = 1)

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

png('plots/boxplots_codex.png', width = 750, height = 250)
plot_grid(g1, g2, g3, nrow=1, labels='AUTO', scale=0.95)
dev.off()

# Downstream effects

tbl_decompr<-table(spe_8$cell_type[spe$cell_category=='Epithelial'], spe$unique_region[spe$cell_category=='Epithelial'])
tbl_decompr<-t(scale(tbl_decompr, center=F, scale=colSums(tbl_decompr)))
tbl_decompr<-cbind(rownames(tbl_decompr), tbl_decompr)
colnames(tbl_decompr)[1]<-'unique_region'
images<-unique(colData(spe)[,c('unique_region', 'tissue')])
df<-inner_join(as.data.frame(tbl_decompr), as.data.frame(images), by = 'unique_region')

df$Goblet<-as.numeric(df$Goblet)
df$Paneth<-as.numeric(df$Paneth)

t_test(df,Goblet~tissue, alternative = 'greater')
t_test(df,Paneth~tissue, alternative = 'less')

tbl_decompr<-table(spe_8$cell_type[spe$cell_category=='Immune'], spe$unique_region[spe$cell_category=='Immune'])
tbl_decompr<-t(scale(tbl_decompr, center=F, scale=colSums(tbl_decompr)))
tbl_decompr<-cbind(rownames(tbl_decompr), tbl_decompr)
colnames(tbl_decompr)[1]<-'unique_region'
df<-inner_join(as.data.frame(tbl_decompr), as.data.frame(images), by = 'unique_region')

df$DC<-as.numeric(df$DC)
df$`CD8+ T`<-as.numeric(df$`CD8+ T`)

t_test(df, DC~tissue, alternative = 'greater')
t_test(df,`CD8+ T`~tissue, alternative = 'less')

spe$Tissue_location<-factor(spe$Tissue_location, levels=c("Duodenum" ,  "Proximal Jejunum",    "Mid-jejunum" ,       "Ileum",                      
                                                          "Ascending",            "Transverse",  "Descending" ,         "Descending - Sigmoid"   ))
spe_8$Tissue_location<-factor(spe_8$Tissue_location, levels=c("Duodenum" ,  "Proximal Jejunum",    "Mid-jejunum" ,       "Ileum",                      
                                                              "Ascending",            "Transverse",  "Descending" ,         "Descending - Sigmoid"   ))
cell_types<-c('Enterocyte', 'Goblet', 'Paneth', 'Neuroendocrine', 'TA', 
              'B', 'CD4+ T', 'CD8+ T',  'NK', 'Plasma', 'DC', 'Neutrophil', 'M1 Macrophage', 'M2 Macrophage',
              'Smooth muscle', 'Endothelial', 'Lymphatic', 'Nerve', 'ICC', 'Stroma')

spe$cell_type<-factor(spe$cell_type, levels=cell_types)
spe$stellar_type<-factor(spe$stellar_type, levels=cell_types)
spe_8$cell_type<-factor(spe_8$cell_type, levels=cell_types)

a<-ggplot(as.data.frame(colData(spe)[spe$cell_category=='Stromal',]), aes(x = Tissue_location, fill = cell_type)) +
  geom_bar(position='fill')+
  labs(
    x = "Segment",
    y = "Stromal cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(name='Cell Type',palette = "Paired")+
  theme(legend.position = "none",
        text = element_text(size = 20))
b<-ggplot(as.data.frame(colData(spe)[spe$cell_category=='Epithelial',]), aes(x = Tissue_location, fill = cell_type)) +
  geom_bar(position='fill')+
  labs(
    x = "Segment",
    y = "Epithelial cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(name='Cell Type', palette = "Set2")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        text = element_text(size = 20))
c<-ggplot(as.data.frame(colData(spe)[spe$cell_category=='Immune',]), aes(x = Tissue_location, fill = cell_type)) +
  geom_bar(position='fill')+
  labs(
    x = "Segment",
    y = "Immune cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(name='Cell Type', palette = "Set1")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        text = element_text(size = 20))

d<-ggplot(as.data.frame(colData(spe)[spe$stellar_category=='Stromal',]), aes(x = Tissue_location, fill = stellar_type)) +
  geom_bar(position='fill')+
  labs(
    x = "Segment",
    y = "Stromal cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(palette = "Paired")+
  theme(legend.position = "none",
        text = element_text(size = 20))
e<-ggplot(as.data.frame(colData(spe)[spe$stellar_category=='Epithelial',]), aes(x = Tissue_location, fill = stellar_type)) +
  geom_bar(position='fill')+
  labs(
    x = "Segment",
    y = "Epithelial cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        text = element_text(size = 20))
f<-ggplot(as.data.frame(colData(spe)[spe$stellar_category=='Immune',]), aes(x = Tissue_location, fill = stellar_type)) +
  geom_bar(position='fill')+
  labs(
    x = "Segment",
    y = "Immune cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        text = element_text(size = 20))

g<-ggplot(as.data.frame(colData(spe_8)[spe$cell_category=='Stromal',]), aes(x = Tissue_location, fill = cell_type)) +
  geom_bar(position='fill')+
  labs(
    x = "Segment",
    y = "Stromal cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(palette = "Paired")+
  theme(legend.position = "none",
        text = element_text(size = 20))
h<-ggplot(as.data.frame(colData(spe_8)[spe$cell_category=='Epithelial',]), aes(x = Tissue_location, fill = cell_type)) +
  geom_bar(position='fill')+
  labs(
    x = "Segment",
    y = "Epithelial cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        text = element_text(size = 20))
i<-ggplot(as.data.frame(colData(spe_8)[spe$cell_category=='Immune',]), aes(x = Tissue_location, fill = cell_type)) +
  geom_bar(position='fill')+
  labs(legend.position = "none",
       x = "Segment",
       y = "Immune cells (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        text = element_text(size = 20))

title1 <- ggdraw() + draw_label("STELLAR", fontface = 'bold', hjust=0.5, size=25) 
title2 <- ggdraw() + draw_label("Uncompressed", fontface = 'bold', hjust=0.5, size=25)
title3 <- ggdraw() + draw_label("Compressed", fontface = 'bold', hjust=0.5, size=25)

titles <- plot_grid(title1, title2, title3, NULL, nrow=1, rel_widths = c(1,1,1,.75))

legend1 <- get_legend( b +theme(legend.position = 'right'))
legend2 <- get_legend( c +theme(legend.position = 'right'))
legend3 <- get_legend( a +theme(legend.position = 'right'))

legends <- plot_grid(legend1, legend2, legend3, NULL, ncol=1, rel_heights = c(1,1,1,.7))

grid <- plot_grid(e,b,h, f,c,i, d,a,g, rel_heights = c(1,1,1.7), nrow=3)
grid_legend<-plot_grid(grid, legends,  rel_widths = c(3,.75))

png('plots/downstream_effects.png' ,width=1000, height=1000)
plot_grid(titles, grid_legend, ncol=1, rel_heights = c(.5,3.5))
dev.off()

