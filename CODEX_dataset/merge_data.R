pacman::p_load(plyr, readr, SpatialExperiment)

targets<-c( "MUC2","SOX9","MUC1","CD31","Synapto","CD49f","CD15","CHGA","CDX2","ITLN1","CD4","CD127","Vimentin","HLADR","CD8","CD11c","CD44","CD16","BCL2","CD3","CD123","CD38","CD90","aSMA","CD21","NKG2D","CD66","CD57","CD206","CD68","CD34","aDef5","CD7","CD36","CD138","CD45RO","Cytokeratin","CD117","CD19","Podoplanin","CD45","CD56","CD69","Ki67","CD49a","CD163","CD161")

B004_CL_1<-read_csv('intensities/B004_CL_reg_1.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B004_CL_1$donor<-'B004'
B004_CL_1$array<-'B004_CL'
B004_CL_2<-read_csv('intensities/B004_CL_reg_2.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B004_CL_2$donor<-'B004'
B004_CL_2$array<-'B004_CL'
B004_CL_3<-read_csv('intensities/B004_CL_reg_3.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B004_CL_3$donor<-'B004'
B004_CL_3$array<-'B004_CL'
B004_CL_4<-read_csv('intensities/B004_CL_reg_4.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B004_CL_4$donor<-'B004'
B004_CL_4$array<-'B004_CL'

B004_SB_1<-read_csv('intensities/B004_SB_reg_1.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B004_SB_1$donor<-'B004'
B004_SB_1$array<-'B004_SB'
B004_SB_2<-read_csv('intensities/B004_SB_reg_2.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B004_SB_2$donor<-'B004'
B004_SB_2$array<-'B004_SB'
B004_SB_3<-read_csv('intensities/B004_SB_reg_3.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B004_SB_3$donor<-'B004'
B004_SB_3$array<-'B004_SB'
B004_SB_4<-read_csv('intensities/B004_SB_reg_4.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B004_SB_4$donor<-'B004'
B004_SB_4$array<-'B004_SB'

B010A_1<-read_csv('intensities/B010A_reg_1.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B010A_1$donor<-'B010'
B010A_1$array<-'B010A'
B010A_2<-read_csv('intensities/B010A_reg_2.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B010A_2$donor<-'B010'
B010A_2$array<-'B010A'
B010A_3<-read_csv('intensities/B010A_reg_3.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B010A_3$donor<-'B010'
B010A_3$array<-'B010A'
B010A_4<-read_csv('intensities/B010A_reg_4.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B010A_4$donor<-'B010'
B010A_4$array<-'B010A'

B010B_1<-read_csv('intensities/B010B_reg_1.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B010B_1$donor<-'B010'
B010B_1$array<-'B010B'
B010B_2<-read_csv('intensities/B010B_reg_2.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B010B_2$donor<-'B010'
B010B_2$array<-'B010B'
B010B_3<-read_csv('intensities/B010B_reg_3.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B010B_3$donor<-'B010'
B010B_3$array<-'B010B'
B010B_4<-read_csv('intensities/B010B_reg_4.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B010B_4$donor<-'B010'
B010B_4$array<-'B010B'

B011A_1<-read_csv('intensities/B011A_reg_1.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B011A_1$donor<-'B011'
B011A_1$array<-'B011A'
B011A_2<-read_csv('intensities/B011A_reg_2.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B011A_2$donor<-'B011'
B011A_2$array<-'B011A'
B011A_3<-read_csv('intensities/B011A_reg_3.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B011A_3$donor<-'B011'
B011A_3$array<-'B011A'
B011A_4<-read_csv('intensities/B011A_reg_4.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B011A_4$donor<-'B011'
B011A_4$array<-'B011A'

B011B_1<-read_csv('intensities/B011B_reg_1.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B011B_1$donor<-'B011'
B011B_1$array<-'B011B'
B011B_2<-read_csv('intensities/B011B_reg_2.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B011B_2$donor<-'B011'
B011B_2$array<-'B011B'
B011B_3<-read_csv('intensities/B011B_reg_3.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B011B_3$donor<-'B011'
B011B_3$array<-'B011B'
B011B_4<-read_csv('intensities/B011B_reg_4.csv', col_select=c(targets, 'Hoechst1', 'DRAQ5', 'Reg', 'Absolute X', 'Absolute Y', 'Cell Size'))
B011B_4$donor<-'B011'
B011B_4$array<-'B011B'



all_cells<-rbind(B004_CL_1, B004_CL_2, B004_CL_3, B004_CL_4,
                 B004_SB_1, B004_SB_2, B004_SB_3, B004_SB_4,
                # B005_CL_1, B005_CL_2, B005_CL_3, B005_CL_4, 
                # B005_SB_2, B005_SB_3, B005_SB_4, # B005_SB_1,
               B010A_1, B010A_2, B010A_3, B010A_4,
               B010B_1, B010B_2, B010B_3, B010B_4,
              B011A_1, B011A_2, B011A_3, B011A_4,
              B011B_1, B011B_2, B011B_3, B011B_4)
                

spe<-SpatialExperiment(
  colData=all_cells[,c('donor', 'array', 'Reg', 'Cell Size')],
  spatialCoords=as.matrix(all_cells[,c('Absolute X','Absolute Y')]),
  rowData=data.frame(channel=1:49, name=colnames(all_cells)[1:49]),
  assays=list(counts=t(all_cells[,1:49]))
)

cells<-as.data.frame(read_csv("CODEX_HuBMAP.csv", lazy=T))
cells$unique_region<-mapvalues(cells$unique_region, from=c('B010_Mid jejunum', 'B011_Mid jejunum', 'B010_Proximal jejunum', 'B011_Proximal jejunum',
'B010_Sigmoid', 'B011_Sigmoid', 'B010_Trans', 'B011_Trans', 'B010_Right', 'B011_Right', 'B010_Left', 'B011_Left'),
to=c('B010_Mid-jejunum', 'B011_Mid-jejunum', 'B010_Proximal Jejunum', 'B011_Proximal Jejunum',
     'B010_Descending - Sigmoid',  'B011_Descending - Sigmoid', 'B010_Transverse','B011_Transverse', 
     'B010_Ascending', 'B011_Ascending', 'B010_Descending', 'B011_Descending' ))

region_tissue_map <- unique(cells[,c('array', 'region', 'Tissue_location')])
tissues <- region_tissue_map$`Tissue_location`
names(tissues) <- paste0(region_tissue_map$array, region_tissue_map$region)

spe$Tissue_location<- as.character(tissues[paste0(spe$array, spe$Reg)] )
spe$unique_region <-paste0(spe$donor, '_', spe$Tissue_location)
spe$tissue <- ifelse(spe$Tissue_location%in%c('Duodenum', 'Ileum', 'Mid-jejunum', 'Proximal Jejunum'), 'SB', 'CL')

colnames(spatialCoords(spe))<-c('x','y')

saveRDS(spe, 'spe_raw.rds')