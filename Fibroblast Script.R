#Fibroblast by subtype UMAP 
#Subset for only fibroblasts
#Note: This was the code used to subset and sub-cluster fibroblasts from the full 
#scRNA Seurat object. We have also included the exact fibrobalst object used in the manuscript figures
#which preserve the cluster assignemnts and UMAP coordinates. 
fibro <- subset(full, subset = celltype == "Fibroblasts")

fibro_cols <- c(
  "APOE Fibroblasts"         = "#4E79A7",  
  "CCL19 Fibroblasts"        = "#E15759", 
  "CLDN Fibroblasts"         = "#59A14F",  
  "COL11A1 Fibroblasts"      = "#F28E2B",  
  "EN1 Fibroblasts"          = "#76B7B2",  
  "FMO2 Fibroblasts"         = "#B07AA1",  
  "Proliferating Fibroblasts"= "#FF9DA7",  
  "SFRP2 Fibroblasts"        = "#9C755F",  
  "TNN Fibroblasts"          = "#EDC948"   
)

p_fibro <- DimPlot(fibro, reduction = "umap",
                   group.by = "subtype", cols = fibro_cols,    
                   label = FALSE, pt.size = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 4))) + ggtitle(NULL)

p_fibro

#Proportion of each fibroblast subtype by sample type 
df <- fibro@meta.data %>%
  mutate(
    disease = factor(
      disease,
      levels = c("Control", "Nonles Lupus", "Les Lupus")
    )
  ) %>%
  count(disease, subtype) %>%
  complete(disease, subtype, fill = list(n = 0)) %>%
  group_by(disease) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup()

p_disease <- ggplot(df, aes(x = disease, y = percent, fill = subtype)) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black",
    linewidth = 0.5
  ) + 
  scale_fill_manual(values = fibro_cols) +
  labs(
    y = "% of Fibroblast Cells",
    x = NULL,
    fill = "Fibroblast subtype"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15),
    legend.position = "right"
  )

p_disease

#Fibroblasts by sample type UMAP 
p_fibro_disease <- DimPlot(fibro, reduction = "umap",
                           group.by = "disease", cols = disease_cols,    
                           label = FALSE, pt.size = 0.5) + ggtitle(NULL) +
  guides(color = guide_legend(override.aes = list(size = 4))) 

p_fibro_disease

#Proportion of sample type within fibroblast subtypes 
df_d <- fibro@meta.data %>%
  count(subtype, disease) %>%
  complete(
    subtype,
    disease,
    fill = list(n = 0)
  ) %>%
  group_by(subtype) %>%
  mutate(percent = n / sum(n) * 100)

p_d <- ggplot(df_d, aes(x = subtype, y = percent, fill = disease)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", size = 0.2) +
  scale_fill_manual(
    values = c(
      "Control"       =  "#4CAF50",
      "Nonles Lupus"  =  "#4A90C2", 
      "Les Lupus"     = "#D62728" 
    ),
    name = ""
  ) +
  labs(y = "% of Each Subset", x = NULL) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_d

#IFNA and IFNG response score UMAP 

IFNA_response_genes <- c(
  "ADAR","B2M","BATF2","BST2","C1S","CASP1","CASP8","CCRL2","CD47","CD74",
  "CMPK2","CNP","CSF1","CXCL10","CXCL11","DDX60","DHX58","EIF2AK2","ELF1",
  "EPSTI1","FAM125A","FAM46A","FTSJD2","GBP2","GBP4","GMPR","HERC6","HLA-C",
  "IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT2","IFIT3","IFITM1",
  "IFITM2","IFITM3","IL15","IL4R","IL7","IRF1","IRF2","IRF7","IRF9","ISG15",
  "ISG20","LAMP3","LAP3","LGALS3BP","LPAR6","LY6E","MOV10","MX1","NCOA7","NMI",
  "NUB1","OAS1","OASL","OGFR","PARP12","PARP14","PARP9","PLSCR1","PNPT1",
  "PRIC285","PROCR","PSMA3","PSMB8","PSMB9","PSME1","PSME2","RIPK2","RNF31",
  "RSAD2","RTP4","SAMD9","SAMD9L","SELL","SLC25A28","SP110","STAT2","TAP1",
  "TDRD7","TMEM140","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TRIM5",
  "TXNIP","UBA7","UBE2L6","USP18","WARS"
)

IFNG_response_genes <- c(
  "ADAR","APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2",
  "BTG1","C1R","C1S","CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5",
  "CCL7","CD274","CD38","CD40","CD69","CD74","CD86","CDKN1A","CFB","CFH",
  "CIITA","CMKLR1","CMPK2","CSF2RB","CXCL10","CXCL11","CXCL9","DDX58","DDX60",
  "DHX58","EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1","FTSJD2",
  "GBP4","GBP6","GCH1","GPR18","GZMA","HERC6","HIF1A","HLA-A","HLA-B","HLA-DMA",
  "HLA-DQA1","HLA-DRB1","HLA-G","ICAM1","IDO1","IFI27","IFI30","IFI35","IFI44",
  "IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFITM2","IFITM3","IFNAR2","IL10RA",
  "IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF2","IRF4",
  "IRF5","IRF7","IRF8","IRF9","ISG15","ISG20","ISOC1","ITGB7","JAK2","KLRK1",
  "LAP3","LATS2","LCP2","LGALS3BP","LY6E","LYSMD2","MARCH1","METTL7B","MT2A","MTHFD2",
  "MVP","MX1","MX2","MYD88","NAMPT","NCOA3","NFKB1","NFKBIA","NLRC5","NMI",
  "NOD1","NUP93","OAS2","OAS3","OASL","OGFR","P2RY14","PARP12","PARP14","PDE4B",
  "PELI1","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1","PRIC285","PSMA2",
  "PSMA3","PSMB10","PSMB2","PSMB8","PSMB9","PSME1","PSME2","PTGS2","PTPN1","PTPN2",
  "PTPN6","RAPGEF6","RBCK1","RIPK1","RIPK2","RNF213","RNF31","RSAD2","RTP4","SAMD9L",
  "SAMHD1","SECTM1","SELP","SERPING1","SLAMF7","SLC25A28","SOCS1","SOCS3","SOD2","SP110",
  "SPPL2A","SRI","SSPN","ST3GAL5","ST8SIA4","STAT1","STAT2","STAT3","STAT4","TAP1",
  "TAPBP","TDRD7","TNFAIP2","TNFAIP3","TNFAIP6","TNFSF10","TOR1B","TRAFD1","TRIM14","TRIM21",
  "TRIM25","TRIM26","TXNIP","UBE2L6","UPP1","USP18","VAMP5","VAMP8","VCAM1","WARS",
  "XAF1","XCL1","ZBP1","ZNFX1"
)

fibro <- AddModuleScore(
  object = fibro,
  features = list(IFNA_response_genes), 
  name = "IFNA_Response_Score"
)

fibro <- AddModuleScore(
  object = fibro,
  features = list(IFNG_response_genes), 
  name = "IFNG_Response_Score"
)

p_IFNA_fibro <- FeaturePlot(
  fibro,
  features  = "IFNA_Response_Score1",  
  reduction = "umap",
  pt.size   = 1,
  raster    = FALSE,
  cols      = c("white", "red"),   
  order     = TRUE,
  min.cutoff = "q05", max.cutoff = "q95"
)

p_IFNA_fibro <- p_IFNA_fibro + ggtitle("IFN-α Response Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + scale_color_gradientn(
    colors = c("grey95", "#fdd0a2", "#fb6a4a", "#cb181d")
  ) +  coord_equal() +
  theme_classic(base_size = 16) +
  labs(x = "umap_1", y = "umap_2") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )


p_IFNG_fibro <- FeaturePlot(
  fibro,
  features  = "IFNG_Response_Score1",   
  reduction = "umap",
  pt.size   = 1,
  raster    = FALSE,
  cols      = c("white", "red"),   
  order     = TRUE,
  min.cutoff = "q05", max.cutoff = "q95"
)

p_IFNG_fibro <- p_IFNG_fibro + ggtitle("IFN-y Response Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + scale_color_gradientn(
    colors = c("grey95", "#fdd0a2", "#fb6a4a", "#cb181d")
  ) +  coord_equal() +
  theme_classic(base_size = 16) +
  labs(x = "umap_1", y = "umap_2") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )

#Violin plots of IFNA and IFNG response scores across fibroblast subtypes and sample types 

plot_df <- fibro@meta.data %>%
  dplyr::select(
    subtype,
    disease,
    IFNA_Response_Score1,
    IFNG_Response_Score1
  ) %>%
  dplyr::filter(!is.na(subtype), !is.na(disease))


fibro_IFNA_violin <- ggplot(plot_df, aes(x = subtype, y = IFNA_Response_Score1, fill = disease)) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = TRUE,
    scale = "width",
    width = 0.9,
    alpha = 0.7
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 1.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#00BA38",
      "Nonles Lupus" = "#619CFF",
      "Les Lupus" = "#F8766D"
    )
  ) +
  theme_classic(base_size = 17) +
  labs(
    x = "Fibroblast subtype",
    y = "IFN-α module score",
    fill = "Disease"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

fibro_IFNA_violin

fibro_IFNG_violin <- ggplot(plot_df, aes(x = subtype, y = IFNG_Response_Score1, fill = disease)) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = TRUE,
    scale = "width",
    width = 0.9,
    alpha = 0.7
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 1.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#00BA38",
      "Nonles Lupus" = "#619CFF",
      "Les Lupus" = "#F8766D"
    )
  ) +
  theme_classic(base_size = 17) +
  labs(
    x = "Fibroblast subtype",
    y = "IFN-y module score",
    fill = "Disease"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

fibro_IFNG_violin

#Fibroblast cellchat chord plots

library(CellChat)
library(circlize)
library(RColorBrewer)

cellchat <- createCellChat(object = full, group.by = "subtype")

CellChatDB <- CellChatDB.human  
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)   

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat_weight <- cellchat@net$weight

cellchat_weight_df <- as.data.frame(cellchat_weight) %>%
  rownames_to_column("source") %>%
  pivot_longer(-source, names_to = "target", values_to = "weight")

#Incoming to CCL19 Fibroblasts
df_chord_in_ccl19 <- cellchat_weight_df %>%
  filter(target == "CCL19 Fibroblasts") %>%
  group_by(source, target) %>%
  summarise(total_weight = sum(weight), .groups = "drop")

top_sources <- df_chord_in_ccl19 %>%
  arrange(desc(total_weight)) %>%
  slice_head(n = 8) %>%
  pull(source)

df_chord_in_ccl19 <- df_chord_in_ccl19 %>%
  filter(source %in% top_sources)

incoming_cells <- unique(df_chord_in_ccl19$source)

pal <- brewer.pal(min(length(incoming_cells), 8), "Set2")
pal <- pal[seq_along(incoming_cells)]
names(pal) <- incoming_cells

grid.col <- c(
  "CCL19 Fibroblasts" = "#228B22",  
  pal
)

tiff("Fig_CCL19_incoming_chord.tiff",
     width = 8,
     height = 8,
     units = "in",
     res = 600,
     compression = "lzw")

circos.clear()
circos.par(start.degree = 90, gap.degree = 4)

chordDiagram(
  x = df_chord_in_ccl19[, c("source", "target", "total_weight")],
  grid.col = grid.col,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow",
  transparency = 0.35
)

dev.off()

#Outgoing from CCL19 Fibroblasts
df_chord_out_ccl19 <- cellchat_weight_df %>%
  filter(source == "CCL19 Fibroblasts") %>%   
  filter(weight > 0) %>%
  group_by(source, target) %>%
  summarise(total_weight = sum(weight), .groups = "drop")

top_targets <- df_chord_out_ccl19 %>%
  arrange(desc(total_weight)) %>%
  slice_head(n = 8) %>%                       
  pull(target)

df_chord_out_ccl19 <- df_chord_out_ccl19 %>%
  filter(target %in% top_targets)

outgoing_targets <- unique(df_chord_out_ccl19$target)

pal <- brewer.pal(min(length(outgoing_targets), 8), "Set2")
pal <- pal[seq_along(outgoing_targets)]
names(pal) <- outgoing_targets

grid.col <- c(
  "CCL19 Fibroblasts" = "#228B22",  
  pal
)

tiff("Fig_CCL19_outgoing_chord.tiff",
     width = 10,
     height = 10,
     units = "in",
     res = 600,
     compression = "lzw")

circos.clear()
circos.par(start.degree = 90, gap.degree = 4)

chordDiagram(
  x = df_chord_out_ccl19[, c("source", "target", "total_weight")],
  grid.col = grid.col,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow",
  transparency = 0.35
)

dev.off()


#Fibroblast and myeloids cellchat bubble plots 
library(CellChat)

cellchat <- createCellChat(object = full, group.by = "subtype")

CellChatDB <- CellChatDB.human  
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)   

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

myeloid <- c("pDCs", "IFN_high_pDC", "IFN_low_pDC", "pDCs_2","cDC1","cDC2A","cDC2B","Macrophages")

fibros <- c(
  "CCL19 Fibroblasts",
  "SFRP2 Fibroblasts",
  "FBN1 Fibroblasts",
  "Proliferating Fibroblasts",
  "FMO2 Fibroblasts",
  "COL11A1 Fibroblasts",
  "CLDN Fibroblasts",
  "APOE Fibroblasts"
)

ccl19 <- "CCL19 Fibroblasts"

p <- netVisual_bubble(
  cellchat,
  sources.use = myeloid,
  targets.use = "CCL19 Fibroblasts",
  pairLR.use = pairLR_df,
  remove.isolate = TRUE   
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p2 <- netVisual_bubble(
  cellchat,
  targets.use = myeloid,
  sources.use = "CCL19 Fibroblasts",
  pairLR.use = pairLR_df,
  remove.isolate = TRUE   
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
