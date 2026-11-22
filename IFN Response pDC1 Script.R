library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

#pDC subclustering 
full <- readRDS("fullseurat.rds") #read in full scRNA Seurat object here

Idents(full) <- full$celltype
pDC1 <- subset(full,idents=c("pDCs"))

pDC1<- NormalizeData(pDC1)
pDC1 <- FindVariableFeatures(pDC1, selection.method = "vst")
all.genes <- c(rownames(pDC1))

pDC1 <- ScaleData(pDC1, verbose = FALSE,features=all.genes)
pDC1 <- RunPCA(pDC1, npcs = 30, verbose = FALSE)

data.use <- Stdev(object = pDC1, reduction = 'pca')

index <- c()
for (i in 1:(length(data.use)-1)){
  d1 <- data.use[i]
  d2 <- data.use[i+1]
  if (d2 > 0.995*d1) {
    index <- append(index, i+1)
  }
}

pc.dim <- ifelse(length(index) == 0, length(data.use), min(index))

# t-SNE and Clustering / Harmony Embeddings
options(repr.plot.height = 2.5, repr.plot.width = 6)
pDC1 <- pDC1 %>% RunHarmony("batch", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(pDC1, 'harmony')

# t-SNE and Clustering
pDC1 <- FindNeighbors(pDC1, reduction = "harmony",dims = 1:pc.dim)
pDC1 <- FindClusters(pDC1, resolution = 1.0)
pDC1 <- RunUMAP(pDC1, reduction = "harmony",dims = 1:pc.dim)

pDC1 <- FindNeighbors(pDC1, dims = 1:2, reduction="umap")
pDC1 <- FindClusters(pDC1, resolution = 1)

#pDC1 subcluster UMAP 
pDC_subcluster_plot <- DimPlot(pDC1, reduction = "umap",
                               group.by = "seurat_clusters",      # column in meta.data
                               label = FALSE, pt.size = 1.2) +
  guides(color = guide_legend(override.aes = list(size = 4))) + ggtitle(NULL) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )

pDC_subcluster_plot

#pDC marker genes feature plots
pDC1_feats <- c("CLEC4C","LILRA4","TCF4","IRF7","GZMB")
pDC_markers_umap <- FeaturePlot(pDC1, pDC1_feats, pt.size = 0.25) 

pDC_markers_umap

#IFNA and IFNG response module scoring 
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

pDC1 <- AddModuleScore(
  object = pDC1,
  features = list(IFNA_response_genes), 
  name = "IFNA_Response_Score"
)

pDC1 <- AddModuleScore(
  object = pDC1,
  features = list(IFNG_response_genes), 
  name = "IFNG_Response_Score"
)

#pDC interferon module score UMAP
p_IFNA <- FeaturePlot(
  pDC1,
  features  = "IFNA_Response_Score1",   
  reduction = "umap",
  pt.size   = 1.5,
  raster    = FALSE,
  cols      = c("white", "red"),   
  order     = TRUE,
  min.cutoff = "q05", max.cutoff = "q95"
)

p_IFNA <- p_IFNA + ggtitle("IFN-α Response Score") +
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

p_IFNA

p_IFNG <- FeaturePlot(
  pDC1,
  features  = "IFNG_Response_Score1",   
  reduction = "umap",
  pt.size   = 1.5,
  raster    = FALSE,
  cols      = c("white", "red"),   
  order     = TRUE,
  min.cutoff = "q05", max.cutoff = "q95"
)

p_IFNG <- p_IFNG + ggtitle("IFN-y Response Score") +
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

p_IFNG

#IFN-high vs IFN-low pDC annotated UMAP

pDC1$IFN_group2 <- dplyr::case_when(
  pDC1$seurat_clusters %in% c(2, 3, 5, 6, 7, 9, 10, 11, 12) ~ "IFN_high_pDC1",
  pDC1$seurat_clusters %in% c(0, 1, 4, 8, 13) ~ "IFN_low_pDC1",
  TRUE ~ NA_character_
)

IFN_group_plot <- DimPlot(pDC1, reduction = "umap",
                          group.by = "IFN_group2",     
                          label = FALSE, pt.size = 1.2) +
  guides(color = guide_legend(override.aes = list(size = 4))) + ggtitle(NULL) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )

IFN_group_plot

#IFN-high vs IFN-low volcano plots 
Idents(pDC1) <- "IFN_group2"   

markers_IFN <- FindMarkers(
  pDC1,
  ident.1 = "IFN_high_pDC",
  ident.2 = "IFN_low_pDC",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

label_genes_IFN <- c(
  "TLR7", "TLR9", "UBE2L3", "NLRP2", "LYZ", "IFI30", 
  "B2M", "CD74", "PLXNA4", "MT-CO2", "MT-CO3", "MT-ND2", "MT-ND4"
)

plot_v_IFN <- ggplot(volc_df_IFN, aes(x = avg_log2FC, y = negLog10P)) +
  geom_point(aes(color = signif), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "IFN_high_pDC1" = "red",
    "IFN_low_pDC1" = "blue",
    "NS" = "grey70"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(volc_df_IFN, Gene %in% label_genes_IFN),
    aes(label = Gene),
    size = 5,
    box.padding = 0.4,
    max.overlaps = Inf
  ) +
  labs(
    x = "log2 Fold Change (IFN-high pDC1 vs IFN-low pDC1)",
    y = "-log10(p-value)",
    color = "Enrichment"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15),
    legend.position = "right"
  )

#IFN-high pDC1 enriched pathways 

Idents(pDC1) <- "IFN_group2"   

markers_IFN <- FindMarkers(
  pDC1,
  ident.1 = "IFN_high_pDC",
  ident.2 = "IFN_low_pDC",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

pDC_IFN_high_genes <- rownames(markers_IFN)[markers_IFN$avg_log2FC > 1 & markers_IFN$p_val_adj < 0.05]
pDC_IFN_low_genes <- rownames(markers_IFN)[markers_IFN$avg_log2FC < -1 & markers_IFN$p_val_adj < 0.05]

ego_pDC_IFN_high <- enrichGO(
  gene          = pDC_IFN_high_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

GO_pDC_IFN_high <- dotplot(ego_pDC_IFN_high, showCategory = 15) +
  ggtitle("IFN-high pDC1 enriched pathways")

#IFN high/low pDCs pseudotime (monocle3)

library(monocle3)
data <- GetAssayData(pDC1, slot = "counts", assay = "RNA")

cell_metadata <- pDC1@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(data),
  row.names = rownames(data)
)

pDC1 <- new_cell_data_set(
  data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

pDC1 <- preprocess_cds(pDC1, num_dim = 50)
pDC1 <- reduce_dimension(pDC1, reduction_method = "UMAP")
pDC1 <- cluster_cells(pDC1)
pDC1 <- learn_graph(pDC1)

#Define Root
pDC1 <- order_cells(pDC1, root_cells = colnames(pDC1)[pDC1@colData$subtype == "Quiescent"])
plot_cells(pDC1, color_cells_by = "subtype")
plot_cells(pDC1, color_cells_by = "pseudotime")

plot_cells(
  pDC1,
  color_cells_by = "IFN_group2",
  cell_size = 2.5,
  label_cell_groups = FALSE,      
  label_groups_by_cluster = FALSE,
  group_label_size = 6,         
  graph_label_size = 5,           
) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(size = 20),
    axis.text  = ggplot2::element_text(size = 16)
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))


#IPA pathway figures
library(readxl)

ipa <- read_excel("IFN high vs low pdcs pathways.xls")

colnames(ipa) <- c(
  "Pathway",
  "negLogP",
  "Ratio",
  "zscore",
  "Molecules"
)

ipa <- ipa %>%
  mutate(
    negLogP = as.numeric(negLogP),
    zscore  = as.numeric(zscore)
  )

ipa_top <- ipa %>%
  arrange(desc(negLogP)) %>%
  slice_head(n = 35) %>%
  mutate(
    Direction = case_when(
      zscore > 0  ~ "Activated",
      zscore < 0  ~ "Inhibited",
      TRUE        ~ "No prediction"
    )
  ) %>%
  filter(Direction != "No prediction")

ipa_pathway <- ggplot(
  ipa_top,
  aes(
    x = negLogP,
    y = reorder(Pathway, negLogP),
    fill = Direction
  )
) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      "Activated" = "#D73027",
      "Inhibited" = "#4575B4",
      "No prediction" = "grey70"
    )
  ) +
  labs(
    x = expression(-log[10](italic(p))),
    y = NULL,
    fill = "Predicted activity"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

#IPA top predicted upstream regulators 
ipa_up <- read_excel("IFN high vs low pdcs URS.xls")

library(dplyr)
library(ggplot2)

up_df <- ipa_up %>% 
  filter(
    `Molecule Type` %in% c("cytokine", "growth factor"),
    (`Flags` == "" | is.na(`Flags`)),              
    `Predicted Activation State` == "Activated",
    as.numeric(`p-value of overlap`) < 0.05           
  ) %>%  
  transmute(
    Regulator = `Upstream Regulator`,
    negLog10P = -log10(`p-value of overlap`),
    zscore = `Activation z-score`,
    Activity = case_when(
      zscore > 0 ~ "Activated",
      zscore < 0 ~ "Inhibited",
      TRUE ~ "No prediction"
    )
  ) %>%
  filter(!is.na(zscore)) %>%         
  arrange(desc(negLog10P)) %>%
  slice_head(n = 10)                 


IPA_up_plot <- ggplot(up_df, aes(
  x = negLog10P,
  y = reorder(Regulator, negLog10P),
  fill = Activity
)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "Activated" = "#D73027",
      "Inhibited" = "grey60"
    )
  ) +
  labs(
    x = expression(log[2]~"(P value)"),
    y = NULL,
    title = "Upstream regulators of DE genes"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 17)
  )

# Cellchat Chord Diagram of Incoming Signal Weights to IFN-high/low pDCs 
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

cellchat_weight_df_incoming_to_pDC <- cellchat_weight_df %>%
  filter(target %in% c("IFN_high_pDC", "IFN_low_pDC")) %>%
  group_by(target, source) %>%
  summarise(total_weight = sum(weight), .groups = "drop")

cellchat_weight_df_incoming_to_pDC_chord <- cellchat_weight_df_incoming_to_pDC %>%
  filter(total_weight > 0) %>%
  filter(target %in% c("IFN_high_pDC", "IFN_low_pDC"))

top_sources <- cellchat_weight_df_incoming_to_pDC_chord %>%
  group_by(source) %>%
  summarise(total = sum(total_weight)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%   
  pull(source)

cellchat_weight_df_incoming_to_pDC_chord <- cellchat_weight_df_incoming_to_pDC_chord %>%
  filter(source %in% top_sources)

#Assign colors
incoming_cells <- setdiff(unique(cellchat_weight_df_incoming_to_pDC_chord$source),
                          c("IFN_high_pDC", "IFN_low_pDC"))

pal <- brewer.pal(length(incoming_cells), "Set2")
names(pal) <- incoming_cells

grid.col <- c(
  "IFN_high_pDC" = "#F8766D",
  "IFN_low_pDC"  = "#01BFC4",
  pal
)

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 4
)

chordDiagram(
  x = df_chord_in[, c("source", "target", "total_weight")],
  grid.col = grid.col,
  transparency = 0.35,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow"
)
