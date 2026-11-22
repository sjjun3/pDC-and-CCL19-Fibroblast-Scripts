library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

#Celltype UMAP
full <- readRDS("fullseurat.rds") #read in full scRNA Seurat object here

DimPlot(
  full,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,      
  repel = TRUE       
) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

#Myeloids by Subtype UMAP
#Subset for only pDCs and myeloids
mye <- subset(
  full,
  subset = celltype %in% c("pDCs", "Myeloid Cells")
)

mye$subtype[mye$subtype == "pDCs"] <- "pDC1"
mye$subtype[mye$subtype == "pDCs_2"] <- "pDC2"

mye <- NormalizeData(mye)
mye <- FindVariableFeatures(mye, selection.method = "vst")
all.genes <- c(rownames(mye))

#Scaling 
mye <- ScaleData(mye, verbose = FALSE,features=all.genes)
mye <- RunPCA(mye, npcs = 30, verbose = FALSE)

data.use <- Stdev(object = mye, reduction = 'pca')

index <- c()
for (i in 1:(length(data.use)-1)){
  d1 <- data.use[i]
  d2 <- data.use[i+1]
  if (d2 > 0.995*d1) {
    index <- append(index, i+1)
  }
}

pc.dim <- ifelse(length(index) == 0, length(data.use), min(index))
options(repr.plot.height = 2.5, repr.plot.width = 6)
mye <- mye %>% RunHarmony("sample", plot_convergence = TRUE)

#Harmony Embeddings
harmony_embeddings <- Embeddings(mye, 'harmony')

# t-SNE and Clustering
mye <- FindNeighbors(mye, reduction = "harmony",dims = 1:pc.dim)
mye <- FindClusters(mye, resolution = 1.0)
mye <- RunUMAP(mye, reduction = "harmony",dims = 1:pc.dim)

mye <- FindNeighbors(mye, dims = 1:2, reduction="umap")
mye <- FindClusters(mye, resolution = 1)

myeloid_cols <- c(
  "Macrophages"         = "#8E44AD",  
  "Langerhans Cells"    = "#1F77B4",  
  "cDC2B"               = "#2CA02C", 
  "cDC1"                = "#E74C3C",  
  "pDC1"                = "#17BECF", 
  "cDC2A"               = "#FF7F0E",  
  "pDC2"                = "#BCBD22",  
  "Proliferating Cells" = "#FF66CC"   
)

new_order <- names(myeloid_cols)

p1 <- DimPlot(
  mye,
  reduction = "umap",
  group.by = "subtype",
  cols = myeloid_cols,   
  label = FALSE,
  pt.size = 0.5
) + 
  ggtitle(NULL)


#Myeloid subtype percentages across each disease category 
df <- mye@meta.data %>%
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
  scale_fill_manual(values = myeloid_cols, breaks = new_order) +
  labs(
    y = "% of Myeloid Cells",
    x = NULL,
    fill = "Myeloid subtype"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15),
    legend.position = "right"
  )

#Myeloid markers DotPlot 
Idents(mye) <- "subtype"

markers <- FindAllMarkers(
  mye,
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

markers_filt <- markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.5
  )

write.csv(markers_filt, "myeloid_markers.csv", row.names = FALSE)

#Top 7 marker genes with highest log2FC for each subtype
genes_use <- c(
  # Macrophages
  "LILRB5","PLEKHG5","MAMDC2","RNASE1","IGF1","PLTP","FOLR2",
  
  # cDC2A
  "LAMP3","BIRC3","GPR157","CCR7","CD83","TRAF1","NR4A3",
  
  # cDC1
  "CLNK","XCR1","DNASE1L3","CLEC9A","CADM1","SNX22","ENPP1",
  
  # Proliferating Cells
  "MKI67","NCAPG","NUSAP1","CENPM","KIFC1","IQGAP3",
  
  # cDC2B
  "CD1C","CLEC10A","FCER1A","IL1R2","ITGAX",
  
  # Langerhans Cells
  "CD207","FCGBP","CD1A","SLC18A2","TACSTD2","LIMS2","RASAL1",
  
  # pDCs_2
  "SPEG","PPP1R14A","SIGLEC6","KCNK17","CX3CR1","LTK","AXL",
  
  # pDCs
  "PTPRS","CLEC4C","IL3RA","BCL11A","LILRA4","GZMB","IRF7","TCF4"
)

#Ordering top genes within each subtype
avg_exp <- AverageExpression(
  mye,
  assays = "RNA",
  slot = "data"
)$RNA

avg_exp_sub <- avg_exp[genes_use, , drop = FALSE]

gene_order_df <- as.data.frame(avg_exp_sub) |>
  tibble::rownames_to_column("gene") |>
  dplyr::rowwise() |>
  dplyr::mutate(
    dominant_subtype = colnames(avg_exp_sub)[
      which.max(c_across(where(is.numeric)))
    ],
    max_expr = max(c_across(where(is.numeric)))
  ) |>
  dplyr::ungroup()

subtype_order <- c(
  "cDC1",
  "cDC2A",
  "cDC2B",
  "Langerhans Cells",
  "Macrophages",
  "pDCs",
  "pDCs_2",
  "Proliferating Cells"
)

plot_markers <- DotPlot(
  mye,
  features = gene_order,
  group.by = "subtype"
) +
  coord_flip() +
  scale_x_discrete(limits = rev(gene_order)) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#Myeloids by disease status UMAP 
disease_cols <- c(
  "Control"       =  "#4CAF50",
  "Nonles Lupus"  =  "#4A90C2", 
  "Les Lupus"     = "#D62728"  
)

p2 <- DimPlot(mye, reduction = "umap",
              group.by = "disease", cols = disease_cols,     # column in meta.data
              label = FALSE, pt.size = 0.5) + ggtitle(NULL) +
  guides(color = guide_legend(override.aes = list(size = 4))) 

p2

#Proportion of sample type represented across each myeloid subtype 
df_disease <- mye@meta.data %>%
  count(subtype, disease) %>%
  complete(
    subtype,
    disease,
    fill = list(n = 0)
  ) %>%
  group_by(subtype) %>%
  mutate(percent = n / sum(n) * 100)

p_d <- ggplot(df_disease, aes(x = subtype, y = percent, fill = disease)) +
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

#pDC1 vs pDC2 DEG Volcano Plot 
Idents(mye) <- "subtype"  
markers_pDC <- FindMarkers(
  mye,
  ident.1 = "pDC1",
  ident.2 = "pDC2",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

label_genes_pDC1_pDC2 <- c(
  "OAS1", "DNAS1L3", "MPEG1", "TRAF4", "PPP1R14A", 
  "SIGLEC6", "AXL", "SIGLEC1", "CX3CR1", "FGR", "CLEC10A", "CLEC7A",
  "MS4A7", "IFI30", "LYZ", "SERPINB89"
)

volc_df <- markers_pDC %>%
  tibble::rownames_to_column("Gene") %>%
  mutate(
    negLog10P = -log10(p_val),
    signif = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0 ~ "pDC1",
      p_val_adj < 0.05 & avg_log2FC < 0 ~ "pDC2",
      TRUE ~ "NS"
    )
  )

plot_volcano <- ggplot(volc_df, aes(x = avg_log2FC, y = negLog10P)) +
  geom_point(aes(color = signif), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "pDC1" = "red",
    "pDC2" = "blue",
    "NS" = "grey70"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = subset(volc_df, Gene %in% label_genes_pDC1_pDC2),
    aes(label = Gene),
    size = 5,
    box.padding = 0.4,
    max.overlaps = Inf
  ) +
  labs(
    x = "log2 Fold Change (pDC1 vs pDC2)",
    y = "-log10(p-value)",
    color = "Enrichment"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 15),
    legend.position = "right"
  )

plot_volcano

#pDC2 Gene Set Enrichment Analyses - Gene Ontology Pathways 
library(clusterProfiler)
library(org.Hs.eg.db)

pDC2_genes <- rownames(markers_pDC)[markers_pDC$avg_log2FC < -1 & markers_pDC$p_val_adj < 0.05]

ego_pDC2 <- enrichGO(
  gene          = pDC2_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

GO_pDC2 <- GO_pDC2 <- dotplot(ego_pDC2, showCategory = 15) +
  ggtitle("pDC2-enriched pathways")

GO_pDC2
