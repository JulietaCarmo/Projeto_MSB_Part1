library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggsci)
library(DESeq2)


#barcodes <- read.table("C:/Users/anaan/OneDrive - Universidade de Aveiro/Attachments/Desktop/OneDrive - Universidade de Coimbra/Biologia de sistemas/Parte 1/projeto/dados/barcodes.tsv", 
#                       header = FALSE)

#features <- read.table("C:/Users/anaan/OneDrive - Universidade de Aveiro/Attachments/Desktop/OneDrive - Universidade de Coimbra/Biologia de sistemas/Parte 1/projeto/dados/features.tsv", 
#                       header = FALSE)

#matrix <- readMM("C:/Users/anaan/OneDrive - Universidade de Aveiro/Attachments/Desktop/OneDrive - Universidade de Coimbra/Biologia de sistemas/Parte 1/projeto/dados/matrix.mtx")

barcodes <- read.table("/home/julieta-carmo/Documents/BCM/Part1/Project_Irina/Data/barcodes.tsv", 
                       header = FALSE)

features <- read.table("/home/julieta-carmo/Documents/BCM/Part1/Project_Irina/Data/features.tsv", 
                       header = FALSE)

matrix <- readMM("/home/julieta-carmo/Documents/BCM/Part1/Project_Irina/Data/matrix.mtx")


# Lables da matrix
rownames(matrix) <- features[, 1]
colnames(matrix) <- barcodes[, 1]


seurat_obj <- CreateSeuratObject(counts = matrix, project = "BAL", min.cells = 3, #genes expressos em pelo menos 3 células
                                 min.features = 200)
seurat_obj #21347 features across 20188 samples within 1 assay 

#metadata = read.table("C:/Users/anaan/OneDrive - Universidade de Aveiro/Attachments/Desktop/OneDrive - Universidade de Coimbra/Biologia de sistemas/Parte 1/projeto/dados/BAL_alexandria_structured_metadata3.txt", 
#                      header = TRUE, row.names = 1, sep = "\t",)
metadata = read.table("/home/julieta-carmo/Documents/BCM/Part1/Project_Irina/Data/RBAL_alexandria_structured_metadata3.txt", 
                      header = TRUE, row.names = 1, sep = "\t",)

seurat_obj <- AddMetaData(
  object = seurat_obj,
  metadata = metadata[, c("CellTypeAnnotations", "Smoking_Status", "donor_id")]
)

View(seurat_obj@meta.data)

length(unique(seurat_obj@meta.data$CellTypeAnnotations)) #temos 16 tipos celulares diferentes
# Remover o underscore 
#seurat_obj@meta.data$CellTypeAnnotations <- gsub("_", " ", seurat_obj@meta.data$CellTypeAnnotations)
#seurat_obj@meta.data$Smoking_Status <- gsub("_", " ", seurat_obj@meta.data$Smoking_Status)



unique(seurat_obj@meta.data$CellTypeAnnotations)

seurat_obj@meta.data %>%
  ggplot(aes(x= Smoking_Status, fill = CellTypeAnnotations)) + 
  geom_bar(position = "fill", color = "black", linewidth = 0.3) +
  theme_minimal() +  
  labs(title = "Cell Type Proportions", y = "Proportion", x = "Smoking_status",
       fill = "Cell Type")+
  scale_fill_d3(palette = "category20")

##################coisas da aula
#controlo de qualidade 
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") #não têm genes com este padrão
head(seurat_obj@meta.data)

VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

seurat_obj = subset(seurat_obj, subset = nCount_RNA< 150000 & nFeature_RNA < 7700)
seurat_obj # 21347 features across 20042 samples within 1 assay
#removemos cerca de 140 células

#normalização
seurat_obj = NormalizeData(seurat_obj, 
                           normalization.method = "LogNormalize", 
                           scale.factor = 1e4)
head(seurat_obj@assays$RNA@layers$data) # mostra a expressão normalizada/relativa de cada gene por célula
head(seurat_obj@assays$RNA@layers$counts) # mostra a expressão absoluta de cada gene por célula

#genes que variam mais
seurat_obj = FindVariableFeatures(seurat_obj, 
                                  selection.method = "vst", 
                                  nfeatures = 2000)
VariableFeatures(seurat_obj)
VariableFeaturePlot(seurat_obj)
top10 <- head(VariableFeatures(seurat_obj), 10)
top10
VariableFeaturePlot(seurat_obj) %>%
  LabelPoints(points = top10, repel= T, xnudge = 0, ynudge = 0)

#info sobre o top 3 de genes mais variados 
#MKI67 expression reflects the cellular proliferation rate
#Encondes Interleukin-8 (IL8) is a member of the CXC chemokine family These small basic heparin-binding proteins are proinflammatory and primarily mediate the activation and migration of neutrophils into tissue from peripheral blood
# top2A Topoisomerase II from eukaryotic cells catalyzes the relaxation of supercoiled DNA molecules, catenation, decatenation, knotting, and unknotting of circular DNA.

#Redução das dimensões

#scaling 

# STEP 4 - Scale
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

#PCA apenas nos features que variam mais, para ter menos ruido
seurat_obj = RunPCA(seurat_obj, features= VariableFeatures(seurat_obj))
DimPlot(seurat_obj, reduction = "pca")
DimHeatmap(seurat_obj, dims= 1, cells = 200)

#escolher o numero de componentes
ElbowPlot(seurat_obj)
seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20) # 20 vectors are chosen
JackStrawPlot(seurat_obj, dims = 1:20) #plot and choose based on p-value
#manter 20 pcs

saveRDS(seurat_obj, file = "3_04_26.rds")

####

seurat_obj = FindNeighbors(seurat_obj, dims= 1:20)
seurat_obj= FindClusters(seurat_obj, resolution = 0.9)
Idents(seurat_obj)
#vizualiar os cluster obtidos com o UMAP

seurat_obj = RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, reduction = "umap", group.by = c("CellTypeAnnotations", "seurat_clusters"), 
        split.by = "Smoking_Status", cols = DiscretePalette(n = 23, palette = "polychrome"))
seurat_obj@meta.data[c("CellTypeAnnotations", "seurat_clusters")]
DimPlot(seurat_obj, group.by = c("donor_id", "Smoking_Status"), cols = DiscretePalette(9, palette = "polychrome"))


#####################
# 1. Split by donor
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$donor_id)

# 2. Process
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# 3. Integrate
seurat_obj <- IntegrateLayers(seurat_obj, 
                              method = CCAIntegration,
                              orig.reduction = "pca",
                              new.reduction = "integrated.cca")

# 4. JOIN LAYERS HERE ← before clustering!
seurat_obj <- JoinLayers(seurat_obj)

# 5. Cluster on integrated reduction
seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated.cca", dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.9)
seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.cca", dims = 1:20)


DimPlot(seurat_obj, group.by = c("donor_id", "Smoking_Status"))


DimPlot(seurat_obj, group.by = "CellTypeAnnotations", split.by = "Smoking_Status", 
                 cols = DiscretePalette(16, palette = "polychrome"))


saveRDS(seurat_obj, file = "4_04_26.rds")
seurat_obj <- readRDS("4_04_26.rds")

############### Expressão diferencial entre tipo de células para as condições
Idents(seurat_obj) <- "Smoking_Status"
FindMarkers(seurat_obj, ident.1 = "SMOKER", ident.2 = "NON_SMOKER")

#para considerar a expressão dos individuos, e não de células isoladamente
#junta se a expressão de cada tipo de célula para cada individuo
aggr <- AggregateExpression(seurat_obj, assays = "RNA", return.seurat = T, 
                            group.by = c("CellTypeAnnotations", "Smoking_Status", "donor_id"))

View(aggr@assays$RNA$data)
tail(Cells(aggr))



Idents(aggr) <- "Smoking_Status"

condicoes <- FindMarkers(object = aggr, 
                            ident.1 = "SMOKER", 
                            ident.2 = "NON-SMOKER",
                            test.use = "wilcox_limma")
head(condicoes, n = 5)

aggr$celltype.smoke <- paste(aggr$CellTypeAnnotations, aggr$Smoking_Status, sep = "_")

Idents(aggr) <- "CellTypeAnnotations"

talvez <- FindAllMarkers(object = aggr, 
                    group_by = "Smoking_Status",
                    test.use = "wilcox_limma", 
                    logfc.threshold = 0.50,
                    min.pct = 0.25)
head(talvez, n = 5)

#acho que é igual. penso que está a comparar contra todas as outras as possibilidades de células e nnão só do mesmo tipo
Idents(aggr) = "celltype.smoke"
levels(Idents(aggr))

talvez.2 = FindAllMarkers(object = aggr,
                          test.use = "wilcox_limma", 
                          logfc.threshold = 0.50, 
                          min.pct = 0.25)

head(talvez.2, n = 5)


#quero comparar cada tipo de célula entre as duas combinações

Idents(aggr) = "CellTypeAnnotations"

levels(Idents(aggr))

tipo_cel = unique(aggr@meta.data$CellTypeAnnotations)
tipo_smo = unique(aggr@meta.data$Smoking_Status)

markers.celulas = list()
for (cel in tipo_cel) {
  #ident.1 = paste(cel, tipo_smo[1])
  #ident.2 = paste(cel, tipo_smo[2])

  markers.celulas[[cel]] = FindMarkers(object = aggr,
                            test.use = "DESeq2", #"wilcox_limma", 
                            logfc.threshold = 0.50, 
                            min.pct = 0.25, 
                            ident.1 = "SMOKER",
                            ident.2= "NON-SMOKER",
                            group.by = "Smoking_Status", 
                            subset.ident = cel
                            )
}

names(markers.celulas)

head(markers.celulas$"BAL-Monocytes")



sig_results <- lapply(markers.celulas, function(x) {
  x[x$p_val_adj < 0.05, ]
})
sig_results$"BAL-Monocytes"


sapply(sig_results, nrow)


sig_results[["T-Cells"]]
sig_results[["Macrophage1"]]

saveRDS(seurat_obj, file = "1_05_26.rds")

###############################
################ GSE for one cell type
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(europepmc)
library(org.Hs.eg.db)

organism = "org.Hs.eg.db"
# Escolhi este tipo de celulas pq foi o que tinhas usado antes
celltype <- "BAL-Monocytes"

# differential-expression table for BAL-Monocytes
de <- markers.celulas[[celltype]]
de_5 <- head(de, n=5)

# inspect columns first
head(de)
colnames(de)

# build ranked gene list
# for DESeq2 results in Seurat, avg_log2FC may differ by version;
# sometimes you may have avg_log2FC, sometimes log2FoldChange-like output
gene_list <- de_5$avg_log2FC
names(gene_list) <- rownames(de_5)

# remove NA and sort
gene_list <- na.omit(gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

head(gene_list)
head(names(gene_list))

# run GSEA GO
gse <- gseGO(geneList = gene_list, 
             OrgDb = org.Hs.eg.db, 
             keyType = "GO", 
             ont = "BP", 
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             pAdjustMethod = "none")

head(gse@result)

saveRDS(gse, file = "gse_Monocytes.rds")
gse_Mono <- readRDS("gse_Monocytes.rds")

# plots
dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)

x <- pairwise_termsim(gse)
emapplot(x, showCategory = 10)

ridgeplot(gse) + labs(x = "enrichment distribution")

# literature trend for top 3 terms
terms <- gse@result$Description[1:3]
pmcplot(terms, 2010:2025, proportion = FALSE)







#############################################
celltype <- "Macrophage1"

# differential-expression table for BAL-Monocytes
de <- markers.celulas[[celltype]]
de_5 <- head(de, n=5)

# inspect columns first
head(de)
colnames(de)

# build ranked gene list
# for DESeq2 results in Seurat, avg_log2FC may differ by version;
# sometimes you may have avg_log2FC, sometimes log2FoldChange-like output
gene_list <- de_5$avg_log2FC
names(gene_list) <- rownames(de_5)

# remove NA and sort
gene_list <- na.omit(gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

head(gene_list)
head(names(gene_list))

# run GSEA GO
gse <- gseGO(geneList = gene_list, 
             OrgDb = org.Hs.eg.db, 
             keyType = "GO", 
             ont = "BP", 
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             pAdjustMethod = "none")

head(gse@result)

# plots
dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)

x <- pairwise_termsim(gse)
emapplot(x, showCategory = 10)

ridgeplot(gse) + labs(x = "enrichment distribution")

# literature trend for top 3 terms
terms <- gse@result$Description[1:3]
pmcplot(terms, 2010:2025, proportion = FALSE)

saveRDS(gse, file = "gse_Macrophage1.rds")

