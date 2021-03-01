library("recount3")

###INITIAL MODIFICATIONS

## Obtaining projects
human_projects <- available_projects()

## Title:
## Single-nucleus RNA sequencing of human cortex affected by multiple sclerosis

## Creating RSE object with specific info
rse_gene_SRP173190	 <- create_rse(
  subset(
    human_projects,
    project == "SRP173190" & project_type == "data_sources"
  )
)


## Extracting read counts
assay(rse_gene_SRP173190, "counts") <- compute_read_counts(rse_gene_SRP173190)

## Making info easier to handle
rse_gene_SRP173190 <- expand_sra_attributes(rse_gene_SRP173190)


## Exploring interest columns

colData(rse_gene_SRP173190)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP173190)))
]

## Saving SE object just in case
rse_gene_SRP173190_unfiltered <- rse_gene_SRP173190


## Exploring differneces between groups
table(rse_gene_SRP173190$sra_attribute.disease_state)
barplot(table(rse_gene_SRP173190$sra_attribute.disease_state), main = "Groups")

## Exploring expression means
summary(rowMeans(assay(rse_gene_SRP173190, "counts")))

## Filtering of low expression genes
expr_means <- rowMeans(assay(rse_gene_SRP173190, "counts"))
rse_gene_SRP173190 <- rse_gene_SRP173190[expr_means > 0.1, ]

## Exploring assign percentage
hist(rse_gene_SRP173190$`recount_qc.gene_fc.all_%`)
## No filtering because assignment percentage is quite good

## How much did we lost
round(nrow(rse_gene_SRP173190) / nrow(rse_gene_SRP173190_unfiltered) * 100, 2)

### NORMALIZING DATA

library("edgeR")
dge <- DGEList(
  counts = assay(rse_gene_SRP173190, "counts"),
  genes = rowData(rse_gene_SRP173190)
)
dge <- calcNormFactors(dge)

## Exploring our proposed statistic model

library("ExploreModelMatrix")

data <- data.frame(disease_state = colData(rse_gene_SRP173190)$sra_attribute.disease_state, source = colData(rse_gene_SRP173190)$sra_attribute.source_name)
vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = data,
  designFormula = ~ disease_state + source,
  textSizeFitted = 4
)

## Veamos las imágenes
cowplot::plot_grid(plotlist = vd$plotlist)

ExploreModelMatrix(data)
## Proposed model (~disease_state + source) is full rank


### DE ANALYSIS

## Exploring differences between groups

library("iSEE")
iSEE::iSEE(rse_gene_SRP173190)

## At first sight there are no big differences


library("ggplot2")
library("ggsignif")

## Checking if differences are significant in total % of assignment

df <- data.frame(disease_state = rse_gene_SRP173190$sra_attribute.disease_state, gene_assign = rse_gene_SRP173190$`recount_qc.gene_fc.all_%`)

ggplot(df, aes(x = disease_state, y = gene_assign)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assignment %") +
  xlab("Disease state") +
  geom_signif(comparisons=list(c("MS", "healthy control")),map_signif_level=TRUE)

## Assigning our already defined model to variable

mod <- model.matrix(~ sra_attribute.disease_state + sra_attribute.source_name,
                    data = colData(rse_gene_SRP173190)
)


library("limma")

## Visualizing SD of the expression of the genes
vGene <- voom(dge, mod, plot = FALSE)

## Creates linear regression model and calculates p-values
eb_results <- eBayes(lmFit(vGene))


## Summary of expression results
de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP173190),
  sort.by = "none"
)
dim(de_results)

head(de_results)

## how many genes are actually expressing differentially
table(de_results$adj.P.Val < 0.05)

## Visualizing expression differences in both groups
## valores positivos en eje Y significan mayor expresion en healthy y negativos mayor expresion en MS
## Hay más valores positivos
plotMA(eb_results, coef = 2)

## Highlighting 5 genes with lower p-values
volcanoplot(eb_results, coef = 2, highlight = 5, names = de_results$gene_name)


##Extracting most important genes

exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

## creating data.frame
df <- as.data.frame(colData(rse_gene_SRP173190)[, c("sra_attribute.disease_state", "sra_attribute.source_name")])
colnames(df) <- c("DiseaseState", "Tissue")


## Creating pheatmap
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df,
  fontsize_row=4
)


library("RColorBrewer")

## MDS por estado de enfermedad
df <- as.data.frame(colData(rse_gene_SRP173190)[, c("sra_attribute.disease_state", "sra_attribute.source_name")])
plotMDS(vGene$E, labels = df$sra_attribute.disease_state, col = rep(c("red", "blue"), each = 25))

## MDS por tejido

col.tissue <- df$sra_attribute.source_name

col.tissue[col.tissue==unique(df$sra_attribute.source_name)] <- c("blue", "red", "purple", "orange", "black")

col.tissue[45:49] <- c("blue", "red", "purple", "orange", "black")

plotMDS(vGene$E, labels = df$sra_attribute.source_name, col = col.tissue)
