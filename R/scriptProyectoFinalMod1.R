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
barplot(table(rse_gene_SRP173190$sra_attribute.disease_state))

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


### DE ANALYSIS

## Exploring differences between groups

library("ggplot2")

df <- data.frame(disease_state = rse_gene_SRP173190$sra_attribute.disease_state, gene_assign = rse_gene_SRP173190$`recount_qc.gene_fc.all_%`)

ggplot(df, aes(x = disease_state, y = gene_assign)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assignment %") +
  xlab("Disease state")

##
