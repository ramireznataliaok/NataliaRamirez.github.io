data <- read.delim("Data/GSE11414_series_matrix.txt.gz",
                   header = TRUE, 
                   sep = "\t",
                   stringsAsFactors = FALSE)

readLines("Data/GSE11414_series_matrix.txt.gz", n = 20)
lines <- readLines("Data/GSE11414_series_matrix.txt.gz")
grep("!series_matrix_table_begin", lines)
grep("!series_matrix_table_end", lines)
lines[grep("!series_matrix_table_begin", lines) + 1]

# Path to the GEO series matrix
file_path <- "Data/GSE11414_series_matrix.txt.gz"

# Read the complete compressed file
lines <- readLines(file_path)

# Find the expression table boundaries
table_start <- grep("!series_matrix_table_begin", lines) + 1
table_end <- grep("!series_matrix_table_end", lines) - 1

# Read only the expression table
expression_data <- read.delim(
  text = lines[table_start:table_end],
  header = TRUE,
  sep = "\t",
  quote = "\"",
  check.names = FALSE
)

# Inspect the imported data
dim(expression_data)
head(expression_data)

grep("!Sample_title", lines, value = TRUE)

# Install packages
install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("GEOquery")
install.packages(c(
  "ggplot2",
  "pheatmap",
  "EnhancedVolcano"))


# Confirm BiocManager is installed
library(BiocManager)

# Set the Bioconductor release compatible with R 4.6
BiocManager::install(version = "3.23", ask = FALSE)

# Install the Bioconductor packages
BiocManager::install(
  c("limma", "GEOquery", "EnhancedVolcano"),
  ask = FALSE,
  update = FALSE
)
force = TRUE

# Load packages
library(limma)
library(GEOquery)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Convert probe IDs to row names
rownames(expression_data) <- expression_data$ID_REF

# Remove the probe ID column and create a numeric matrix
expression_matrix <- as.matrix(expression_data[, -1])

# Confirm the matrix is numeric
storage.mode(expression_matrix) <- "numeric"

# Check the dimensions
dim(expression_matrix)

# Preview part of the matrix
expression_matrix[1:5, 1:3]

# Create the sample metadata
sample_metadata <- data.frame(
  sample = colnames(expression_matrix),
  cell_type = factor(
    c("HOB", "HOB", "MG63", "MG63", "U2OS", "U2OS"),
    levels = c("HOB", "MG63", "U2OS")
  ),
  replicate = c(1, 2, 1, 2, 1, 2)
)

sample_metadata

# Check normalization of data
boxplot(
  expression_matrix,
  las = 2,
  main = "Expression distributions before additional normalization",
  ylab = "Expression intensity"
)

summary(expression_matrix)

# Transpose the matrix so samples are rows
pca_result <- prcomp(
  t(expression_matrix),
  center = TRUE,
  scale. = FALSE
)

# Calculate percentage of variance explained
percent_variance <- round(
  100 * pca_result$sdev^2 / sum(pca_result$sdev^2),
  1
)

# Create a PCA data frame
pca_data <- data.frame(
  sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  cell_type = sample_metadata$cell_type,
  replicate = sample_metadata$replicate
)

# Plot PCA
pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = cell_type,
    label = sample
  )
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.8,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA of GSE11414 expression data",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    color = "Cell type"
  ) +
  theme_minimal()
pca_plot

ggsave(
  filename = "Figures/PCA_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

save(
  expression_data,
  expression_matrix,
  sample_metadata,
  pca_result,
  pca_data,
  pca_plot,
  file = "Data/analysis_objects.RData"
)

# Confirm the sample order
colnames(expression_matrix)
sample_metadata

identical(
  colnames(expression_matrix),
  as.character(sample_metadata$sample)
)

# Create the design matrix
design <- model.matrix(
  ~0 + cell_type,
  data = sample_metadata
)
colnames(design) <- levels(sample_metadata$cell_type)

design

# Fit linear model independently to each probe
fit <- lmFit(expression_matrix, design)

# Define the comparisons
contrast_matrix <- makeContrasts(
  MG63_vs_HOB = MG63 - HOB,
  U2OS_vs_HOB = U2OS - HOB,
  MG63_vs_U2OS = MG63 - U2OS,
  levels = design
)
contrast_matrix

#Apply contrasts and empirical Bayes statistics
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

#Extract complete results
results_MG63 <- topTable(
  fit2,
  coef = "MG63_vs_HOB",
  number = Inf,
  adjust.method = "BH",
  sort.by = "p"
)
results_U2OS <- topTable(
  fit2,
  coef = "U2OS_vs_HOB",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
results_MG63_vs_U2OS <- topTable(
  fit2,
  coef = "MG63_vs_U2OS",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Inspect first few rows
head(results_MG63)
head(results_U2OS)

#Define statistically significant probes
#adjusted p-value < 0.05
#absolute log2 fold change >= 1
sig_MG63 <- subset(
  results_MG63,
  adj.P.Val < 0.05 & abs(logFC) >= 1
)
sig_U2OS <- subset(
  results_U2OS,
  adj.P.Val < 0.05 & abs(logFC) >= 1
)
#Count the significant probes
nrow(sig_MG63)
nrow(sig_U2OS)

#Count upregulated and downregulated probes separately
table(
  MG63_direction = ifelse(
    sig_MG63$logFC > 0,
    "Upregulated",
    "Downregulated"
  )
)

table(
  U2OS_direction = ifelse(
    sig_U2OS$logFC > 0,
    "Upregulated",
    "Downregulated"
  )
)

#Add probe IDs as normal column
results_MG63$Probe_ID <- rownames(results_MG63)
results_U2OS$Probe_ID <- rownames(results_U2OS)

sig_MG63$Probe_ID <- rownames(sig_MG63)
sig_U2OS$Probe_ID <- rownames(sig_U2OS)

#Save the tables
write.csv(
  results_MG63,
  "Data/MG63_vs_HOB_all_results.csv",
  row.names = FALSE
)

write.csv(
  results_U2OS,
  "Data/U2OS_vs_HOB_all_results.csv",
  row.names = FALSE
)

write.csv(
  sig_MG63,
  "Data/MG63_vs_HOB_significant_results.csv",
  row.names = FALSE
)

write.csv(
  sig_U2OS,
  "Data/U2OS_vs_HOB_significant_results.csv",
  row.names = FALSE
)
nrow(sig_MG63)
nrow(sig_U2OS)

#Next step is to annotate the probe IDs with gene symbols
# Install and load the annotation package
BiocManager::install(
  c("AnnotationDbi", "hugene10sttranscriptcluster.db"),
  ask = FALSE,
  update = FALSE
)

library(AnnotationDbi)
library(hugene10sttranscriptcluster.db)

# Annotate your probe IDs
probe_ids <- as.character(results_MG63$Probe_ID)

annotation <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "PROBEID"
)
results_MG63$SYMBOL <- AnnotationDbi::mapIds(
  hugene10sttranscriptcluster.db,
  keys = results_MG63$Probe_ID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)
results_MG63$GENENAME <- AnnotationDbi::mapIds(
  hugene10sttranscriptcluster.db,
  keys = results_MG63$Probe_ID,
  column = "GENENAME",
  keytype = "PROBEID",
  multiVals = "first"
)
results_MG63$ENTREZID <- AnnotationDbi::mapIds(
  hugene10sttranscriptcluster.db,
  keys = results_MG63$Probe_ID,
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)
results_U2OS$SYMBOL <- AnnotationDbi::mapIds(
  hugene10sttranscriptcluster.db,
  keys = results_U2OS$Probe_ID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)
results_U2OS$GENENAME <- AnnotationDbi::mapIds(
  hugene10sttranscriptcluster.db,
  keys = results_U2OS$Probe_ID,
  column = "GENENAME",
  keytype = "PROBEID",
  multiVals = "first"
)
results_U2OS$ENTREZID <- AnnotationDbi::mapIds(
  hugene10sttranscriptcluster.db,
  keys = results_U2OS$Probe_ID,
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

# Define the annotated tables
annotated_MG63 <- results_MG63
annotated_U2OS <- results_U2OS

#Confirm the row counts
nrow(annotated_MG63)
nrow(annotated_U2OS)

head(annotation)
dim(annotation)

# Merge the annotation with both comparisons
colnames(results_MG63)
head(rownames(results_MG63))
results_MG63$Probe_ID <- rownames(results_MG63)
results_U2OS$Probe_ID <- rownames(results_U2OS)

head(results_MG63$Probe_ID)
head(results_U2OS$Probe_ID)

length(results_MG63$Probe_ID)
nrow(results_MG63)

#Annotate using the IDs
probe_ids <- as.character(results_MG63$Probe_ID)

annotation <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "PROBEID"
)

#Ensure the ID columns are characters
results_MG63$Probe_ID <- as.character(results_MG63$Probe_ID)
results_U2OS$Probe_ID <- as.character(results_U2OS$Probe_ID)
annotation$PROBEID <- as.character(annotation$PROBEID)

#Merge
annotated_MG63 <- merge(
  results_MG63,
  annotation,
  by.x = "Probe_ID",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)

annotated_U2OS <- merge(
  results_U2OS,
  annotation,
  by.x = "Probe_ID",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)

# Check whether gene symbols were successfully added
head(annotated_MG63[, c("Probe_ID", "SYMBOL", "GENENAME", "logFC", "adj.P.Val")])

# Check annotation coverage
sum(!is.na(annotated_MG63$SYMBOL))
sum(is.na(annotated_MG63$SYMBOL))

# Recreate the significant annotated tables
sig_annotated_MG63 <- subset(
  annotated_MG63,
  adj.P.Val < 0.05 & abs(logFC) >= 1
)

sig_annotated_U2OS <- subset(
  annotated_U2OS,
  adj.P.Val < 0.05 & abs(logFC) >= 1
)

# Check how many significant probes have gene symbols
nrow(sig_annotated_MG63)
sum(!is.na(sig_annotated_MG63$SYMBOL))

nrow(sig_annotated_U2OS)
sum(!is.na(sig_annotated_U2OS$SYMBOL))

#Count upregulated and downregulated probes
table(
  ifelse(
    sig_annotated_MG63$logFC > 0,
    "Upregulated",
    "Downregulated"
  )
)

table(
  ifelse(
    sig_annotated_U2OS$logFC > 0,
    "Upregulated",
    "Downregulated"
  )
)
# Store the direction in the tables
sig_annotated_MG63$Direction <- ifelse(
  sig_annotated_MG63$logFC > 0,
  "Upregulated",
  "Downregulated"
)

sig_annotated_U2OS$Direction <- ifelse(
  sig_annotated_U2OS$logFC > 0,
  "Upregulated",
  "Downregulated"
)
# Remove missing gene symbols before counting
length(unique(na.omit(sig_annotated_MG63$SYMBOL)))
length(unique(na.omit(sig_annotated_U2OS$SYMBOL)))

#View the strongest differentially expressed genes
# Highest expression in MG63 relative to HOB:
head(
  sig_annotated_MG63[
    order(sig_annotated_MG63$logFC, decreasing = TRUE),
    c("Probe_ID", "SYMBOL", "GENENAME", "logFC", "adj.P.Val")
  ],
  10
)
# Lowest expression in MG63 relative to HOB:
head(
  sig_annotated_MG63[
    order(sig_annotated_MG63$logFC),
    c("Probe_ID", "SYMBOL", "GENENAME", "logFC", "adj.P.Val")
  ],
  10
)
#Repeat for U2OS:
head(
  sig_annotated_U2OS[
    order(sig_annotated_U2OS$logFC, decreasing = TRUE),
    c("Probe_ID", "SYMBOL", "GENENAME", "logFC", "adj.P.Val")
  ],
  10
)

head(
  sig_annotated_U2OS[
    order(sig_annotated_U2OS$logFC),
    c("Probe_ID", "SYMBOL", "GENENAME", "logFC", "adj.P.Val")
  ],
  10
)
#Prepare labels for volcano plots
#Use gene symbols where available
annotated_MG63$Plot_Label <- ifelse(
  is.na(annotated_MG63$SYMBOL),
  "",
  annotated_MG63$SYMBOL
)

annotated_U2OS$Plot_Label <- ifelse(
  is.na(annotated_U2OS$SYMBOL),
  "",
  annotated_U2OS$SYMBOL
)

#Create volcano plots
library(EnhancedVolcano)

volcano_MG63 <- EnhancedVolcano(
  annotated_MG63,
  lab = annotated_MG63$Plot_Label,
  x = "logFC",
  y = "adj.P.Val",
  title = "MG63 versus HOB",
  subtitle = "GSE11414 differential expression",
  xlab = expression(Log[2] ~ "fold change"),
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.5,
  labSize = 3,
  drawConnectors = TRUE,
  max.overlaps = 15
)

volcano_MG63

volcano_U2OS <- EnhancedVolcano(
  annotated_U2OS,
  lab = annotated_U2OS$Plot_Label,
  x = "logFC",
  y = "adj.P.Val",
  title = "U2OS versus HOB",
  subtitle = "GSE11414 differential expression",
  xlab = expression(Log[2] ~ "fold change"),
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 1.5,
  labSize = 3,
  drawConnectors = TRUE,
  max.overlaps = 15
)
volcano_U2OS

# Save the plots
ggsave(
  "Figures/volcano_MG63_vs_HOB.png",
  plot = volcano_MG63,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "Figures/volcano_U2OS_vs_HOB.png",
  plot = volcano_U2OS,
  width = 9,
  height = 7,
  dpi = 300
)
# Save final annotated tables
write.csv(
  annotated_MG63,
  "Data/MG63_vs_HOB_annotated_results.csv",
  row.names = FALSE
)

write.csv(
  annotated_U2OS,
  "Data/U2OS_vs_HOB_annotated_results.csv",
  row.names = FALSE
)

write.csv(
  sig_annotated_MG63,
  "Data/MG63_vs_HOB_significant_annotated.csv",
  row.names = FALSE
)

write.csv(
  sig_annotated_U2OS,
  "Data/U2OS_vs_HOB_significant_annotated.csv",
  row.names = FALSE
)
