# -------------------------------
# 0. Cargar librerías necesarias
# -------------------------------
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("sva", quietly = TRUE)) install.packages("sva")
if (!requireNamespace("limma", quietly = TRUE)) install.packages("limma")
library(dplyr)
library(readr)
library(sva)
library(limma)

# -------------------------------
# 1. Leer archivos de expresión
# -------------------------------
gse1 <- read.delim("GSE56421_series_matrix.txt", comment.char = "!", header = TRUE, stringsAsFactors = FALSE)
gse2 <- read.delim("GSE38190_series_matrix.txt", comment.char = "!", header = TRUE, stringsAsFactors = FALSE)
setwd('/sers/katiaavinapadilla/Downloads')
# -------------------------------
# 2. Leer anotaciones de plataformas
# -------------------------------
gpl6244 <- read_csv("GPL6244.csv", show_col_types = FALSE)
colnames(gpl6244)[1:2] <- c("ID", "RefSeq_ID")

gpl570 <- read_csv("GPL570_Clean.csv", show_col_types = FALSE)
colnames(gpl570) <- c("ID", "Symbol", "RefSeq_ID", "Entrez_ID")

# -------------------------------
# 3. Función de anotación y colapso
# -------------------------------
annotate_by_refseq <- function(expr_mat, gpl_df, probe_col = "ID", refseq_col = "RefSeq_ID") {
  gpl_df <- gpl_df[!is.na(gpl_df[[refseq_col]]) & gpl_df[[refseq_col]] != "", ]
  common_probes <- intersect(rownames(expr_mat), gpl_df[[probe_col]])
  expr_mat <- expr_mat[common_probes, , drop = FALSE]
  gpl_df <- gpl_df[gpl_df[[probe_col]] %in% common_probes, ]
  refseq_map <- gpl_df[match(rownames(expr_mat), gpl_df[[probe_col]]), refseq_col]
  refseq_map <- as.character(refseq_map)
  refseq_map <- sapply(strsplit(refseq_map, " /// "), `[`, 1)
  expr_df <- as.data.frame(expr_mat)
  expr_df$RefSeq <- refseq_map
  expr_df <- expr_df[!is.na(expr_df$RefSeq), ]
  expr_df <- expr_df %>%
    group_by(RefSeq) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    as.data.frame()
  rownames(expr_df) <- expr_df$RefSeq
  expr_df$RefSeq <- NULL
  return(expr_df)
}

# -------------------------------
# 4. Preprocesamiento de matrices
# -------------------------------
# GSE56421
rownames(gse1) <- gse1[[1]]
expr1 <- gse1[, -1]
expr1 <- apply(expr1, 2, as.numeric)
rownames(expr1) <- rownames(gse1)
expr1 <- annotate_by_refseq(expr1, gpl6244)

# GSE38190
rownames(gse2) <- gse2[[1]]
expr2 <- gse2[, -1]
expr2 <- apply(expr2, 2, as.numeric)
rownames(expr2) <- rownames(gse2)
expr2 <- annotate_by_refseq(expr2, gpl570)

# -------------------------------
# 5. Combinar y aplicar ComBat
# -------------------------------
common_genes <- intersect(rownames(expr1), rownames(expr2))
expr1 <- expr1[common_genes, ]
expr2 <- expr2[common_genes, ]
combined_expr <- cbind(expr1, expr2)
batch <- c(rep("GSE56421", ncol(expr1)), rep("GSE38190", ncol(expr2)))

# Aplicar ComBat con mod=NULL
combat_expr <- ComBat(dat = as.matrix(combined_expr), batch = batch, par.prior = TRUE, mod = NULL)

# -------------------------------
# 6. Exportar resultado
# -------------------------------
final_df <- as.data.frame(combat_expr)
final_df <- tibble::rownames_to_column(final_df, var = "RefSeq_ID")
write_csv(final_df, "combined_expression_matrix_refseq2.csv")

library(sva)
library(ggplot2)
library(ggfortify)

# Leer matriz combinada
setwd("/Users/katiaavinapadilla/Downloads")
expr <- read.csv("Combat_Corrected_Expression_Matrix_Sept.csv", row.names = 1)

# Crear vector de batches
batch <- c(rep("GSE56421", 9), rep("GSE38190", 11))

# Aplicar ComBat
combat_expr <- ComBat(dat = as.matrix(expr), batch = batch, par.prior = TRUE)

# PCA
pca <- prcomp(t(combat_expr), scale. = TRUE)
autoplot(pca, data = data.frame(batch = batch), colour = "batch")
# Instala paquetes si es necesario
if (!requireNamespace("limma")) install.packages("limma")
if (!requireNamespace("edgeR")) install.packages("edgeR")
if (!requireNamespace("ggplot2")) install.packages("ggplot2")

library(limma)
library(edgeR)
library(ggplot2)

# === 1. Leer matriz corregida ===
expr <- read.csv("Combat_Corrected_Expression_Matrix_Sept.csv", row.names = 1)

# === 2. Definir condiciones manualmente ===
group <- c(
  # GSE56421 (orden correcto)
  "Cultured_Limbal_Epithelium", "Cultured_Limbal_Epithelium", "Cultured_Limbal_Epithelium",
  "Native_Cornea", "Native_Cornea", "Native_Cornea",
  "Cultured_Limbal_Fibroblasts", "Cultured_Limbal_Fibroblasts", "Cultured_Limbal_Fibroblasts",
  # GSE38190
  "Conjunctiva", "Conjunctiva", "Conjunctiva", "Cornea", "Cornea", "Cornea",
  "Conjunctiva", "Limbus", "Limbus", "Limbus", "Limbus"
)
"Limbus_vs_Fibro" = c("GSM932370", "GSM932371", "GSM932372", "GSM932373",        # Limbus
                      "GSM1361191", "GSM1361192", "GSM1361193")   

group <- factor(group)
colnames(expr) <- paste0("S", 1:ncol(expr))  # Nombres simplificados para evitar errores

# === 3. Crear objeto DGEList y normalizar con voom ===
dge <- DGEList(counts = expr)
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

v <- voom(dge, design, plot = TRUE)

# === 4. Ajustar modelo lineal ===
fit <- lmFit(v, design)

# --- 5. Definir contrastes (corrigiendo Limbus_vs_Fibro) ---
contrast.matrix <- makeContrasts(
  Limbus_vs_Cornea = Limbus- Native_Cornea,
  Limbus_vs_CulturedEpi = Limbus - Cultured_Limbal_Epithelium,
  Limbus_vs_Conjuctiva        = Limbus - Conjunctiva,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# --- 6. Exportar resultados y volcanos (con FDR en Y y líneas guía) ---
alpha <- 0.05
lfc_cut <- 0.5

for (contrast_name in colnames(contrast.matrix)) {
  res <- topTable(fit2, coef = contrast_name, number = Inf, adjust.method = "fdr")
  write.csv(res, paste0("DEA_", contrast_name, ".csv"), row.names = TRUE)
  
  # Etiqueta de significancia
  res$signif <- with(res, adj.P.Val < alpha & abs(logFC) >= lfc_cut)
  
  # Guardar solo significativos
  sig <- res[res$signif, , drop = FALSE]
  write.csv(sig, paste0("DEA_", contrast_name, "_SIGNIFICATIVOS.csv"), row.names = TRUE)
  
  # Volcán: usa FDR en eje Y
  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = signif)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("gray", "red")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    labs(title = paste("Volcano Plot:", contrast_name),
         x = "log2 Fold Change", y = "-log10(FDR)", color = "Signif") +
    theme_minimal()
  
  ggsave(paste0("Volcano_", contrast_name, ".png"), p, width = 6, height = 5, dpi = 300)
}
############

for (contrast_name in colnames(contrast.matrix)) {
  # Resultados (con FDR)
  res <- topTable(fit2, coef = contrast_name, number = Inf, adjust.method = "fdr")
  write.csv(res, paste0("DEA_", contrast_name, ".csv"), row.names = TRUE)
  
  # Clasificar: Up / Down / NS
  res$status <- "NS"
  res$status[res$adj.P.Val < alpha & res$logFC >=  lfc_cut] <- "Up"
  res$status[res$adj.P.Val < alpha & res$logFC <= -lfc_cut] <- "Down"
  
  # Contar
  counts <- table(res$status)
  up_n   <- ifelse("Up"   %in% names(counts), counts[["Up"]],   0)
  down_n <- ifelse("Down" %in% names(counts), counts[["Down"]], 0)
  ns_n   <- ifelse("NS"   %in% names(counts), counts[["NS"]],   0)
  
  # Guardar listas y conteos
  write.csv(res[res$status=="Up",   , drop=FALSE],
            paste0("DEA_", contrast_name, "_UP.csv"),   row.names = TRUE)
  write.csv(res[res$status=="Down", , drop=FALSE],
            paste0("DEA_", contrast_name, "_DOWN.csv"), row.names = TRUE)
  write.csv(as.data.frame(counts),
            paste0("DEA_", contrast_name, "_counts_up_down_ns.csv"),
            row.names = FALSE)
  
  # Volcán: azul (Down), gris (NS), rojo (Up), con FDR en Y
  cols <- c("Down" = "#1f77b4", "NS" = "gray80", "Up" = "#d62728")
  
  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
    geom_point(alpha = 0.7, size = 1.2) +
    scale_color_manual(values = cols, breaks = c("Down","NS","Up")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      subtitle = paste0("Up: ", up_n, " | Down: ", down_n,
                        " | FDR < ", alpha, " & |log2FC| ≥ ", lfc_cut),
      x = "log2 Fold Change", y = "-log10(FDR)", color = "Estado"
    ) +
    theme_minimal()
  
  ggsave(paste0("Volcano_", contrast_name, ".png"), p, width = 6, height = 5, dpi = 300)
  
  # Acumular resumen
  summary_counts <- rbind(summary_counts,
                          data.frame(contrast = contrast_name, up = up_n, down = down_n, ns = ns_n))
}

# Guardar resumen global de conteos
write.csv(summary_counts, "DEA_counts_summary.csv", row.names = FALSE)

# Mostrar en consola
print(summary_counts)
###########


realizar_enriquecimiento <- function(file, nombre_contraste) {

  res <- read.csv("DEA_Cornea_vs_Conjunctiva.csv", header = FALSE)
  colnames(res) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  head(res)
  
  res <- read.csv("DEA_Cornea_vs_Conjunctiva.csv", header = FALSE)
  colnames(res) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  
  # Forzar columnas numéricas
  res$logFC <- as.numeric(as.character(res$logFC))
  res$adj.P.Val <- as.numeric(as.character(res$adj.P.Val))
  
  degs <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)

  
  entrez_ids <- bitr(degs$RefSeq_ID, fromType = "REFSEQ", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
    cat("⚠️ No gene could be mapped for", nombre_contraste, "\n")
    return(NULL)
  }
  
  ego <- enrichGO(gene = entrez_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff = 0.05,
                  readable = TRUE)
  
  write.csv(as.data.frame(ego), paste0("GO_enrichment_", nombre_contraste, ".csv"))
  
  pdf(paste0("GO_barplot.pdf"))
  print(barplot(ego, showCategory = 15, title = paste("GO BP:", nombre_contraste)))
  dev.off()
  
  pdf(paste0("GO_dotplot_", nombre_contraste, ".pdf"))
  print(dotplot(ego, showCategory = 15, title = paste("GO BP:", nombre_contraste)))
  dev.off()
  
  return(ego)
}

# 1. Leer matriz de expresión y archivo DEA
expr <- read.csv("Combat_Corrected_Expression_Matrix_Sept.csv", row.names = 1)
dea <- read.csv("DEA_Cornea_vs_Conjunctiva.csv", header = FALSE)
colnames(dea) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

# 2. Filtrar DEGs
dea$logFC <- as.numeric(as.character(dea$logFC))
dea$adj.P.Val <- as.numeric(as.character(dea$adj.P.Val))
degs <- subset(dea, adj.P.Val < 0.05 & abs(logFC) > 1)

# 3. Subconjunto de expresión sólo para esos genes
expr_deg <- expr[rownames(expr) %in% degs$RefSeq_ID, ]

# 4. PCA
expr_t <- t(expr_deg)  # Transponer: filas = muestras
pca <- prcomp(expr_t, scale. = TRUE)

# 5. Crear vector de condiciones para las muestras
group <- c(
  # GSE56421
  "Cultured_Limbal_Epithelium", "Cultured_Limbal_Epithelium", "Cultured_Limbal_Epithelium",
  "Native_Cornea", "Native_Cornea", "Native_Cornea",
  "Cultured_Limbal_Fibroblasts", "Cultured_Limbal_Fibroblasts", "Cultured_Limbal_Fibroblasts",
  # GSE38190
  "Conjunctiva", "Conjunctiva", "Conjunctiva", "Cornea", "Cornea", "Cornea",
  "Conjunctiva", "Limbus", "Limbus", "Limbus", "Limbus"
)

# 6. Gráfico
library(ggplot2)
df_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = group)

ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  labs(title = "PCA based on DEGs: Cornea vs Conjunctiva") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# 1. Leer archivos
expr <- read.csv("combined_expression_matrix_refseq.csv", row.names = 1)
dea <- read.csv("DEA_Cornea_vs_Conjunctiva.csv", header = FALSE)
colnames(dea) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

# 2. Filtrar DEGs
dea$logFC <- as.numeric(as.character(dea$logFC))
dea$adj.P.Val <- as.numeric(as.character(dea$adj.P.Val))
degs <- subset(dea, adj.P.Val < 0.05 & abs(logFC) > 1)

# 3. Subset de expresión solo con los DEGs
expr_deg <- expr[rownames(expr) %in% degs$RefSeq_ID, ]

# 4. Identificar muestras de Cornea y Conjunctiva
# Asegúrate que estos IDs estén bien escritos y ordenados como en la matriz
samples_cornea <- c("GSM1361188", "GSM1361189", "GSM1361190")   # Cornea
samples_conjunctiva <- c("GSM724093", "GSM724094", "GSM724095", "GSM932369")  # Conjunctiva
samples_to_use <- c(samples_cornea, samples_conjunctiva)

# 5. Subset de la matriz solo con esas muestras
expr_sub <- expr_deg[, samples_to_use]

# 6. PCA
expr_t <- t(expr_sub)  # Transponer
pca <- prcomp(expr_t, scale. = TRUE)

# 7. Metadata para condiciones
condition <- factor(c(
  rep("Cornea", length(samples_cornea)),
  rep("Conjunctiva", length(samples_conjunctiva))
))

# 8. Visualización
library(ggplot2)
df_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = condition)

ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  labs(title = "PCA based on DEGs (only Cornea vs Conjunctiva)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Leer matriz y resultados
expr <- read.csv("combined_expression_matrix_refseq.csv", row.names = 1)
dea <- read.csv("DEA_Cornea_vs_Conjunctiva.csv", header = FALSE)
colnames(dea) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

# Filtrar DEGs
dea$logFC <- as.numeric(as.character(dea$logFC))
dea$adj.P.Val <- as.numeric(as.character(dea$adj.P.Val))
degs <- subset(dea, adj.P.Val < 0.05 & abs(logFC) > 1)

# Subset expresión
expr_deg <- expr[rownames(expr) %in% degs$RefSeq_ID, ]
samples <- c("GSM1361188", "GSM1361189", "GSM1361190", "GSM724093", "GSM724094", "GSM724095", "GSM932369")
expr_sub <- expr_deg[, samples]
condition <- factor(c(rep("Cornea", 3), rep("Conjunctiva", 4)))

# PCA
pca <- prcomp(t(expr_sub), scale. = TRUE)
df_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = condition)
ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  labs(title = "PCA: Cornea vs Conjunctiva") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
dea <- read.csv("DEA_Epi_vs_Fibro.csv", header = FALSE)
colnames(dea) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
dea$logFC <- as.numeric(as.character(dea$logFC))
dea$adj.P.Val <- as.numeric(as.character(dea$adj.P.Val))
degs <- subset(dea, adj.P.Val < 0.05 & abs(logFC) > 1)

expr_deg <- expr[rownames(expr) %in% degs$RefSeq_ID, ]
samples <- c("GSM1361185", "GSM1361186", "GSM1361187", "GSM1361191", "GSM1361192", "GSM1361193")
expr_sub <- expr_deg[, samples]
condition <- factor(c(rep("Cultured_Epi", 3), rep("Fibroblasts", 3)))

pca <- prcomp(t(expr_sub), scale. = TRUE)
df_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = condition)
ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  labs(title = "PCA: Cultured Epithelium vs Fibroblasts") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
###
dea <- read.csv("DEA_Limbus_vs_CulturedEpi.csv", header = FALSE)
colnames(dea) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
dea$logFC <- as.numeric(as.character(dea$logFC))
dea$adj.P.Val <- as.numeric(as.character(dea$adj.P.Val))
degs <- subset(dea, adj.P.Val < 0.05 & abs(logFC) > 1)

expr_deg <- expr[rownames(expr) %in% degs$RefSeq_ID, ]
samples <- c("GSM932370", "GSM932371", "GSM932372", "GSM932373", "GSM1361185", "GSM1361186", "GSM1361187")
expr_sub <- expr_deg[, samples]
condition <- factor(c(rep("Limbus", 4), rep("Cultured_Epi", 3)))

pca <- prcomp(t(expr_sub), scale. = TRUE)
df_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Condition = condition)
ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 4) +
  labs(title = "PCA: Limbus vs Cultured Epithelium") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
library(pheatmap)
library(dplyr)
library(readr)
library(pheatmap)
library(dplyr)
library(readr)
library(tibble)
library(pheatmap)
library(dplyr)
library(readr)
library(tibble)

# Diccionario de muestras por grupo
grupo_muestras <- list(
  "Cornea_vs_Conjunctiva" = c("GSM724093","GSM724094","GSM724095","GSM932369",  # Conjunctiva
                              "GSM724096","GSM724097","GSM724098"),             # Cornea
  "Limbus_vs_CulturedEpi" = c("GSM932370","GSM932371","GSM932372","GSM932373",   # Limbus
                              "GSM1361185","GSM1361186","GSM1361187"),           # Cultured Epithelium
  "Epi_vs_Fibro" = c("GSM1361185","GSM1361186","GSM1361187",                     # Cultured Epithelium
                     "GSM1361191","GSM1361192","GSM1361193")                     # Fibroblasts
)

generar_heatmap_top50 <- function(expr_path, dea_path, contrast_name) {
  # 1. Cargar datos
  expr <- read_csv(expr_path) %>% column_to_rownames(var = colnames(.)[1])
  dea <- read_csv(dea_path, col_names = FALSE)
  colnames(dea) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  
  dea$logFC <- as.numeric(dea$logFC)
  dea$adj.P.Val <- as.numeric(dea$adj.P.Val)
  
  # 2. Seleccionar top 25 up y down
  degs <- dea %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)
  top_up <- degs %>% arrange(desc(logFC)) %>% head(25)
  top_down <- degs %>% arrange(logFC) %>% head(25)
  top50 <- bind_rows(top_up, top_down) %>% distinct(RefSeq_ID)
  
  # 3. Filtrar matriz de expresión
  expr_top <- expr[rownames(expr) %in% top50$RefSeq_ID, ]
  muestras_usar <- grupo_muestras[[contrast_name]]
  expr_top <- expr_top[, colnames(expr_top) %in% muestras_usar]
  
  # 4. Z-score y heatmap
  expr_scaled <- t(scale(t(as.matrix(expr_top))))
  pheatmap(expr_scaled,
           cluster_rows = TRUE, cluster_cols = TRUE,
           fontsize_row = 7, show_rownames = TRUE,
           main = paste("Top 50 DEGs -", contrast_name),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}
expr_file <- "combined_expression_matrix_refseq.csv"

generar_heatmap_top50(expr_file, "DEA_Cornea_vs_Conjunctiva.csv", "Cornea_vs_Conjunctiva")
generar_heatmap_top50(expr_file, "DEA_Limbus_vs_CulturedEpi.csv", "Limbus_vs_CulturedEpi")
generar_heatmap_top50(expr_file, "DEA_Epi_vs_Fibro.csv", "Epi_vs_Fibro")

grupo_muestras <- list(
  "Cornea_vs_Conjunctiva" = c("GSM724093","GSM724094","GSM724095","GSM932369",  # Conjunctiva
                              "GSM724096","GSM724097","GSM724098"),             # Cornea
  "Limbus_vs_CulturedEpi" = c("GSM932370","GSM932371","GSM932372","GSM932373",   # Limbus
                              "GSM1361185","GSM1361186","GSM1361187"),           # Cultured Epithelium
  "Epi_vs_Fibro" = c("GSM1361185","GSM1361186","GSM1361187",                     # Cultured Epithelium
                     "GSM1361191","GSM1361192","GSM1361193")                     # Fibroblasts
)
# 1. Diccionario de muestras por grupo (debes ejecutarlo antes)
grupo_muestras <- list(
  "Cornea_vs_Conjunctiva" = c("GSM724093","GSM724094","GSM724095","GSM932369",
                              "GSM724096","GSM724097","GSM724098"),
  "Limbus_vs_CulturedEpi" = c("GSM932370","GSM932371","GSM932372","GSM932373",
                              "GSM1361185","GSM1361186","GSM1361187"),
  "Epi_vs_Fibro" = c("GSM1361185","GSM1361186","GSM1361187",
                     "GSM1361191","GSM1361192","GSM1361193")
)

# 2. Librerías necesarias
library(pheatmap)
library(dplyr)
library(readr)

# 3. Función corregida
generar_heatmap_top50 <- function(expr_path, dea_path, contrast_name) {
  # Leer y preparar expresión
  expr <- read_csv(expr_path)
  rownames(expr) <- expr[[1]]
  expr <- expr[ , -1]
  
  # Leer DEA
  dea <- read_csv(dea_path, col_names = FALSE)
  colnames(dea) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  
  # Filtrar DEGs
  degs <- dea %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1)
  
  top_up <- degs %>% arrange(desc(logFC)) %>% head(25)
  top_down <- degs %>% arrange(logFC) %>% head(25)
  top50 <- bind_rows(top_up, top_down) %>% distinct(RefSeq_ID)
  
  # Filtrar matriz de expresión
  expr_top <- expr[rownames(expr) %in% top50$RefSeq_ID, ]
  muestras_usar <- grupo_muestras[[contrast_name]]
  expr_top <- expr_top[, colnames(expr_top) %in% muestras_usar]
  
  # Z-score y heatmap
  expr_scaled <- t(scale(t(as.matrix(expr_top))))
  pheatmap(expr_scaled,
           cluster_rows = TRUE, cluster_cols = TRUE,
           fontsize_row = 7, show_rownames = TRUE,
           main = paste("Top 50 DEGs -", contrast_name),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
}

#######################################
# -------------------------------
# Cargar librerías necesarias
# -------------------------------
library(pheatmap)
library(dplyr)
library(readr)

# -------------------------------
# Definir muestras por contraste
# -------------------------------
grupo_muestras <- list(
  "Cornea_vs_Conjunctiva" = c("GSM724093","GSM724094","GSM724095","GSM932369",
                              "GSM724096","GSM724097","GSM724098"),
  "Limbus_vs_CulturedEpi" = c("GSM932370","GSM932371","GSM932372","GSM932373",
                              "GSM1361185","GSM1361186","GSM1361187"),
  "Epi_vs_Fibro" = c("GSM1361185","GSM1361186","GSM1361187",
                     "GSM1361191","GSM1361192","GSM1361193")
)

# -------------------------------
# Función para generar heatmaps
# -------------------------------
generar_heatmap_top50 <- function(expr_path, dea_path, contrast_name, output_pdf = FALSE) {
  # Leer matriz de expresión y establecer rownames
  expr <- read_csv(expr_path, show_col_types = FALSE)
  rownames(expr) <- expr[[1]]
  expr <- expr[ , -1]
  
  # Leer archivo de expresión diferencial
  dea <- read_csv(dea_path, show_col_types = FALSE)
  
  # Asignar nombres si están mal o ausentes
  if (!"RefSeq_ID" %in% colnames(dea)) {
    colnames(dea) <- c("RefSeq_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  }
  
  # Asegurar que columnas clave sean numéricas
  dea$logFC <- suppressWarnings(as.numeric(dea$logFC))
  dea$adj.P.Val <- suppressWarnings(as.numeric(dea$adj.P.Val))
  
  # Filtrar genes significativamente expresados
  degs <- dea %>%
    filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1)
  
  # Si no hay suficientes DEGs, salir
  if (nrow(degs) == 0) {
    message("No DEGs found with adj.P.Val < 0.05 and |logFC| > 1.")
    return(NULL)
  }
  
  # Seleccionar top 25 up y 25 down
  top_up <- degs %>% arrange(desc(logFC)) %>% head(25)
  top_down <- degs %>% arrange(logFC) %>% head(25)
  top50 <- bind_rows(top_up, top_down) %>% distinct(RefSeq_ID)
  
  # Filtrar expresión
  expr_top <- expr[rownames(expr) %in% top50$RefSeq_ID, , drop = FALSE]
  muestras_usar <- grupo_muestras[[contrast_name]]
  expr_top <- expr_top[, colnames(expr_top) %in% muestras_usar, drop = FALSE]
  
  # Verificar que haya datos suficientes
  if (nrow(expr_top) < 2 | ncol(expr_top) < 2) {
    message("No hay suficientes genes o muestras para generar heatmap.")
    return(NULL)
  }
  
  # Estandarizar expresión (Z-score por gen)
  expr_scaled <- t(scale(t(as.matrix(expr_top))))
  
  # Guardar a PDF si se indica
  if (output_pdf) {
    pdf(paste0("Heatmap_", contrast_name, ".pdf"), width = 8, height = 10)
  }
  
  # Generar heatmap
  pheatmap(expr_scaled,
           cluster_rows = TRUE, cluster_cols = TRUE,
           fontsize_row = 7, show_rownames = TRUE,
           main = paste("Top 50 DEGs -", contrast_name),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  
  if (output_pdf) {
    dev.off()
    message("Heatmap guardado como: ", paste0("Heatmap_", contrast_name, ".pdf"))
  }
}

expr_file <- "combined_expression_matrix_refseq.csv"

generar_heatmap_top50(expr_file, "DEA_Cornea_vs_Conjunctiva.csv", "Cornea_vs_Conjunctiva", output_pdf = TRUE)
generar_heatmap_top50(expr_file, "DEA_Limbus_vs_CulturedEpi.csv", "Limbus_vs_CulturedEpi", output_pdf = TRUE)
generar_heatmap_top50(expr_file, "DEA_Epi_vs_Fibro.csv", "Epi_vs_Fibro", output_pdf = TRUE)
# Verifica primeros IDs del DEA
dea <- read_csv("DEA_Cornea_vs_Conjunctiva.csv", show_col_types = FALSE)
head(dea[[1]])

# Verifica los rownames del archivo de expresión
expr <- read_csv("combined_expression_matrix_refseq.csv", show_col_types = FALSE)
head(expr[[1]])
expr <- read_csv("combined_expression_matrix_refseq.csv", show_col_types = FALSE)
head(expr[[1]])

# Leer DEA
dea <- read_csv("DEA_Cornea_vs_Conjunctiva.csv", show_col_types = FALSE)
colnames(dea)[1] <- "RefSeq_ID"

# Leer expresión
expr <- read_csv("combined_expression_matrix_refseq.csv", show_col_types = FALSE)
colnames(expr)[1] <- "RefSeq_ID"

# Ver cuántos genes DEA están en expresión
intersectados <- dea$RefSeq_ID[dea$adj.P.Val < 0.05 & abs(as.numeric(dea$logFC)) > 1]
intersectados <- intersect(intersectados, expr$RefSeq_ID)

length(intersectados)  # ¿cuántos de los significativos están en la matriz?
#########################
library(readr)
library(dplyr)
library(pheatmap)

# Define los grupos de muestras por contraste
grupo_muestras <- list(
  "Cornea_vs_Conjunctiva" = c("GSM724093", "GSM724094", "GSM724095", "GSM932369",  # Conjunctiva
                              "GSM724096", "GSM724097", "GSM724098"),             # Cornea
  "Limbus_vs_CulturedEpi" = c("GSM932370", "GSM932371", "GSM932372", "GSM932373",  # Limbus
                              "GSM1361185", "GSM1361186", "GSM1361187"),            # CulturedEpi
  "Epi_vs_Fibro" = c("GSM1361185", "GSM1361186", "GSM1361187",                     # CulturedEpi
                     "GSM1361191", "GSM1361192", "GSM1361193")                     # Fibroblasts
)

# Función para generar heatmaps
generar_heatmap_top50 <- function(expr_path, dea_path, contrast_name, output_pdf = TRUE, guardar_csv = TRUE) {
  # Leer matriz de expresión
  expr <- read_csv(expr_path, show_col_types = FALSE)
  colnames(expr)[1] <- "RefSeq_ID"
  rownames(expr) <- expr$RefSeq_ID
  expr <- expr[, -1]
  
  # Leer DEA
  dea <- read_csv(dea_path, show_col_types = FALSE)
  colnames(dea)[1] <- "RefSeq_ID"
  dea$logFC <- suppressWarnings(as.numeric(dea$logFC))
  dea$adj.P.Val <- suppressWarnings(as.numeric(dea$adj.P.Val))
  
  # Filtrar DEGs significativos y presentes
  degs <- dea %>%
    filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
    filter(RefSeq_ID %in% rownames(expr))
  
  if (nrow(degs) < 10) {
    message("Muy pocos genes para hacer un heatmap.")
    return(NULL)
  }
  
  # Top 25 up y 25 down
  top_up <- degs %>% arrange(desc(logFC)) %>% head(25)
  top_down <- degs %>% arrange(logFC) %>% head(25)
  top50 <- bind_rows(top_up, top_down) %>% distinct(RefSeq_ID)
  
  # Submatriz de expresión
  expr_top <- expr[rownames(expr) %in% top50$RefSeq_ID, , drop = FALSE]
  muestras_usar <- grupo_muestras[[contrast_name]]
  expr_top <- expr_top[, colnames(expr_top) %in% muestras_usar, drop = FALSE]
  
  if (nrow(expr_top) < 2 | ncol(expr_top) < 2) {
    message("No hay suficientes genes o muestras para graficar.")
    return(NULL)
  }
  
  # Anotaciones de columna
  condiciones <- case_when(
    contrast_name == "Cornea_vs_Conjunctiva" ~ c(rep("Conjunctiva", 4), rep("Cornea", 3)),
    contrast_name == "Limbus_vs_CulturedEpi" ~ c(rep("Limbus", 4), rep("CulturedEpi", 3)),
    contrast_name == "Epi_vs_Fibro" ~ c(rep("CulturedEpi", 3), rep("Fibroblasts", 3))
  )
  col_annot <- data.frame(Group = condiciones)
  rownames(col_annot) <- colnames(expr_top)
  
  # Escalar expresión por gen
  expr_scaled <- t(scale(t(as.matrix(expr_top))))
  
  # Guardar CSV
  if (guardar_csv) {
    write.csv(expr_top, paste0("ExpressionMatrix_", contrast_name, "_Top50.csv"))
  }
  
  # Exportar como PDF
  if (output_pdf) {
    pdf(paste0("Heatmap_", contrast_name, ".pdf"), width = 8, height = 10)
  }
  
  # Graficar heatmap
  pheatmap(expr_scaled,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = col_annot,
           fontsize_row = 7,
           show_rownames = TRUE,
           main = paste("Top 50 DEGs -", contrast_name),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  
  if (output_pdf) dev.off()
}
expr_file <- "combined_expression_matrix_refseq.csv"

generar_heatmap_top50(expr_file, "DEA_Cornea_vs_Conjunctiva.csv", "Cornea_vs_Conjunctiva")
generar_heatmap_top50(expr_file, "DEA_Limbus_vs_CulturedEpi.csv", "Limbus_vs_CulturedEpi")
generar_heatmap_top50(expr_file, "DEA_Epi_vs_Fibro.csv", "Epi_vs_Fibro")

# Cargar DEA
dea1 <- read_csv("DEA_Cornea_vs_Conjunctiva.csv", show_col_types = FALSE)
dea2 <- read_csv("DEA_Limbus_vs_CulturedEpi.csv", show_col_types = FALSE)
dea3 <- read_csv("DEA_Epi_vs_Fibro.csv", show_col_types = FALSE)

colnames(dea1)[1] <- "RefSeq_ID"
colnames(dea2)[1] <- "RefSeq_ID"
colnames(dea3)[1] <- "RefSeq_ID"

# Filtrar significativos con logFC y adj.P.Val válidos
dea1$logFC <- as.numeric(dea1$logFC); dea1$adj.P.Val <- as.numeric(dea1$adj.P.Val)
dea2$logFC <- as.numeric(dea2$logFC); dea2$adj.P.Val <- as.numeric(dea2$adj.P.Val)
dea3$logFC <- as.numeric(dea3$logFC); dea3$adj.P.Val <- as.numeric(dea3$adj.P.Val)

sum(dea1$adj.P.Val < 0.05 & abs(dea1$logFC) > 1)  # Cornea vs Conjunctiva
sum(dea2$adj.P.Val < 0.05 & abs(dea2$logFC) > 1)  # Limbus vs CulturedEpi
sum(dea3$adj.P.Val < 0.05 & abs(dea3$logFC) > 1)  # Epi vs Fibro

expr <- read_csv(expr_file, show_col_types = FALSE)
colnames(expr)[1] <- "RefSeq_ID"
rownames(expr) <- expr$RefSeq_ID
expr_ids <- rownames(expr)

# Intersección
sum(dea1$RefSeq_ID %in% expr_ids & dea1$adj.P.Val < 0.05 & abs(dea1$logFC) > 1)
sum(dea2$RefSeq_ID %in% expr_ids & dea2$adj.P.Val < 0.05 & abs(dea2$logFC) > 1)
sum(dea3$RefSeq_ID %in% expr_ids & dea3$adj.P.Val < 0.05 & abs(dea3$logFC) > 1)
colnames(expr)[1:20]  # o head(colnames(expr))
############################

library(readr)
library(dplyr)
library(pheatmap)

grupo_muestras <- list(
  "Cornea_vs_Conjunctiva" = c("GSM724093", "GSM724094", "GSM724095", "GSM932369",  # Conjunctiva
                              "GSM724096", "GSM724097", "GSM724098"),             # Cornea
  "Limbus_vs_CulturedEpi" = c("GSM932370", "GSM932371", "GSM932372", "GSM932373",  # Limbus
                              "GSM1361185", "GSM1361186", "GSM1361187"),            # CulturedEpi
  "Epi_vs_Fibro" = c("GSM1361185", "GSM1361186", "GSM1361187",                     # CulturedEpi
                     "GSM1361191", "GSM1361192", "GSM1361193")                     # Fibroblasts
)

generar_heatmap_top50 <- function(expr_path, dea_path, contrast_name, output_pdf = TRUE, guardar_csv = TRUE) {
  expr <- read_csv(expr_path, show_col_types = FALSE)
  colnames(expr)[1] <- "RefSeq_ID"
  rownames(expr) <- expr$RefSeq_ID
  expr <- expr[, -1]
  
  dea <- read_csv(dea_path, show_col_types = FALSE)
  colnames(dea)[1] <- "RefSeq_ID"
  dea$logFC <- suppressWarnings(as.numeric(dea$logFC))
  dea$adj.P.Val <- suppressWarnings(as.numeric(dea$adj.P.Val))
  
  # Filtrar genes DE con criterio clásico
  degs <- dea %>%
    filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
    filter(RefSeq_ID %in% rownames(expr))
  
  if (nrow(degs) < 10) {
    message("Muy pocos genes para hacer un heatmap.")
    return(NULL)
  }
  
  top_up <- degs %>% arrange(desc(logFC)) %>% head(25)
  top_down <- degs %>% arrange(logFC) %>% head(25)
  top50 <- bind_rows(top_up, top_down) %>% distinct(RefSeq_ID)
  
  expr_top <- expr[rownames(expr) %in% top50$RefSeq_ID, , drop = FALSE]
  
  muestras_solicitadas <- grupo_muestras[[contrast_name]]
  muestras_disponibles <- colnames(expr_top)
  muestras_usar <- intersect(muestras_solicitadas, muestras_disponibles)
  
  if (length(muestras_usar) < 2) {
    message("No hay suficientes muestras disponibles para el contraste: ", contrast_name)
    return(NULL)
  }
  
  expr_top <- expr_top[, muestras_usar, drop = FALSE]
  
  condiciones <- case_when(
    contrast_name == "Cornea_vs_Conjunctiva" ~ c(rep("Conjunctiva", 4), rep("Cornea", 3)),
    contrast_name == "Limbus_vs_CulturedEpi" ~ c(rep("Limbus", 4), rep("CulturedEpi", 3)),
    contrast_name == "Epi_vs_Fibro" ~ c(rep("CulturedEpi", 3), rep("Fibroblasts", 3))
  )
  
  # Alinear anotación a muestras disponibles
  col_annot <- data.frame(Group = condiciones[match(muestras_usar, muestras_solicitadas)])
  rownames(col_annot) <- muestras_usar
  
  expr_scaled <- t(scale(t(as.matrix(expr_top))))
  
  if (guardar_csv) {
    write.csv(expr_top, paste0("ExpressionMatrix_", contrast_name, "_Top50.csv"))
  }
  
  if (output_pdf) {
    pdf(paste0("Heatmap_", contrast_name, ".pdf"), width = 8, height = 10)
  }
  
  pheatmap(expr_scaled,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           annotation_col = col_annot,
           fontsize_row = 7,
           main = paste("Top 50 DEGs -", contrast_name),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
  
  if (output_pdf) dev.off()
}
expr_file <- "Combat_Corrected_Expression_Matrix_Sept.csv"

generar_heatmap_top50(expr_file, "DEA_Cornea_vs_Conjunctiva.csv", "Cornea_vs_Conjunctiva")
generar_heatmap_top50(expr_file, "DEA_Limbus_vs_CulturedEpi.csv", "Limbus_vs_CulturedEpi")
generar_heatmap_top50(expr_file, "DEA_Epi_vs_Fibro.csv", "Epi_vs_Fibro")


######
generar_heatmap_top50(expr_file, "DEA_Limbus_vs_Fibro.csv", "Limbus_vs_Fibro")

