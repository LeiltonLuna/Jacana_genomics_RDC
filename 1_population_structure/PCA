library(SNPRelate)
library(gdsfmt)
library(ggplot2)

# Define input files
vcf_file <- "jacana__N27_filtered.vcf.gz"
gds_file <- "jacana_N27_filtered.gds"

# Convert VCF to GDS
if (!file.exists(gds_file)) {
  snpgdsVCF2GDS(vcf_file, gds_file)
}

# Open GDS file
genofile <- snpgdsOpen(gds_file)

# Load metadata containing Sample and taxon columns
metadata <- read.csv("jacana_pca_N27.csv",
                     stringsAsFactors = FALSE)

# Identify non-hybrid samples present in the GDS
all_samples <- read.gdsn(
  index.gdsn(genofile, "sample.id")
)

samples_keep <- all_samples[
  all_samples %in% metadata$Sample
]

# Perform PCA using ONLY non-hybrid individuals
jacana_pca <- snpgdsPCA(
  genofile,
  sample.id = samples_keep,
  autosome.only = FALSE
)

# Create PCA dataframe
pca_df <- data.frame(
  Sample = jacana_pca$sample.id,
  PC1 = jacana_pca$eigenvect[, 1],
  PC2 = jacana_pca$eigenvect[, 2]
)

# Add taxon information
pca_df$taxon <- metadata$taxon[
  match(pca_df$Sample, metadata$Sample)
]

# Calculate percentage of variance explained
pc1_var <- round(jacana_pca$varprop[1] * 100, 1)
pc2_var <- round(jacana_pca$varprop[2] * 100, 1)

# Define shapes and colors
shape_values <- c(
  "jacana_spinosa" = 21,
  "jacana_jacana_hypo" = 24,
  "jacana_jacana_jac" = 24
)

fill_colors <- c(
  "jacana_spinosa" = "yellow",
  "jacana_jacana_hypo" = "black",
  "jacana_jacana_jac" = "red"
)

# Plot PCA
ggplot(pca_df, aes(
  x = PC1,
  y = PC2,
  fill = taxon,
  shape = taxon
)) +
  geom_hline(
    yintercept = 0,
    color = "gray40",
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    color = "gray40",
    linetype = "dashed"
  ) +
  geom_point(size = 4, color = "black") +
  scale_shape_manual(values = shape_values) +
  scale_fill_manual(values = fill_colors) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    fill = "Taxon",
    shape = "Taxon"
  ) +
  theme_classic() +
  theme(
    panel.background = element_rect(
      fill = "white",
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA
    ),
    axis.line = element_line(color = "black")
  )

# Close GDS file
snpgdsClose(genofile)
