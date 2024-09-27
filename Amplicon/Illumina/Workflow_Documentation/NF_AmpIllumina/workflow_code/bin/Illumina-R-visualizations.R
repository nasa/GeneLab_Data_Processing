#!/usr/bin/env Rscript
pdf(file = NULL)
library(vegan)
library(tidyverse)
library(dendextend)
library(phyloseq)
library(DESeq2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(grid)

##################################################################################
## R visualization script for Illumina paired-end amplicon data                 ##
##################################################################################
# This script is automatically executed as part of the Snakemake workflow when the run_workflow.py --visualizations TRUE argument is used.
# This script can also be manually executed using processed data from the workflow 

# Store command line args as variables #
args <- commandArgs(trailingOnly = TRUE)
runsheet_file <- paste0(args[1])
sample_info <- paste0(args[2])
counts <- paste0(args[3])
taxonomy <- paste0(args[4])
plots_dir <- paste0(args[5])
output_prefix <- paste0(args[6])
assay_suffix <- paste(args[7])
########################################

RColorBrewer_Palette <- "Set1"

# Runsheet read1 path/filename column name
read1_path_colname <- 'read1_path'
# Runsheet read1 suffix column name
raw_R1_suffix_colname <- 'raw_R1_suffix'
# Runsheet groups column name
groups_colname <- 'groups'
# Runsheet colors column name
color_colname <- 'color'

####################
# Helper functions #
####################

# Identify the matching rows by removing suffix from basename of file
remove_suffix <- function(path, suffix) {
  file_name <- basename(path)
  sub(suffix, "", file_name)
}

# Remove the longest common prefix from the sample names (only used for visualizations)
longest_common_prefix <- function(strs) {
  if (length(strs) == 1) return(strs)
  
  prefix <- strs[[1]]
  for (str in strs) {
    while (substring(str, 1, nchar(prefix)) != prefix) {
      prefix <- substr(prefix, 1, nchar(prefix) - 1)
    }
  }
  
  return(prefix)
}
remove_common_prefix <- function(strs) {
  prefix <- longest_common_prefix(strs)
  sapply(strs, function(x) substr(x, nchar(prefix) + 1, nchar(x)))
}

# Adust cex based on number of samples
adjust_cex <- function(num_samples, start_samples = 40, end_samples = 150, default_cex = 1, min_cex = 0.6) {
  slope <- (min_cex - default_cex) / (end_samples - start_samples)
  
  new_cex <- default_cex + slope * (num_samples - start_samples)
  
  adjusted_cex <- max(min(new_cex, default_cex), min_cex)
  
  return(adjusted_cex)
}


# Extract legend from a plot
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend 
} 

###########################################
# Read in data, create output directories #
###########################################

# Assign the directory paths to variables
beta_diversity_out_dir <- file.path(plots_dir, "beta_diversity")
alpha_diversity_out_dir <- file.path(plots_dir, "alpha_diversity")
taxonomy_out_dir <- file.path(plots_dir, "taxonomy")
de_out_dir <- file.path(plots_dir, "da")

abundance_out_dir <- file.path(de_out_dir, "differential_abundance")
volcano_out_dir <- file.path(de_out_dir, "volcano")

# List of all directory variables
out_dirs <- list(plots_dir, beta_diversity_out_dir, alpha_diversity_out_dir, taxonomy_out_dir, de_out_dir, abundance_out_dir, volcano_out_dir)

# Loop through each directory path to check and create if necessary
for (dir_path in out_dirs) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

# Read in processed data
runsheet <- as.data.frame(read.table(file = runsheet_file, 
                                     header = TRUE, sep = ","))
row.names(runsheet) <- runsheet$'Sample.Name'
runsheet$'Sample.Name' <- NULL

count_tab <- read.table(file = counts, 
                        header = TRUE, row.names = 1, sep = "\t")
tax_tab <- read.table(file = taxonomy, 
                      header = TRUE, row.names = 1, sep = "\t")
# Use only samples listed in sample_info, which should correspond to the file names
sample_names <- readLines(sample_info)
deseq2_sample_names <- make.names(sample_names, unique = TRUE)


# Check if the runsheet uses links instead of local file paths
# Extract what is after 'file=' and before any '&' or other URL parameters
uses_links <- any(grepl("^(http|https)://|genelab-data.ndc.nasa.gov", runsheet[[read1_path_colname]]))
if (uses_links) {
    # Use rownames as basenames if links are used
    runsheet$basename <- rownames(runsheet)
} else {
    # Remove extensions from filenames in runsheet for local file paths
    runsheet$basename <- mapply(remove_suffix, runsheet[[read1_path_colname]], runsheet[[raw_R1_suffix_colname]])
}

# Make the basenames DESeq2 compatible, add temporary s_ prefix to fix bugs caused by basenames starting w/ number
runsheet$basename <- paste0("s_", runsheet$basename)
runsheet$basename <- make.names(runsheet$basename, unique = TRUE)
runsheet$basename <- sub("^s_", "", runsheet$basename)

# Subset runsheet and count tab to only include samples in sample_info
runsheet <- runsheet[runsheet$basename %in% deseq2_sample_names, ]
count_tab <- count_tab[, colnames(count_tab) %in% runsheet$basename]

# Order runsheet based on the groups column
runsheet <- runsheet[order(runsheet[[groups_colname]]), ]

# Reorder count_tab columns to match the order in the runsheet
count_tab <- count_tab[, runsheet$basename]

# Rename runsheet row names
rownames(runsheet) <- runsheet$basename

if (!identical(rownames(runsheet), colnames(count_tab))) {
  stop("The read file names in the runsheet do not match the colnames of count_tab.")
}

# Keep only genes with at least 1 count
count_tab <- count_tab[rowSums(count_tab) > 0, ]
count_tab_vst <- count_tab

# Check if every gene has a 0 in the row, add +1 pseudocount for VST, not ideal but fixes VST for sparse counts if needed
if (all(apply(count_tab_vst, 1, any))) {
  count_tab_vst <- count_tab_vst + 1
}

# Create VST normalized counts matrix
deseq_counts <- DESeqDataSetFromMatrix(countData = count_tab_vst, 
                                       colData = runsheet, 
                                       design = ~1)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)


###########################################
# Create plots, save them in output dirs #
###########################################

# Add colors to runsheet
num_colors <- length(unique(runsheet[[groups_colname]]))

# List of RColorBrewer Palette lengths
RColorBrewer_Palette_lengths <- c(Accent=8, Dark2=8, Paired=12, Pastel1=9, Pastel2=8, Set1=9, Set2=8, Set3=12)

# Check if number of colors exceeds the limit of the selected RColorBrewer_Palette palette
if (num_colors > RColorBrewer_Palette_lengths[RColorBrewer_Palette]) {
    # If so, reate a custom palette with more colors
    custom_palette <- colorRampPalette(brewer.pal(RColorBrewer_Palette_lengths[RColorBrewer_Palette], RColorBrewer_Palette))(num_colors)
    colors <- custom_palette
} else {
    # Else just use the standard RColorBrewer_Palette palette
    colors <- brewer.pal(num_colors, RColorBrewer_Palette)
}



group_colors <- setNames(colors, unique(runsheet[[groups_colname]]))
runsheet <- runsheet %>%
  mutate(!!color_colname := group_colors[.data[[groups_colname]]])


########
## Save original par settings
##   Par may be temporarily changed for plotting purposes and reset once the plotting is done
original_par <- par(no.readonly = TRUE)
options(preferRaster=TRUE) # use Raster when possible to avoid antialiasing artifacts in images


width_in_inches <- 11.1
height_in_inches <- 8.33
dpi <- 300
width_in_pixels <- width_in_inches * dpi
height_in_pixels <- height_in_inches * dpi


# Hierarchical Clustering
sample_info_tab <- runsheet[, c(groups_colname, color_colname)]

# Add short group names and legend text column to sample_info
sample_info_tab$short_groups <- as.integer(factor(sample_info_tab[[groups_colname]], levels = unique(sample_info_tab[[groups_colname]])))

group_levels <- unique(sample_info_tab[[groups_colname]])
short_group_labels <- sprintf("%d: %s", seq_along(group_levels), group_levels)
names(short_group_labels) <- group_levels
sample_info_tab$short_group_labels <- short_group_labels[sample_info_tab[[groups_colname]]]
colors_vector <- unique(setNames(sample_info_tab[[color_colname]], sample_info_tab$short_group_labels))

euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(d = euc_dist, method = "ward.D2")

# Color the dendrogram sample labels by group using dendxtend
euc_dend <- as.dendrogram(euc_clust, h = .1)
dend_cols <- sample_info_tab[[color_colname]][order.dendrogram(euc_dend)]
labels_colors(euc_dend) <- dend_cols


##########
default_cex = 1
# Lower cex if over 40 samples to prevent names from crashing on plot
default_cex <- adjust_cex(length(rownames(sample_info_tab)))


# Set for 11x8 plot margins, else try ggdendrogram
space_available <- height_in_inches/5.4

longest_name <- rownames(sample_info_tab)[which.max(nchar(rownames(sample_info_tab)))]

calculate_max_cex <- function(longest_name, space_avail) {
  # Define weights for lower case letters, periods, else
  lower_case_weight <- 0.10
  other_char_weight <- 0.15
  dot_weight <- 0.02  # Weight for the period character
  
  # Calculate weights in longest sample name
  char_weights <- sapply(strsplit(longest_name, "")[[1]], function(char) {
    if (char == ".") {
      return(dot_weight)
    } else if (grepl("[a-z]", char)) {
      return(lower_case_weight)
    } else {
      return(other_char_weight)
    }
  })
  
  average_weight <- mean(char_weights)
  
  # Calculate the maximum cex that fits the space using the average weight
  n = nchar(longest_name)
  max_cex <- space_avail / (n * average_weight)
  
  return(max_cex)
}

max_cex <- calculate_max_cex(longest_name, space_available)
dendro_cex <- min(max_cex, default_cex)


legend_groups <- unique(sample_info_tab$groups)
legend_colors <- unique(sample_info_tab$color)
num_unique_groups <- length(legend_groups)
legend_cex <- ifelse(num_unique_groups > 5, 1 / (num_unique_groups / 5), 1)

png(file.path(beta_diversity_out_dir, paste0(output_prefix, "dendrogram_by_group", assay_suffix, ".png")),
    width = width_in_pixels,
    height = height_in_pixels,
    res = dpi)
par(mar = c(10.5, 4.1, 0.6 , 2.1))
euc_dend %>% set("labels_cex", dendro_cex) %>% plot(ylab = "VST Euc. dist.") 
par(xpd=TRUE)
legend("bottom", inset = c(0, -.34), legend = legend_groups, fill = legend_colors, bty = 'n', cex = legend_cex)
dev.off()
par(original_par)






# making a phyloseq object with our transformed table
vst_count_phy <- otu_table(object = vst_trans_count_tab, taxa_are_rows = TRUE)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
vst_physeq

# generating a PCoA with phyloseq
vst_pcoa <- ordinate(physeq = vst_physeq, method = "PCoA", distance = "euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues

# Calculate the percentage of variance
percent_variance <- eigen_vals / sum(eigen_vals) * 100

betadisper(d = euc_dist, group = sample_info_tab$groups) %>% anova()
adonis_res <- adonis2(formula = euc_dist ~ sample_info_tab$groups)
r2_value <- adonis_res$R2[1]
prf_value <- adonis_res$`Pr(>F)`[1]

label_PC1 <- sprintf("PC1 [%.1f%%]", percent_variance[1])
label_PC2 <- sprintf("PC2 [%.1f%%]", percent_variance[2])


# Save unlabeled PCoA plot
ordination_plot_u <- plot_ordination(vst_physeq, vst_pcoa, color = "groups") + 
  geom_point(size = 1) + 
  labs( 
    x = label_PC1,
    y = label_PC2,
    col = "Groups"
  ) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  scale_color_manual(values = unique(sample_info_tab[[color_colname]][order(sample_info_tab[[groups_colname]])]),
                     labels = unique(sample_info_tab$short_group_labels[order(sample_info_tab[[groups_colname]])])) +
  theme_bw() + theme(legend.position = "bottom",  text = element_text(size = 15, ),
                     legend.direction = "vertical",
                     legend.justification = "center",
                     legend.box.just = "center",
                     legend.title.align = 0.5) +
  annotate("text", x = Inf, y = -Inf, label = paste("R2:", toString(round(r2_value, 3))), hjust = 1.1, vjust = -2, size = 4)+
  annotate("text", x = Inf, y = -Inf, label = paste("Pr(>F)", toString(round(prf_value,4))), hjust = 1.1, vjust = -0.5, size = 4)+ ggtitle("PCoA")
ggsave(filename=file.path(beta_diversity_out_dir, paste0(output_prefix, "PCoA_without_labels", assay_suffix, ".png")), plot=ordination_plot_u, width = 11.1, height = 8.33, dpi = 300)
# Save labeled PCoA plot
ordination_plot <- plot_ordination(vst_physeq, vst_pcoa, color = "groups") + 
  geom_point(size = 1) + 
  labs(
    col = "Groups", 
    x = label_PC1,
    y = label_PC2
  ) + 
  geom_text(aes(label = rownames(sample_info_tab)), show.legend = FALSE, hjust = 0.3, vjust = -0.4, size = 4) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  scale_color_manual(values = unique(sample_info_tab[[color_colname]][order(sample_info_tab[[groups_colname]])]),
                     labels = unique(sample_info_tab$short_group_labels[order(sample_info_tab[[groups_colname]])])) +
  theme_bw() + theme(legend.position = "bottom",  text = element_text(size = 15, ),
                     legend.direction = "vertical",
                     legend.justification = "center",
                     legend.box.just = "center",
                     legend.title.align = 0.5) +
  annotate("text", x = Inf, y = -Inf, label = paste("R2:", toString(round(r2_value, 3))), hjust = 1.1, vjust = -2, size = 4)+
  annotate("text", x = Inf, y = -Inf, label = paste("Pr(>F)", toString(round(prf_value,4))), hjust = 1.1, vjust = -0.5, size = 4)+ ggtitle("PCoA")
ggsave(filename=file.path(beta_diversity_out_dir, paste0(output_prefix, "PCoA_w_labels", assay_suffix, ".png")), plot=ordination_plot, width = 11.1, height = 8.33, dpi = 300)
########################

#4. Alpha diversity

# 4a. Rarefaction curves

p <- rarecurve(x = t(count_tab), step = 100, col = sample_info_tab[[color_colname]], 
               lwd = 2, ylab = "ASVs", label = FALSE, tidy = TRUE)

sample_info_tab_names <- tibble::rownames_to_column(sample_info_tab, var = "Site")
p <- p %>%
  left_join(sample_info_tab_names, by = "Site")

rareplot <- ggplot(p, aes(x = Sample, y = Species, group = Site, color = groups)) + 
  geom_line() + 
  scale_color_manual(values = unique(sample_info_tab[[color_colname]][order(sample_info_tab[[groups_colname]])]),
                   labels = unique(sample_info_tab$short_group_labels[order(sample_info_tab[[groups_colname]])]),
                   breaks = unique(sample_info_tab[[groups_colname]])) +
  labs(x = "Number of Sequences", y = "Number of ASVs", col = "Groups") + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 15),
        legend.direction = "vertical",
        legend.justification = "center",
        legend.box.just = "center",
        legend.title.align = 0.5,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")) +
  guides(color = guide_legend(title = "Groups"))
ggsave(filename = file.path(alpha_diversity_out_dir, paste0(output_prefix, "rarefaction_curves", assay_suffix, ".png")), plot=rareplot, width = 8.33, height = 8.33, dpi = 300)

# 4b. Richness and diversity estimates

# create a phyloseq object similar to how we did above in step 3B, only this time also including our taxonomy table:
count_tab_phy <- otu_table(count_tab, taxa_are_rows = TRUE)
tax_tab_phy <- tax_table(as.matrix(tax_tab))
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)


calculate_text_size <- function(num_samples, start_samples = 25, min_size = 3) {
  max_size = 11  # Maximum size for up to start_samples
  slope = -0.15

  if (num_samples <= start_samples) {
    return(max_size)
  } else {
    # Calculate the current size with the hardcoded slope
    current_size = max_size + slope * (num_samples - start_samples)
    
    # Ensure the size doesn't go below the minimum
    return(max(current_size, min_size))
  }
}

richness_sample_label_size <- calculate_text_size(length(rownames(sample_info_tab)))

richness_plot <- plot_richness(ASV_physeq, color = "groups", measures = c("Chao1", "Shannon")) + 
  scale_color_manual(values = unique(sample_info_tab[[color_colname]][order(sample_info_tab[[groups_colname]])]),
                     labels = unique(sample_info_tab$short_group_labels[order(sample_info_tab[[groups_colname]])])) + 
  theme_bw() +labs(x = "Samples",
                   color = "Groups") +
  theme(
    text = element_text(size = 15),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title.align = 0.5,
    axis.text.x = element_text(angle = 90,
                               size = richness_sample_label_size,
                               vjust = 0.5,  # Vertically center the text
                               hjust = 1)
  )
ggsave(filename = file.path(alpha_diversity_out_dir, paste0(output_prefix, "richness_and_diversity_estimates_by_sample", assay_suffix, ".png")), plot=richness_plot, width = 11.1, height = 8.33, dpi = 300)

richness_by_group <- plot_richness(ASV_physeq, x = "groups", color = "groups", measures = c("Chao1", "Shannon")) +
  scale_color_manual(values = unique(sample_info_tab[[color_colname]][order(sample_info_tab[[groups_colname]])]),
                     labels = unique(sample_info_tab$short_group_labels[order(sample_info_tab[[groups_colname]])])) + 
  scale_x_discrete(labels = unique(sample_info_tab$short_groups[order(sample_info_tab[[groups_colname]])])) +
  theme_bw() + labs(colors = "Groups",
                    x = "Groups") +
  theme(
    text = element_text(size = 15),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title.align = 0.5,
    legend.title = element_blank()
  ) 
ggsave(filename = file.path(alpha_diversity_out_dir, paste0(output_prefix, "richness_and_diversity_estimates_by_group", assay_suffix, ".png")), plot=richness_by_group, width = 11.1, height = 8.33, dpi = 300)

# Extract legend from unlabeled pca plot, also save it as its own plot
legend <- g_legend(ordination_plot)
grid.newpage()
grid.draw(legend)
legend_filename <- file.path(plots_dir, paste0(output_prefix, "color_legend", assay_suffix, ".png"))
increment <- ifelse(length(unique(sample_info_tab$groups)) > 9, ceiling((length(unique(sample_info_tab$groups)) - 9) / 3), 0)
legend_height <- 3 + increment
ggsave(legend_filename, plot = legend, device = "png", width = 11.1, height = legend_height, dpi = 300)


# 5. Taxonomic summaries
# Calculate new plot height of legend is taller than ~ default / 4 

# Get height of legend
heights_in_cm <- sapply(legend$heights, function(h) {
  if (is.null(h)) {
    # Assign 0 to null heights
    return(unit(0, "cm"))
  } else {
    # Convert the unit to centimeters
    return(convertUnit(h, "cm", valueOnly = TRUE))
  }
})
legend_height_in_cm <- sum(heights_in_cm)
legend_height_in_inches <- legend_height_in_cm / 2.54

taxonomy_plots_height <- width_in_inches
required_height_in_inches <- 4 * legend_height_in_inches + 0.5
if (legend_height_in_inches < required_height_in_inches) {
  # Adjust width_in_inches to make sure the legend fits
  taxonomy_plots_height <- required_height_in_inches
}

proportions_physeq <- transform_sample_counts(ASV_physeq, function(ASV) ASV / sum(ASV))
proportions_physeq@sam_data$short_groups <- as.character(proportions_physeq@sam_data$short_groups)

relative_phyla <- plot_bar(proportions_physeq, x = "short_groups", fill = "phylum") + 
  theme_bw() + theme(text = element_text(size = 9)) + labs(x = "Groups")
plot_layout <- grid.layout(nrow = 2, heights = unit(c(3, 1), "null"))
grid.newpage()
pushViewport(viewport(layout = plot_layout))
print(relative_phyla, vp = viewport(layout.pos.row = 1))
pushViewport(viewport(layout.pos.row = 2))
grid.draw(legend)
upViewport(0)
grid_image <- grid.grab()
ggsave(filename = file.path(taxonomy_out_dir, paste0(output_prefix, "relative_phyla", assay_suffix, ".png")), grid_image, width = height_in_inches, height = taxonomy_plots_height, dpi = 500)

relative_classes <- plot_bar(proportions_physeq, x = "short_groups", fill = "class") + 
  theme_bw() + theme(text = element_text(size = 9)) + labs(x = "Groups")
plot_layout <- grid.layout(nrow = 2, heights = unit(c(3, 1), "null"))
grid.newpage()
pushViewport(viewport(layout = plot_layout))
print(relative_classes, vp = viewport(layout.pos.row = 1))
pushViewport(viewport(layout.pos.row = 2))
grid.draw(legend)
upViewport(0)
grid_image <- grid.grab()
ggsave(filename = file.path(taxonomy_out_dir, paste0(output_prefix, "relative_classes", assay_suffix, ".png")), plot=grid_image, width = height_in_inches, height = taxonomy_plots_height, dpi = 500)


#samplewise taxonomy

proportions_physeq <- transform_sample_counts(ASV_physeq, function(ASV) ASV / sum(ASV))
# Calculate the number of samples from the count_tab
num_samples <- ncol(count_tab)
# Scaling function: starts at 1x width for 25 samples
scaling_factor <- (num_samples - 40) / (200 - 40) * (5 - 1) + 1
scaling_factor <- max(1, min(scaling_factor, 5))

# samplewise phyla
samplewise_phylum <- plot_bar(proportions_physeq, fill = "phylum") +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 6)) # Rotate and resize x-axis labels

ggsave(filename = file.path(taxonomy_out_dir, paste0(output_prefix, "samplewise_relative_phyla", assay_suffix, ".png")), 
       plot = samplewise_phylum, 
       width = height_in_inches * scaling_factor,
       height = taxonomy_plots_height)

samplewise_classes <- plot_bar(proportions_physeq, fill = "class") +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5, size = 6)) # Rotate and resize x-axis labels

ggsave(filename = file.path(taxonomy_out_dir, paste0(output_prefix, "samplewise_relative_classes", assay_suffix, ".png")), 
       plot = samplewise_classes, 
       width = height_in_inches * scaling_factor,
       height = taxonomy_plots_height)

# 6 Statistically testing for differences

#### pairwise comparisons
unique_groups <- unique(runsheet$groups)
deseq_obj <- phyloseq_to_deseq2(physeq = ASV_physeq, design = ~groups)

# add pseudocount if any 0 count samples are present
if (sum(colSums(counts(deseq_obj)) == 0) > 0) {
  count_data <- counts(deseq_obj) + 1 
  
  count_data <- as.matrix(apply(count_data, 2, as.integer))
  rownames(count_data) <- rownames(counts(deseq_obj))
  colnames(count_data) <- colnames(counts(deseq_obj))
  counts(deseq_obj) <- count_data
}
# https://rdrr.io/bioc/phyloseq/src/inst/doc/phyloseq-mixture-models.R 
deseq_modeled <- tryCatch({
  # Attempt to run DESeq
  DESeq(deseq_obj)
}, error = function(e) {
  message("Error encountered in DESeq, applying alternative method for size factor estimation...")
  
  # Define the geometric mean function
  gm_mean = function(x, na.rm=TRUE) {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(deseq_obj), 1, gm_mean)
  
  # Apply the alternative size factor estimation method
  deseq_obj <- estimateSizeFactors(deseq_obj, geoMeans=geoMeans)
  
  # Call DESeq again with alternative geom mean size est
  DESeq(deseq_obj)
})

# save final differential abundance counts, individual group comparison results

write.table(counts(deseq_modeled, normalized=TRUE), file = file.path(de_out_dir, paste0(output_prefix, "normalized_counts", assay_suffix, ".tsv")), sep="\t", row.names=TRUE, quote=FALSE)
# make the volcanoplot
plot_comparison <- function(group1, group2) {
  plot_width_inches = 11.1
  plot_height_inches = 8.33
  
  deseq_res <- results(deseq_modeled, contrast = c("groups", group1, group2))
  norm_tab <- counts(deseq_modeled, normalized = TRUE) %>% data.frame()
  
  volcano_data <- as.data.frame(deseq_res)
  
  p_val <- 0.1
  volcano_data <- volcano_data[!is.na(volcano_data$padj), ]
  volcano_data$significant <- volcano_data$padj <= p_val #also logfc cutoff?
  
  ######Long x-axis label adjustments##########
  x_label <- paste("Log2 Fold Change\n(",group1," vs ",group2,")")
  label_length <- nchar(x_label)
  max_allowed_label_length = plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- paste("Log2 Fold Change\n\n(", group1, "\n vs \n", group2, ")", sep="")
  }
  #######################################

  # ASVs promoted in space on right, reduced on left
  p <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("black", "red"), labels=c(paste0("padj > ", p_val), paste0("padj \u2264 ", p_val))) +
    theme_bw() +
    labs(title="Volcano Plot",
         x=x_label,
         y="-Log10 P-value",
         color=paste0("")) +
    theme(legend.position="top")
  
  # label points and plot
  top_points <- volcano_data %>%
    arrange(padj) %>%
    filter(significant) %>%
    head(10)
  
  volcano_plot <- p + geom_text_repel(data=top_points, aes(label=row.names(top_points)), size=3)
  ggsave(filename=file.path(volcano_out_dir, paste0(output_prefix, "volcano_", gsub(" ", "_", group1), "_vs_", gsub(" ", "_", group2), ".png")),
         plot=volcano_plot,
         width = plot_width_inches, height = plot_height_inches, dpi = 300)
  
  write.csv(deseq_res, file = file.path(abundance_out_dir, paste0(output_prefix, gsub(" ", "_", group1), "_vs_", gsub(" ", "_", group2), ".csv")))
}


# setting up pairwise comparisons and running
comparisons <- expand.grid(group1 = unique_groups, group2 = unique_groups)
comparisons <- subset(comparisons, group1 != group2)

apply(comparisons, 1, function(pair) plot_comparison(pair['group1'], pair['group2']))


dev.off()
