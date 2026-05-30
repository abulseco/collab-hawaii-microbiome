# HAWAII MICROBIOME PROJECT
# Sediment samples from Derek Esibill
# Summer 2024
# Phyloseq Analysis (Exploratory)
# Last updated: 12/13/24 by ANB

# PLEASE READ: This code assumes you have run the raw sequences through
# HAWAII-MICROBIOME-DADA2.R. This is a preliminary analysis only as the samples likely
# need to be resequenced. Denoising on forward reads only.

## LIBRARIES-------------------------------------------
library(phyloseq); library(data.table); library(ggplot2);  library(dplyr); library(RColorBrewer)

## DATA PREPARATION------------------------------------
# Set working directory
setwd("data-for-phyloseq")
list.files()

## FORMATTING FUNCTIONS--------------------------------
pretty.theme <- function(){
  theme_bw() +
    theme(axis.text.x=element_text(size = 18, color="black"),
          axis.text.y=element_text(size = 18, color="black"),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          text=element_text(size=18),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA),
          axis.line = element_line(colour = "black"))
}

mat = read.table("ASVs_counts.txt", header = TRUE, sep = "\t", row.names = 1)
tax = read.table("ASVs_taxonomy.txt", header = TRUE, sep = "\t", row.names = 1)
meta = read.csv("mapping-file.csv", header = TRUE)

mat = as.matrix(mat)
tax = as.matrix(tax)

OTU = otu_table(mat, taxa_are_rows = TRUE)
TAX = tax_table(tax)
META = sample_data(meta)
sample_names(OTU)
sample_names(TAX)
tail(META)
head(TAX)
head(OTU)

## Troubleshooting-------
# Issues with getting the names to match up
# Check sample names in each component
otu_sample_names <- sample_names(OTU)
meta_sample_names <- sample_names(META)
tax_sample_names <- sample_names(TAX)  # If you have a taxonomic table

# Print sample names for inspection
print(otu_sample_names)
print(meta_sample_names)


phy = phyloseq(OTU,TAX,META)
sample_names(phy)

## Get to know your phyloseq object-----------
ntaxa(phy) # 741
nsamples(phy) # Check that it is what you expect - 5
sample_variables(phy) # Should match your metadata file
taxa_names(phy)
sample_names(phy)
sample_sums(phy)

# Look at statistics related to your data
num.reads = sum(sample_sums(phy))
lowest.sam = sort(sample_sums(phy)) 
mean.seq = mean(sample_sums(phy))  
std.seq = sd(sample_sums(phy))/sqrt(30) 
med.seq = median(sample_sums(phy))
phy.summary <- data.frame(num.reads, mean.seq, std.seq, med.seq)
phy.summary     

seq.dt = data.table(as(sample_data(phy), "data.frame"),
                    TotalReads = sample_sums(phy), keep.rownames = TRUE)
seq.dt # information per sample after filtering out taxa 

# Visualize dataset properties
readsumsdf = data.frame(nreads = sort(taxa_sums(phy), TRUE), sorted = 1:ntaxa(phy), type = "Nodes")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(phy),TRUE),sorted = 1:nsamples(phy), type = "Samples"))
sample_sums(phy)
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x=sorted, y = nreads)) + geom_bar(stat="identity")
p + ggtitle(title) + 
  facet_wrap(~type, 1, scales = "free") + 
  scale_y_log10() 

# Filtering-----------------------------------------
# Get rid of Mitochondria  
phy %>%
  subset_taxa(Family != "Mitochondria") -> phy.f
phy # BEFORE filtering: 741 and 5 samples
phy.f # AFTER filtering: 417 and 5 samples 

# Filter out low abundances
lowabundnames = filter_taxa(phy.f, function(x) sum(x) > 0.1)
phy.f.nolow = prune_taxa(lowabundnames, phy.f) 
ntaxa(phy.f.nolow) # still 417 remain - can continue with phy.f 

# Calculate relative abundance
per.f = transform_sample_counts(phy.f, function (x) x/sum(x)*100)
per.f.nolow = transform_sample_counts(phy.f.nolow, function (x) x/sum(x)*100)

## Plotting----------
# Ordination
# Tip, use the following to see different color schemes in R Color Brewer
display.brewer.all(colorblindFriendly = TRUE)

# Ordination--------
BC_distance <- phyloseq::distance(phy.f, "bray") 
bcOrd <- ordinate(phy.f.nolow, "PCoA", BC_distance)
plot_scree(bcOrd)

# Difference between stream and riparian
p1 <- plot_ordination(phy.f, bcOrd, type = "samples") +
  geom_point(aes(fill = Location), shape = 24, color = "black", stroke = 1.5, size = 6, alpha = 0.9) +
  # geom_text(aes(label = Location), vjust = 1, size = 3, nudge_x = 0.025) + 
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")) +
  pretty.theme() +
  labs(fill = "Collection Site")
p1

# Barplot of all phyla--------
physeq_phylum <- tax_glom(phy.f, taxrank = "Phylum")

# Transform to relative abundance
physeq_phylum_rel <- transform_sample_counts(physeq_phylum, function(x) x / sum(x))

# Melt the phyloseq object
physeq_phylum_melt <- psmelt(physeq_phylum_rel)

# Optional: Sort the data by abundance for better visualization
physeq_phylum_melt <- physeq_phylum_melt[order(physeq_phylum_melt$Abundance, decreasing = TRUE),]

# Create a stacked barplot
barplot_all <- ggplot(physeq_phylum_melt, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme_minimal() +
  labs(title = "Relative Abundance at Phylum Level",
       x = "Sample",
       y = "Relative Abundance",
       fill = "Phylum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
barplot_all


plot_bar(per.f, x="Location", fill="Phylum")


# Barplot of 20 most abundant---------
# Top 20 most abundant
asv_totals <- taxa_sums(per.f)

# Get the names of the top 20 most abundant ASVs
top_20_asvs <- names(sort(asv_totals, decreasing = TRUE)[1:20])

# Subset the phyloseq object to include only the top 20 ASVs
physeq_top20 <- prune_taxa(top_20_asvs, per.f)

# Melt the phyloseq object to long format
physeq_melt <- psmelt(physeq_top20)

## Make a color palette----
# Get 8 colors from Dark2 and 5 colors from Paired
colors_dark2 <- brewer.pal(10, "Spectral")      # 8 colors from Dark2
colors_paired <- brewer.pal(4, "Paired") # First 5 colors from Paired

# Combine the colors into a 13-color palette
colors_14 <- c(colors_dark2, colors_paired)

## Family level----
p <- ggplot(physeq_melt, aes(x = Location, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(x = "Location of Sampling", y = "Relative Abundance (%)") +
  pretty.theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors_14) 
p

# Genus level----
p1 <- ggplot(physeq_melt, aes(x = Location, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(x = "Location of Sampling", y = "Relative Abundance (%)") +
  pretty.theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # scale_fill_manual(values = colors_14) 
p1

p2 <- ggplot(physeq_melt, aes(x = Location, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(x = "Location of Sampling", y = "Relative Abundance (%)") +
  pretty.theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# scale_fill_manual(values = colors_14) 
p2

# Alpha diversity----------
set.seed(123)  # For reproducibility
physeq_rare <- rarefy_even_depth(
  phy.f,
  sample.size = min(sample_sums(phy.f)),  # Minimum sample size
  rngseed = TRUE,                          # Random number seed
  replace = FALSE,                         # Rarefaction without replacement
  verbose = TRUE                           # Show information
)

alpha_rarefied <- plot_richness(physeq_rare, x="Location", measures=c("Shannon")) +
  geom_bar(aes(fill = Location)) +
  # geom_point(aes(fill = Location), shape = 24, color = "black", stroke = 1.5, size = 6, alpha = 0.9) +
  xlab(NULL) +
  ylab("Shannon Diversity") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")) +
  pretty.theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
alpha_rarefied

alpha_rarefied <- plot_richness(physeq_rare, x = "Location", measures = c("Shannon")) +
  geom_bar(aes(x = Location, y = value, fill = Location), stat = "identity") +
  xlab(NULL) +
  ylab("Shannon Diversity") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")) +
  pretty.theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
alpha_rarefied

# Exporting a table with ASV names and their relative abundances
# Export the table----
# This is to export the ASV table to include the taxonomy and ASV count
# Step 1: Extract ASV table and taxonomy table
otu_df <- as.data.frame(otu_table(per.f))         # Convert OTU table to a data frame
tax_df <- as.data.frame(tax_table(per.f))         # Convert taxonomy table to a data frame

# Step 2: Add ASV IDs to taxonomy table (rownames are the ASV IDs)
tax_df$ASV <- rownames(tax_df)

# Step 3: Merge ASV abundances with taxonomy
asv_tax_table <- cbind(otu_df, tax_df)

# Step 4: Save to CSV
write.csv(asv_tax_table, file = "ASVs_with_Taxonomy.csv", row.names = TRUE)



