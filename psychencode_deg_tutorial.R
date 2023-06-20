library(edgeR)
library(tidyverse)
library(cowplot)
library(ggrepel)

# Load pseudobulk expression counts from PsychEncode
pseudo_fine = read_tsv('/external/rprshnas01/netdata_kcni/stlab/SZBDMulticohort/pseudobulk_fine_grained.counts.tsv')
pseudo_fine[1:5, 1:10] # These are counts!

# Load pseudobulk metadata from PsychEncode
pseudo_meta = read_tsv('/external/rprshnas01/netdata_kcni/stlab/SZBDMulticohort/pseudobulk_fine_grained.metadata.tsv')

# Create a new object, pseudo_fine_cpm, with the pseudobulk counts normalized as counts per million
pseudo_fine_cpm = edgeR::cpm(y = t(pseudo_fine[, 3:17660]))

# Add the metadata about cell types and individuals back to the pseudobulk counts
pseudo_fine_cpm_trans = cbind(pseudo_fine[, 1:2,], t(pseudo_fine_cpm))

# Merge pseudobulk metadata and cpm data frame
pseudo_fine_cpm_trans = left_join(pseudo_meta, pseudo_fine_cpm_trans)

# Make a plot that shows gene expression of a gene (SST below) vs all cell types
pseudo_fine_cpm_trans %>% 
  ggplot(aes(x = Celltype, y = log2(SST+1), fill = Phenotype)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


#### Set up information that we'll use for our stats model in limma-voom

# Set things up to use a pseudobulk matrix from just SST cells from the McLean cohort
USE_CELL_TYPE = 'In-PV_Basket'
USE_COHORT = 'McLean'

# Identify subjects with information for disease, sex, PMI, and age
use_subjects = pseudo_meta %>% 
  filter(Celltype == USE_CELL_TYPE, Cohort == USE_COHORT, num_cells > 50) %>% 
  pull(unique_donor_ID)

disease = factor(pseudo_meta %>% 
                   filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>% 
                   pull(Phenotype), levels = c('CON', 'SZ'))
sex = pseudo_meta %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>% 
  pull(Gender)
pmi = pseudo_meta %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>% 
  pull(PMI)
age = pseudo_meta %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>% 
  pull(Age)
cells_per_donor = pseudo_meta %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>% 
  pull(num_cells)

# Count the number of subjects for each disease group
pseudo_meta %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>% 
  group_by(Phenotype) %>% 
  tally

# Extract the pseudobulk counts for SST cells of the selected subjects
pseudobulk_sst_counts = pseudo_fine %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE)

pseudobulk_sst_counts = pseudobulk_sst_counts[, 3:17660] %>% t()
rownames(pseudobulk_sst_counts) = colnames(pseudo_fine[, 3:17660])
colnames(pseudobulk_sst_counts) = pseudo_meta %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>% 
  pull(unique_donor_ID)

### Set up DGEList object 
dge0 = DGEList(pseudobulk_sst_counts, genes = row.names(pseudo_fine_cpm))

min_samples_expressing_gene <- pseudo_meta %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>% 
  nrow * 0.8 # The 0.8 here refers to the fraction of total samples that needs to express the gene
dge0 = dge0[rowSums(dge0$counts >= 1) >= min_samples_expressing_gene, ] # Filter genes that are detected in at least 80% of samples

dge0 = calcNormFactors(dge0, method = "TMM")

### Set up design based on the factors defined above
design = model.matrix(~ age + pmi + sex + log10(cells_per_donor) + disease) 

# Perform voom transformation
vm = voom(dge0, design, plot = TRUE)

# Perform lmFit and eBayes steps
fit = lmFit(vm, design)
fit = eBayes(fit)

### Plot genes associated with sex
deg_table_sex = topTable(fit, coef = "sexMale",  n = Inf, sort = "none", adjust.method = "BH")

deg_table_sex %>% arrange(adj.P.Val) %>% head()

# Plot volcano plot of genes associated with sex
sex_volcano = deg_table_sex %>% 
  ggplot(aes(x = logFC, y = -log10(P.Value), label = genes)) + 
  geom_point() +
  geom_text_repel(data = subset(deg_table_sex, -log10(P.Value) > 10), 
                  aes(label = genes), 
                  vjust = 1.5) +
  theme_cowplot()

sex_volcano

### Plot genes associated with disease
deg_table_disease = topTable(fit, coef = "diseaseSZ",  n = Inf, sort = "none", adjust.method = "BH")

deg_table_disease %>% arrange(adj.P.Val) %>% head

# Plot volcano plot of genes associated with disease
disease_volcano = deg_table_disease %>% 
  ggplot(aes(x = logFC, y = -log10(P.Value), label = genes)) + 
  geom_point() +
  geom_text_repel(data = subset(deg_table_disease, adj.P.Val < 0.1), 
                  aes(label = genes), 
                  vjust = 1.5) + 
  geom_point(data = subset(deg_table_disease, adj.P.Val < 0.1), 
             aes(color = 'red')) + 
  theme_cowplot()

disease_volcano

# Make plots that show gene expression of a couple of DE genes
arhgap18_plot = pseudo_fine_cpm_trans %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>%
  ggplot(aes(x = Phenotype, y = log2(ARHGAP18+1))) + 
  geom_boxplot() + 
  theme_cowplot()

arhgap18_plot

znf641_plot = pseudo_fine_cpm_trans %>% 
  filter(unique_donor_ID %in% use_subjects, Celltype == USE_CELL_TYPE) %>%
  ggplot(aes(x = Phenotype, y = log2(ZNF641+1))) + 
  geom_boxplot() + 
  theme_cowplot()

znf641_plot