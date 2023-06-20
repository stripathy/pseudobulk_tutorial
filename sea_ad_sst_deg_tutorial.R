library(tidyverse)
library(Seurat)
library(readxl)
library(edgeR)
library(ggrepel)
library(cowplot)

### Read in pseudobulked gene expression for SST cells and metadata
pseudobulk_sst_counts = read_csv('/external/rprshnas01/netdata_kcni/stlab/SEA-AD/SST_pseudobulk.csv')
sea_ad_meta = read_csv('/external/rprshnas01/netdata_kcni/stlab/SEA-AD/SST_pseudobulk_meta.csv')
gene_info = read_csv('/external/rprshnas01/netdata_kcni/stlab/SEA-AD/SST_pseudobulk_gene_info.csv')
rownames(pseudobulk_sst_counts) = gene_info$gene_symbol

sea_ad_meta = sea_ad_meta %>% mutate(disease = factor(disease, levels = c('normal', 'dementia')))

### First plot: Number of SST cells per person to illustrate that there are fewer SST cells in individuals with dementia than controls
sea_ad_meta %>% 
  ggplot(aes(x = disease, y = cells_per_donor)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  ylab('SST cells per donor')


### Create a cpm object with normalized gene expression counts for each pseudobulk sample
pseudobulk_sst_cpm = edgeR::cpm(y = pseudobulk_sst_counts)
rownames(pseudobulk_sst_cpm) = gene_info$feature_name

pseudobulk_sst_cpm_trans = pseudobulk_sst_cpm %>% 
  t() %>% 
  as.data.frame()

pseudobulk_sst_cpm_trans = left_join(sea_ad_meta, pseudobulk_sst_cpm_trans %>% rownames_to_column(var = 'donor_id')) 

# Plot gene expression of a gene (SST below) vs. all cell types
pseudobulk_sst_cpm_trans %>% 
  ggplot(aes(x = disease, y = log2(SST+1))) + 
  geom_boxplot() + 
  theme_cowplot() + 
  ylab('SST expr (log2 CPM+1)')

# Show how SST normalized expression is related to the number of cells sampled per donor
pseudobulk_sst_cpm_trans %>% 
  ggplot(aes(x = log2(SST+1), y = log10(cells_per_donor + 1), color = disease)) + 
  geom_point() + 
  theme_cowplot() + 
  xlab('SST cells per donor (log10 units)') + 
  ylab('SST expr (log2 CPM+1)')


#### Set up information for the stats model in limma-voom

# Find subjects with information for disease, Sex, PMI, and Age
use_subjects = complete.cases(sea_ad_meta %>% dplyr::select(disease, Sex, PMI, Age_norm))

# Drop subjects with fewer than 50 SST cells per donor
use_subjects = use_subjects & (sea_ad_meta$cells_per_donor > 50)

disease = factor(unlist(sea_ad_meta[use_subjects, 'disease']), levels = c('normal', 'dementia'))
sex = sea_ad_meta[use_subjects, 'Sex'] %>% unlist
pmi = sea_ad_meta[use_subjects, 'PMI'] %>% unlist
age = sea_ad_meta[use_subjects, 'Age_norm'] %>% unlist
cells_per_donor = sea_ad_meta[use_subjects, 'cells_per_donor'] %>% unlist

# Count the number of subjects for each disease group
sea_ad_meta[use_subjects, ] %>% 
  group_by(disease) %>% 
  tally

### Set up DGEList object
dge0 = DGEList(pseudobulk_sst_counts[, use_subjects], genes = gene_info)

min_samples_expressing_gene <- sea_ad_meta[use_subjects, ] %>% nrow * 0.8 # The 0.8 here refers to the fraction of total samples that needs to express the gene
dge0 = dge0[rowSums(dge0$counts >= 1) >= min_samples_expressing_gene, ] # This step filters genes such that they need to be detected in at least 80% of samples

# dge0 = DGEList(pb_counts, group = pb_metadata$Phenotype)
dge0 = calcNormFactors(dge0, method = "TMM")

### Set up design based on the factors defined above
design = model.matrix(~ age + pmi + sex + log10(cells_per_donor) + disease) 

# Perform voom transformation
vm = voom(dge0, design, plot = TRUE)

# Perform lmFit and eBayes
fit = lmFit(vm, design)
fit = eBayes(fit)

### Analyze sex and disease terms

# Analyze genes associated with sex
deg_table_sex = topTable(fit, coef = "sexMale",  n = Inf, sort = "none", adjust.method = "BH")

deg_table_sex %>% 
  arrange(adj.P.Val) %>% 
  head()

# Plot volcano plot of genes associated with sex
sex_volcano = deg_table_sex %>% 
  ggplot(aes(x = logFC, y = -log10(P.Value), label = feature_name)) + 
  geom_point() +
  geom_text_repel(data = subset(deg_table_sex, -log10(P.Value) > 10), 
                  aes(label = feature_name), 
                  vjust = 1.5) +
  theme_cowplot()

sex_volcano

# Analyze genes associated with disease
deg_table_disease = topTable(fit, coef = "diseasedementia",  n = Inf, sort = "none", adjust.method = "BH")

deg_table_disease %>% 
  arrange(adj.P.Val) %>% 
  head(20)

# Plot volcano plot of genes associated with disease
disease_volcano = deg_table_disease %>% 
  ggplot(aes(x = logFC, y = -log10(P.Value), label = feature_name)) + 
  geom_point() +
  geom_text_repel(data = subset(deg_table_disease, adj.P.Val < 0.1), 
                  aes(label = feature_name), 
                  vjust = 1.5) + 
  geom_point(data = subset(deg_table_disease, adj.P.Val < 0.1), 
             aes(color = 'red')) + 
  theme_cowplot()

disease_volcano

# Plot gene expression of a couple of DE genes below
ihih5_plot = pseudobulk_sst_cpm_trans %>% 
  ggplot(aes(x = disease, y = log2(ITIH5+1))) + 
  geom_boxplot() + 
  theme_cowplot()

ihih5_plot

aqp4_as1_plot = pseudobulk_sst_cpm_trans %>% 
  ggplot(aes(x = disease, y = log2(`AQP4-AS1`+1))) + 
  geom_boxplot() + 
  theme_cowplot()

aqp4_as1_plot