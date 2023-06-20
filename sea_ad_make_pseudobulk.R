library(tidyverse)
library(Seurat)
library(readxl)
library(edgeR)
library(ggrepel)
library(cowplot)

# ad <- anndata::read_h5ad('/external/rprshnas01/netdata_kcni/stlab/SEA-AD/Sst_Chodl.h5ad')

### this code loads in the data from cell by gene- run this from the terminal
curl -o sst_chodl.rds "https://corpora-data-prod.s3.amazonaws.com/7365e4e1-4711-494b-8171-fda7f5611872/local.rds?AWSAccessKeyId=ASIATLYQ5N5XYBGLKSE3&Signature=FIe8AC%2FDlyTWGp27E3nA5Rdkbk4%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEOj%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJGMEQCIEKJvhtBe%2F48Pfy1LBAF%2B0ScY7jpPNbxsDODkT4QbvvrAiBT%2FMVEIC3e%2B9IXLNmxb6JNdX1sXyKfVWkAwwpgNqoDySrrAwhBEAEaDDIzMTQyNjg0NjU3NSIM4fR6G2kFcno%2F3KtGKsgDo4SuTWcjQbMH0sOHUnW%2B8jyh6cvzdxYpfmk0bwzrpxkc1eHmWQ1GLUgbWQTmSTFPwLBscByjxeqeczUe24aiEQ8%2FfPuUY6ZLhvGabt0E2hTnB6N69lQmZeW2lzLvf0lUf7eaMwQnEytMxnRtljHBJ0run%2Btslkz4Db2l%2B9oTrw5YAkZsiapaumUDBwID7Iw4zhuC%2BSYgSEEg1Y5nwMgplw0QS8fM1lkd%2F3DwFZa32dbRansBt0S8VQKcFilsAFNBtQvVOZA8ZAyt2TIUfaEtt02Aew7YwyrAY8FLWE9DazT5QIku7SjZjCA2gyNelMa%2BaNyHFDiY1Dt8wIBw38rXkhf6sbov4ztZhLElec6bVlOn833L0MoSv1heqlgNrvQ%2BIBBQq%2F6HqCSYlvaaysQGNn43tvhr6uZ84HwT3Ax8fg0tB31BYd2ydNlz%2BU5Kz6YrKL%2BTHY0NicJfWyI2w5YPRQnF3UylDv%2FGlwdfE0r2OMQbmZMUIsXNQPjZxtT3s9cEDlS%2Fum6mDum6r%2BJ9iNY1uvgqLY8bZ8PK%2BM4V%2FqvfiXB0uJJ%2FwaqaVmmYjRY5VWzGQdlyEdX12TX0iA9UrkaKvDd4DAFgOXDkMOGzsKQGOqYBPyST1ZQRVNeaQt85fS8RNl4sek4wdm0Zy0OZXDA5fcO6FJyCxyxAz51Aw6TssSYwHVKeVz55kuaN2uA7XL4kTrOlk6j6m8fR9fku8fTlbwOjsoB9fAwZH4OTw9QC%2FQOMiJGlcovkMC9MKLOjbeooC6IUFpt%2Flr9a%2Fos5lkv9ChAATGzNDMzj96cc5gcuvuXRTrfc1eQXA6eIZUJex8LPIdCLXEELCA%3D%3D&Expires=1687526694"

curl -o sst.rds "https://corpora-data-prod.s3.amazonaws.com/673ff65f-a98b-49f3-a4a9-175f93a36860/local.rds?AWSAccessKeyId=ASIATLYQ5N5XT3F7DM6G&Signature=cGmztMcVPJTLY9FvKfWy831XpkE%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEEoaCXVzLXdlc3QtMiJIMEYCIQCD%2FW29P02USD7aX%2BRjWRNd27uqWWqrAhZIJu3doDCGzAIhAOIGMezdim09WJOrz6gx%2F9IGppRpFj%2FcVrwN3daxhFfAKvQDCKP%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQARoMMjMxNDI2ODQ2NTc1IgwRqX3LoPRTd1SLMyYqyANdn0S7EYUfRiQW1%2BvoWhewtEfsN0OFxjHM7vZQVQnSVyngSMsR2YDin1D5Q7%2BbOFszFB9Y6yHJZ6z2yDGUeOsrA1S4YDtt75aIjhkukH2%2BGnTKvaYTrf0QaXAguF0e%2FVC12iIB3He%2FzMxUF2ywT9qnfOjzMQQomZcxRq691zNXe4LLO0vCqucRh3abRCHOBn0eOiz4Zmq1qndeXKDa8NwxSPGGVovHtbpzl86eubWqqTEEoxifMuov37%2Buno2Sl2cZRfpNDaOVrl6fhadcRc16Fj2EKT4qw%2F1P8e7lSBcHpOK5T62129VardcPB0c779vu8JP2kx%2FocJiUEwj0ed%2FPDqoozpHNiMzLqxf37oBkiW%2F1CzjLi7EfmBmK6ewsqkx84D%2F%2BlmLfNC7L15ECDp%2Bwh%2BPA2RjwCBaUodlNW89B4ayaLVYgqoAt7ynbZfZZsIvzWSn802EnWC%2FTrN16Igi%2FGK0UZOby8CcLT8mUecCHFYhzt7TZnYL36Ny3axhvY9nbMp%2FkNpWjQj3HOUa%2BeOJ6jwUBaXMzJBfTWOnNXIuIsWmRjMCgtqlEZwA6j1criPTWX7wTFGrZ22DAWJPyDU%2F2dReVdGnsYYowjOzFpAY6pAE3FVURQkNA1hHXe1P0JAzlh%2B4xiwRp8kLU6d4KSBKAyQYvCgzBgB%2FcDcRxiTIaaIvS0uTdGBniKk0Xnc5TN%2BqBJmVZUBWRFaXMj%2FVne4FXNf7J43FNszDMe6%2B3JrEd2AJJ9LqhrjV0nKWGGwCstt99%2BIDQRZe6FT%2FZY5vCeFQCJLgysKBBtV4XnTolxvGP4DXecOzQRICnmt3N23qg2%2FCNpmdv6A%3D%3D&Expires=1687871143"

curl -o sst_chodl.rds "https://corpora-data-prod.s3.amazonaws.com/7365e4e1-4711-494b-8171-fda7f5611872/local.rds?AWSAccessKeyId=ASIATLYQ5N5XVSUMPKPK&Signature=uPfBLAEQYTgAeoQJmRXHcqYnmtg%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEEkaCXVzLXdlc3QtMiJHMEUCIQCPzQFJDfK5AulwtmBjc1fQvivMOdUxnMUz4RG0VrpElgIgSkjJzR9wq%2BWSeH6D94RxuEm7mK7f2%2BgCjtMA%2Frj4cs8q9AMIov%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARABGgwyMzE0MjY4NDY1NzUiDCFX0YFkt%2FjrsELaMirIA9p7GbnrvCYMfLxzZ0c90crGA3V1JpDR%2BaUGAwL52U4Gdl6hvRf8%2BxxTqwOpk%2BXZSP8KONss%2FedJjKGUSjpY%2FzxeveogPvvT9tifetiALUz10VGgmX6cIFVawCp54LbtiTvPdqRNFPFfJRkJufb%2FFQ9GaI1MDHG%2Bs9J26zSasXGF0bi5QZKvcAlIwOUMcpd6bqbg4Ni0ElICV8kS%2FkPP17GWql5CwPs7QywJW1%2Fg7xfzScm3PXyQpWxYp0EFggEMKSxxzWb9AbUjuq2f%2FV8SnIq46GG%2FMRHrkwwAuWhiv7EIlc6F9x%2FtuLM7KyF3VQ%2BSTv%2F1So2O%2B5C1P05Ejx%2FH73HjHzXirrgDEIbHMoaPSX6HLjS1QcBgBsSuCCNv9BiQdzTNWu93A%2FuJC0TbIqEOfs4NHhLYpPvU9R%2B79CwTjJsA0qqtifpd7V8qdnq3N9i6vrUQZ6qV2oOqrAlVJjxD%2BqnsVjzYNDvhh4YHx3lPa%2F%2F%2BKxhYrB7705JTgBm24HSYOW9g2iw6WrI9U5i3Na0FKa2XH2lUYzvb1ohgB2KMIUnlFkE0LOUB%2FEmMjXsOP4IZMO0eUwrxg3QRB0LQSeZbe%2BN%2Fu2RKImup%2BTDo2MWkBjqlAQ9JSQ6EstVHfat%2BHuPtxracWiGQF01kca3nA%2BrYUWhuFzag8jyjPCOHjv5z5dJ0VkwvHVOItyDjuzALlVHmvn%2ByamRINGUyZzmtCzrT2aoXIGX5SjQwpg5C9XHq0hHs%2Fg4IXfqiVGhu2yPn8BHzcjMtFSMYUQaekQI5eCvY0JLgRLUmVOCosaYOeVxHhtbPeXLRaVZ%2B0eKVnOFHIQyDz5FyFo%2FWEw%3D%3D&Expires=1687871106"

### read in the downloaded Seurat objects and then pseudobulk them

# load in the single cell dataset from the SCC into Seurat
sst = readRDS('/external/rprshnas01/netdata_kcni/stlab/SEA-AD/sst.rds')

# read in donor metadata
sea_ad_donor_meta = read_excel('/external/rprshnas01/netdata_kcni/stlab/SEA-AD/sea-ad_cohort_donor_metadata.xlsx')

# calculate the number of cells per donor
cells_per_donor_df = sst@meta.data %>% group_by(donor_id) %>% tally(name = 'cells_per_donor')

# generate the pseudobulk of the number of cells per donor
pseudobulk_sst <- Seurat:::PseudobulkExpression(object = sst, group.by = 'donor_id', pb.method = 'aggregate', slot = 'counts')

# this turns the seurat object into a matrix
pseudobulk_sst_counts = pseudobulk_sst[['RNA']]

# this pulls the metadatafrom the seurat object
sea_ad_meta = sst@meta.data %>% distinct(donor_id
                                               , .keep_all = T)
rownames(sea_ad_meta) = sea_ad_meta$donor_id

# this generates a combined metadata object from both the seurat object and donor metadata above
sea_ad_meta = left_join(sea_ad_meta %>% select(donor_id, disease) %>% rownames_to_column(var = 'Donor ID'), sea_ad_donor_meta)

# add in information about numbers of cells from the donor
sea_ad_meta = left_join(sea_ad_meta, cells_per_donor_df)
rownames(sea_ad_meta) = sea_ad_meta$donor_id


# do some munging of subject ages and create a new column, Age_norm
sea_ad_meta = sea_ad_meta %>% mutate(Age_norm = case_when(`Age at Death` == '90+' ~ 90, 
                                                    TRUE ~ as.numeric(`Age at Death`)))

# reorder the counts to match the ordering of subjects in the metadata
pseudobulk_sst_counts = pseudobulk_sst_counts[, rownames(sea_ad_meta)] %>% as.data.frame()

gene_info = sst@assays$RNA@meta.features %>% rownames_to_column(var = 'gene_symbol')


write_csv(x = pseudobulk_sst_counts, '/external/rprshnas01/netdata_kcni/stlab/SEA-AD/SST_pseudobulk.csv')
write_csv(x = sea_ad_meta, '/external/rprshnas01/netdata_kcni/stlab/SEA-AD/SST_pseudobulk_meta.csv')
write_csv(x = gene_info, '/external/rprshnas01/netdata_kcni/stlab/SEA-AD/SST_pseudobulk_gene_info.csv')

                  