library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

setwd('/Users/zkoenig/Documents/hgdp_tgp_local/table_generation')

# Reading in the datasets
# table 1 with 4097 samples (11/16/21)
info <- read.csv('table_x.csv')
# reading in the downsamples dataset for singleton plots
info_downsample <- read.csv('/Users/zkoenig/Documents/hgdp_tgp_local/Singleton_Investigation/hgdp_tgp_bergstrom_metadata_downsample.tsv', sep='\t')
# file with hgdp_tgp metadata
metadata <- read.csv('gnomad_meta_hgdp_tgp_v1.txt', sep='\t')

# Renaming the region SAS to CSA (named differently between 1kG and HGDP)
info$Geographical_region <- sub("SAS", "CSA", info$Geographical_region)

# Changing geographical data to a factor so that the figure is sorted by geographic region
#info$Population <- factor(info$Population, levels = info$Population[order(info$Geographical_region)])

# Fixing the geographic region names from the downsampled dataset
info_downsample <- info_downsample %>% mutate(
  region = recode(info_downsample$region, "Africa" = "AFR",
                  "America" = "AMR", "Central_South_Asia" = "CSA", "East_Asia" = "EAS", 
                  "Europe" = "EUR","Middle_East" = "MID", "Oceania" = "OCE", "SAS" = "CSA"))

# Writing out the dataset for Mary with the correct information
mary_tsv <- read.csv('table_x.csv')
mary_tsv$Geographical_region <- sub("SAS", "CSA", mary_tsv$Geographical_region)

write.csv(mary_tsv, file="table_1.csv")


# Changing geographical data to a factor so that the figure is sorted by geographic region
#info_downsample$Population <- factor(info_downsample$Population, levels = unique(info_downsample$Population[order(info_downsample$region)]))

# Vector which contains the order for the grographical region for the legend
region_order <- c("AFR", "AMR", "CSA", "EAS", "EUR", "MID", "OCE")

# Creating a vector with regional colors from metadata
region_vec <- setNames(unique(metadata$hgdp_tgp_meta.Continent.colors), 
                       unique(metadata$hgdp_tgp_meta.Genetic.region))
# Creating a vector with color per poplation from metadata
pop_vec <- setNames(unique(metadata$hgdp_tgp_meta.Pop.colors), 
                    unique(metadata$hgdp_tgp_meta.Population))

# aggregating region by the mean n_snp for plotting
info_sort <- info %>% 
  group_by(Geographical_region) %>%
  mutate(mean_regional_nsnp = mean(n_snp_stats.mean))

# Ordering the df and applying sort order by setting pop as a factor
info_sort <- info_sort[with(info_sort, order(-mean_regional_nsnp, -n_snp_stats.mean)),]
info_sort$Population <- factor(info_sort$Population, levels = unique(info_sort$Population))

# aggregating region by the mean n_snp for plotting
info_downsample <- info_downsample %>% 
  group_by(region) %>%
  mutate(mean_regional_singleton = mean(n_singleton))


# Ordering the df and applying sort order by setting pop as a factor
info_downsample <- info_downsample[with(info_downsample, order(-mean_regional_singleton, -n_singleton)),]
info_downsample$Population <- factor(info_downsample$Population, levels = unique(info_downsample$Population))

# creating bar plot of the mean n_snp/individual per population
p <- ggplot(info_sort, aes(x = Population, y=n_snp_stats.mean, fill=Geographical_region)) + 
  geom_col() + 
  theme_bw() +
  scale_fill_manual(name="Region", values = region_vec, breaks = region_order) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.text = element_text(color="black")) +
  ggtitle("Mean number of SNPs per Individual per Population") + 
  labs(x="Population", y="Mean SNPs per Individual")

show(p)
save_plot('v1_snps_per_indiv.png', p, base_height=7, base_width=20) 


# plotting the number of singletons/individual per population
p1 <- ggplot(info_sort, aes(x = Population, y=n_singleton, fill=region)) + 
  geom_col() + 
  theme_bw() +
  scale_fill_manual(name="Region", values = region_vec, breaks = region_order) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.text = element_text(color="black")) +
  ggtitle("Mean number of Singletons per Individual per Population") + 
  labs(x="Population", y="Mean Singletons per Individual")

show(p1)
save_plot('singletons_per_indiv.png', p1, base_height=7, base_width=20) 

# plotting the mean coverage per individual comparing HGDP and 1kG
p2 <- ggplot(info, aes(x=cov_stats.mean, group=Project, color=Project, fill=Project)) + 
  geom_density(adjust=1.5, alpha=.3, size=1) + 
  theme_bw() +
  theme(text = element_text(size=20)) +
  ggtitle("Mean Coverage of 1kGP vs HGDP") + 
  labs(x="Mean Coverage per Individual", y="Density")
  
show(p2)
save_plot('v1_cov_1kg_hgdp.png', p2, base_height=7, base_width=14) 

# Dropping the OCE population for the coverage plots since there are only two subpopulations for that global population
info_oce <- info[info$Geographical_region !="OCE",]

# plotting the mean coverage per individual comparing geographical regions
p3 <- ggplot(info_oce, aes(x=cov_stats.mean, group=Geographical_region, 
                       color=Geographical_region, fill=Geographical_region)) + 
  geom_density(adjust=1.5, alpha=.3, size=1) + 
  theme_bw() +
  scale_color_manual(name="Geographical Region", values = region_vec, breaks = region_order) +
  scale_fill_manual(name="Geographical Region", values = region_vec, breaks = region_order) +
  theme(text = element_text(size=20),
        axis.text = element_text(color="black")) +
  ggtitle("Mean Coverage by Continental Region") + 
  labs(x="Mean Coverage per Individual", y="Density")

show(p3)
save_plot('v1_cov_region.png', p3, base_height=7, base_width=14) 
