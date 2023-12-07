library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(readr)
library(tidyverse)
library(stringr)
library(gridExtra)
library(scales)

#options(scipen = 999)

setwd('/Users/zkoenig/Documents/hgdp_tgp_local/revisions/')

# TSV with n_snp, indels, fraction at 10/20x, geographic region, pop, and sample ID
rev_data <- read_tsv(
  "fig1_plot_data.tsv")

# # Reading in the mean number of SVs per population
# mean_sv <- read_csv("/Users/zkoenig/Documents/hgdp_tgp_local/revisions/svs_per_genome.csv")

# Reading in the metadata file which has plotting colors
metadata <- read_tsv('gnomad_meta_for_plots.tsv')

# Reading in the post QC summary info to plot cov/project
post_qc_info <- read.csv(
  '/Users/zkoenig/Documents/hgdp_tgp_local/postQC_summaries/post_qc_summary.tsv', sep = '\t')

# File which matches project to population
pop_to_proj <- read.csv('/Users/zkoenig/Documents/hgdp_tgp_local/postQC_summaries/pop_to_project.tsv', sep='\t')

# Renaming some column names for ease of plotting
names(rev_data)[names(rev_data) == 'hgdp_tgp_meta.Genetic.region'] <- 'genetic_region'
names(rev_data)[names(rev_data) == 'all'] <- 'n_svs'
names(rev_data)[names(rev_data) == 'CNV'] <- 'n_cnvs'

# Removing NA values from the dataset
# for the following two samples there are no number of CNVs or SVs so they were removed:
# NA12546B (CEU) and NA18874A (YRI)
rev_data <- na.omit(rev_data)


# Calculating the mean number of SNVs and indels per population
mean_all <- rev_data %>%
  group_by(population) %>%
  summarize(mean_snp = mean(n_snp), 
            mean_indels = mean(n_indels),
            mean_cnvs = mean(n_cnvs),
            mean_svs = mean(n_svs))

# Merging the mean data onto the other data
mean_data <- merge(x=mean_all, y=metadata, by="population")


# selecting the columns we want and then renaming
merge_select <- mean_data %>% select(genetic.region, population, mean_snp,
                                     mean_indels, mean_svs, mean_cnvs)

# Renaming some column names for ease of plotting
names(merge_select)[names(merge_select) == 'genetic.region'] <- 'genetic_region'
# Renaming some column names for ease of plotting
names(mean_data)[names(mean_data) == 'genetic.region'] <- 'genetic_region'

# Getting only one entry per population (should be 80)
update_data <- distinct(merge_select)




### Plotting the mean number of SNVs, indels, SVs per pop ###

# Vector which contains the order (by mean snv) for the geographical region for the legend
region_order <- c("AFR", "MID", "AMR", "EAS", "OCE", "CSA", "EUR")

# Creating a vector with regional colors from metadata
region_vec <- setNames(unique(metadata$continent.colors),
                       unique(metadata$genetic.region))


# Ordering the mean snvs in the df
mean_snv_sort <- update_data[with(update_data, order(-mean_snp)),]

# Changing geographical data to a factor so that the figure is sorted by geographic region
mean_snv_sort$population <- factor(mean_snv_sort$population,
                               levels = mean_snv_sort$population[order(mean_snv_sort$genetic_region)])

# Arranging data by descenting mean_n_snp
df_sort <- mean_snv_sort %>% arrange(desc(mean_snp))

# Setting the genetic region as a factor with the order as that in region_order
df_sort$genetic_region <- factor(df_sort$genetic_region, levels = region_order)
# Sorting the data based on the order in genetic region
df_sort <- df_sort[order(df_sort$genetic_region),]

# Setting the population as a factor based on genetic region so it shows in the plot
df_sort$population <- factor(df_sort$population,
                             levels = df_sort$population[order(df_sort$genetic_region)])



# Dividing the mean number of metrics by different values to result in better graphing
df_sort$plot_snp <- df_sort$mean_snp/10^6
df_sort$plot_sv <- df_sort$mean_svs/10^3
df_sort$plot_indel <- df_sort$mean_indels/10^5
df_sort$plot_cnv <- df_sort$mean_cnvs/10^2

# creating bar plot of the mean n_snv/individual per population
snv_plot <- ggplot(df_sort, aes(x=population, y=plot_snp, fill=genetic_region)) + 
  geom_col() + 
  theme_classic() +
  scale_fill_manual(name="Region", values = region_vec, breaks = region_order) +
  theme(axis.text.y = element_text(size=20, face="bold", colour = "black"),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(angle=90, hjust=1, size=18),
        axis.text = element_text(color="black")) +
  # ggtitle("Mean number of SNVs per Individual per Population") + 
  labs(x="Population", y=expression(atop(paste("Mean SNVs"),
                                         paste("(*10"^"6",")"))))

show(snv_plot)
# save_plot('snvs_per_pop.png', snv_plot, base_height=7, base_width=20)

# # Ordering the mean_svs in the df
# mean_sv_sort <- mean_data[with(mean_data, order(-mean_svs)),]
# 
# # Changing geographical data to a factor so that the figure is sorted by geographic region
# mean_sv_sort$population <- factor(mean_sv_sort$population,
#                                    levels = mean_sv_sort$population[order(mean_sv_sort$genetic_region)])

# creating bar plot of the mean n_sv/individual per population
sv_plot <- ggplot(df_sort, aes(x=population, y=plot_sv, fill=genetic_region)) + 
  geom_col() + 
  theme_classic() +
  scale_fill_manual(name="Region", values = region_vec, breaks = region_order) +
  theme(axis.text.y = element_text(size=20, face="bold", colour = "black"),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(angle=90, hjust=1, size=18),
        axis.text = element_text(color="black")) +
  # ggtitle("Mean number of SVs per Individual per Population") + 
  labs(x="Population", y=expression(atop(paste("Mean SVs"),
                                         paste("(*10"^"3",")"))))
show(sv_plot)
# save_plot('svs_per_pop.png', sv_plot, base_height=7, base_width=20)

# # Ordering the mean indels in the df
# mean_indel_sort <- mean_data[with(mean_data, order(-mean_indels)),]
# 
# # Changing geographical data to a factor so that the figure is sorted by geographic region
# mean_indel_sort$population <- factor(mean_indel_sort$population,
#                                levels = mean_indel_sort$population[order(mean_indel_sort$genetic_region)])

# creating bar plot of the mean n_sv/individual per population
indel_plot <- ggplot(df_sort, aes(x=population, y=plot_indel, fill=genetic_region)) + 
  geom_col() + 
  theme_classic() +
  scale_fill_manual(name="Region", values = region_vec, breaks = region_order) +
  theme(axis.text.y = element_text(size=20, face="bold", colour = "black"),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(angle=90, hjust=1, size=18),
        axis.text = element_text(color="black"),
        plot.margin = margin(0, 1, 0, 1, "cm")) +
  scale_y_continuous(trans='log10') +
  # ggtitle("Mean number of Indels per Individual per Population") + 
  labs(x="Population", y=expression(atop(paste("Mean Indels"),
                                         paste("(*10"^"5",")"))))

show(indel_plot)
# save_plot('indels_per_pop.png', indel_plot, base_height=7, base_width=20)


# creating bar plot of the mean n_cnv/per individual per population
cnv_plot <- ggplot(df_sort, aes(x=population, y=plot_cnv, fill=genetic_region)) + 
  geom_col() + 
  theme_classic() +
  scale_fill_manual(name="Region", values = region_vec, breaks = region_order) +
  theme(axis.text.y = element_text(size=20, face="bold", colour = "black"),
        axis.title.y = element_text(size=25),
        axis.text.x = element_text(angle=90, hjust=1, size=18),
        axis.text = element_text(color="black"),
        plot.margin = margin(0, 1, 0, 1, "cm")) +
  # ggtitle("Mean number of CNVs per Individual per Population") + 
  labs(x="Population", y=expression(atop(paste("Mean CNVs"),
                                         paste("(*10"^"2",")"))))

show(cnv_plot)
# save_plot('cnvs_per_pop.png', cnv_plot, base_height=7, base_width=20)

# 
# # Grabbing legend from one of the plots to insert in new plot
# legend <- get_legend(indel_plot)


# Creating a stacked plot of all 3 individual plots
# using the y labs from one plot, removing all other x axis info & legends
merge_plot <- plot_grid(snv_plot + theme(axis.title.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.text.x = element_blank(),
                                         legend.position = "none"),
                        indel_plot + theme(axis.title.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.text.x = element_blank(),
                                           legend.position = "none"),
                        sv_plot + theme(axis.title.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.text.x = element_blank(),
                                           legend.position = "none"),
                        cnv_plot  + theme(axis.title.x = element_blank(),
                                         legend.position = "none"),
                        nrow = 4, 
                        rel_heights = c(1,1,1,1.65),
                        align = "v"
                        )
merge_plot

# # Merging the combined 3 plots with the legend, setting width of plots relative to legend
# final_plot <- plot_grid(merge_plot, legend, rel_widths = c(10,1))
# 
# final_plot
save_plot('indel_sv_snv_cnv_stack.png', merge_plot, base_height=15, base_width=24)
save_plot('indel_sv_snv_cnv_stack.pdf', merge_plot, base_height=15, base_width=24)



# ### Creating the 10x and 20x coverage plots ###
# 
# # Gather 10x and 20x coverage fractions into one column
# rev_pivot <- rev_data %>% pivot_longer(c(pct_bases_20x, pct_bases_10x),
#                                        names_to = "coverage",
#                                        values_to = "fraction_of_bases")
# 
# 
# # Creating df with proj mapped to pop
# pop_to_proj <- update_merge %>% select(population.x, project)
# names(pop_to_proj)[names(pop_to_proj) == 'population.x'] <- 'population'
# pop_to_proj <- unique(pop_to_proj)
# 
# # Merging the project ID onto the downsampled dataset
# post_qc_info <- merge(post_qc_info, pop_to_proj, by='population')
# 
# 
# # Plot 10x and 20x coverage with density plot
# p_1020_density <- ggplot(rev_pivot, aes(x=fraction_of_bases, 
#                                         fill=coverage,
#                                         color=coverage,
#                                         group=coverage)) + 
#   geom_density(adjust=1.5, alpha=.3, size=1) +
#   theme_bw() +
#   theme(text = element_text(size=20),
#         axis.text.y = element_text(size=23, colour = "black"),
#         axis.title.y = element_text(size=25),
#         axis.text = element_text(color="black")) +
#   labs(x="Percentage of Bases", y="Density", fill="Coverage", color="Coverage") +
#   guides(fill = guide_legend(), color = guide_legend()) + 
#   ggtitle("Percentage of Bases at \n10 and 20x Coverage") +
#   scale_fill_discrete(labels=c('10x', '20x')) +
#   scale_color_discrete(labels=c("10x", "20x"))
# 
# show(p_1020_density)
# #save_plot('p_1020_density.png', p_1020_density, base_height=10, base_width=13)
# 
# 
# 
#   
# # Plotting the mean coverage per individual comparing HGDP and 1kG
# proj_cov <-  ggplot(post_qc_info, 
#                           aes(x=cov_stats.mean, 
#                               group=project, 
#                               color=project, 
#                               fill=project)) +
#   geom_density(adjust=1.5, alpha=.3, size=1) +
#   theme_bw() +
#   theme(text = element_text(size=20),
#         axis.text.y = element_text(size=23, colour = "black"),
#         axis.title.y = element_text(size=25),
#         axis.text = element_text(color="black")) +
#   ggtitle("Mean Coverage of 1kGP vs HGDP") +
#   labs(x="Mean Coverage per Individual", y="Density", fill="Project",
#        color="Project") +
#   scale_fill_discrete(labels=c('1kGP', 'HGDP')) +
#   scale_color_discrete(labels=c("1kGP", "HGDP"))
# 
# show(proj_cov)
# 
# # y_axis_title <- get_y_axis(proj_cov)
# 
# merge_cov_plot <- plot_grid(proj_cov,# + 
#                             # theme(axis.title.x = element_blank(),
#                             #       axis.ticks.x = element_blank(),
#                             #       axis.text.x = element_blank(),
#                             #       # legend.position = "none",
#                             #       title = element_blank()),
#                             p_1020_density,# +
#                             # theme(axis.title.x = element_blank(),
#                             #       axis.ticks.x = element_blank(),
#                             #       axis.text.x = element_blank()),
#                             #       legend.position = "none"),
#                             # y_axis_title,
#                             ncol = 2, 
#                             rel_heights = c(1,1),
#                             align = "h",
#                             labels="AUTO",
#                             label_size = 18
#                             )
# merge_cov_plot
# 
# 
# save_plot('cov_proj_1020x.png', merge_cov_plot, base_height=7, base_width=14)

