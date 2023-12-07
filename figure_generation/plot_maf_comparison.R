library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(cowplot)
library(grid)
library(gridExtra)
library(patchwork)
library(scales)

# setting the working directory
setwd('/Users/zkoenig/Documents/hgdp_tgp_local/dataset_comparison/comparison_tables')


# Creating a column with the difference between n_var (total var in HGDP+1kGP)
getVarDiff <- function(df) {
  # Create a column with the difference between n_var (total var in HGDP+1kGP)
  # and n_var_in_both to get the n_var for HGDP+1kGP only
  df$hgdp_1kg_only <- df$n_var - df$n_var_in_both
  return(df)
}

pivotLonger <- function(df) {
  # pivot longer the dataset
  pivot_df <- df %>%
    pivot_longer(c(n_var, hgdp_1kg_only, n_var_in_both),
                 names_to = "comparison", values_to="num_var")
  
  pivot_df<- pivot_df[-c(2),]
  
  return(pivot_df)
}

updateFieldNames <- function(pivot_df, dataset_name) {
  # Changing the field names to those which will be used in plotting
  # also creating a 'in comparison only' entry for when n_var_in_both = 0
  pivot_df <- pivot_df %>% mutate(
    comparison = replace(comparison, (comparison == "n_var" & cat == "0%"),
                         "In comparison dataset only"),
    comparison = replace(comparison, comparison == "hgdp_1kg_only", 
                         "In HGDP+1kGP only"),
    comparison = replace(comparison, (comparison == "n_var_in_both" & num_var != 0),
                         "In both HGDP+1kGP and comparison dataset")
  ) %>% 
    as.data.frame()
  
  # removing the n_var rows from the dataset (those are the counts for the total vars in HGDP+1kGP)
  pivot_df <- subset(pivot_df, (comparison !="n_var" & comparison != "n_var_in_both"))
  
  # adding a column with the dataset name so it can be labelled for plotting
  pivot_df$dataset <- dataset_name
  
  return(pivot_df)
}


# reading in annotated comparison datasets (exported from hail)
# they already contain binning information for the histogram plotting
df_berg <- read_tsv('berg_new_hist.tsv')
df_gnomad <- read_tsv('gnomad_new_hist.tsv')
df_nygc <- read_tsv('nygc_new_hist.tsv')
df_phase3 <- read_tsv('phase3_new_hist.tsv')


# adding a column which has the number of vars in HGDP+1kGP only
df_berg <- getVarDiff(df_berg)
df_gnomad <- getVarDiff(df_gnomad)
df_nygc <- getVarDiff(df_nygc)
df_phase3 <- getVarDiff(df_phase3)

# pivoting the datasets longer
pivot_berg <- pivotLonger(df_berg)
pivot_gnomad <- pivotLonger(df_gnomad)
pivot_nygc <- pivotLonger(df_nygc)
pivot_phase3 <- pivotLonger(df_phase3)

# updating the field names for plotting, adding dataset col
pivot_berg <- updateFieldNames(pivot_berg, "Bergstrom HGDP")
pivot_gnomad <- updateFieldNames(pivot_gnomad, "gnomAD")
pivot_nygc <- updateFieldNames(pivot_nygc, "NYGC 1kGP")
pivot_phase3 <- updateFieldNames(pivot_phase3, "Phase3 1kGP")


# merging all of the separate long form datasets into one large dataframe
merge_df <- rbind(pivot_berg, pivot_nygc)
merge_df <- rbind(merge_df, pivot_gnomad)
merge_df <- rbind(merge_df, pivot_phase3)

# this does not work, but do need to get the order correct so that it
# shows up correctly in the key
merge_df$cat <- factor(merge_df$cat, levels = unique(merge_df$cat))


# going to write out a table with the counts for each dataset to add 
# to the manuscript as a supplementary figure for figure 4 (which is this plot)
comp_figure <- merge_df[, c("cat", "comparison", "num_var", "dataset")]


# Setting the comparison as a factor so that it appears as desired in the legend
merge_df$comparison <- factor(merge_df$comparison, 
                              levels = 
                                c("In comparison dataset only",
                                  "In HGDP+1kGP only",
                                  "In both HGDP+1kGP and comparison dataset"
                                ))

# using facet grid to get the multiple plots in a grid format
p <- ggplot(merge_df, 
            aes(y=num_var, 
                x=cat, 
                fill=comparison)) + 
  theme_bw() +
  theme(text = element_text(size=18),
        legend.title = element_blank(),
        legend.position = "bottom") +
  geom_bar(position=position_dodge2(reverse=TRUE, padding = 0), stat="identity") +
  xlab("Minor Allele Frequency in HGDP+1kGP Dataset") +
  ylab("Number of Variants (Note Log Scale)") +
  scale_fill_manual(
    values = c("#377EB8", "#984EA3", "#E41A1C"),
    breaks = c("In comparison dataset only",
               "In both HGDP+1kGP and comparison dataset",
               "In HGDP+1kGP only")) +
  scale_y_log10(label=comma) +
  coord_cartesian(ylim=c(100,NA)) +
  facet_wrap(~dataset)

show(p)

# # Writing out plots as png and pdf
save_plot('figure4_dataset_comparisons.png', p, base_height=8, base_width=14)
save_plot('figure4_dataset_comparisons.pdf', p, base_height=8, base_width=14)


# # Trying a different color scheme
# # using facet grid to get the multiple plots in a grid format
# # Hex for coloring: 
# # red - #DC267F
# # purple - #785EF0
# # blue - #648FFF
# p1 <- ggplot(merge_df, 
#              aes(y=num_var, 
#                  x=cat, 
#                  fill=comparison)) + 
#   theme_bw() +
#   theme(text = element_text(size=12),
#         legend.title = element_blank()) +
#   geom_bar(position=position_dodge2(reverse=TRUE), stat="identity") +
#   ggtitle("Minor Allele Frequency Comparisons") +
#   xlab("Minor Allele Frequency in HGDP+1kGP Dataset") +
#   ylab("Number of Variants") +
#   scale_fill_manual(
#     #labels = c(
#     #"In comparison dataset only",
#     #"In both HGDP+1kGP and\n comparison dataset",
#     #"In HGDP+1kGP only"),
#     values = c("#648FFF", "#785EF0","#DC267F")) +
#   facet_wrap(~dataset)
# 
# show(p1)

# Writing out plots as png and pdf
# save_plot('maf_dataset_comparisons_nolog.png', p1, base_height=10, base_width=13) 
# save_plot('maf_dataset_comparisons_nolog.pdf', p1, base_height=10, base_width=13) 
