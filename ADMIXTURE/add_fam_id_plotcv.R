library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(cowplot)

setwd('/Users/zkoenig/Documents/hgdp_tgp_local/ADMIXTURE')

# reading in hgdp_1kGP metadata file which has metrics needed for annotating
metadata <- read_tsv('gnomad_meta_hgdp_tgp_v1.txt')

fam_file <- read_tsv('hgdp_1kg_plink_data/post_qc_ld_unrel_only.fam',
                         col_names=FALSE)

global_region_id <- metadata %>% select('s', 'hgdp_tgp_meta.Genetic.region')

# Removing the v3.1 from gnomAD samples in the id file so there are no NA
global_region_id$s <- gsub("v3.1::", "", global_region_id$s )

merge_fam <- left_join(x = fam_file, y = global_region_id, by=c('X2'='s'))

merge_fam <- merge_fam[,c(7, 2, 3, 4, 5, 6 )]


ind2pop <- merge_fam[1]

# Getting a table with all of the continent colors for plotting
plot_colors <- metadata %>% select('hgdp_tgp_meta.Continent.colors',
                                   'hgdp_tgp_meta.Genetic.region')

plot_colors <- unique(plot_colors)

write.table(ind2pop, file = 'pong_files/ind2pop.txt', col.names = FALSE,
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(merge_fam, 
            file = 'hgdp_1kg_plink_data/new_post_qc_ld_unrel_only.fam',
            col.names = FALSE,
            row.names = FALSE, sep = "\t", quote = FALSE)


# plotting the cross validation
cv_vals <- read_tsv('cv_error.txt', col_names=FALSE)

cv_vals$X2 <- factor(cv_vals$X2, levels = cv_vals$X2)

p <- ggplot(data=cv_vals, aes(x=X2, y=X3, group=1)) +
  geom_line() +
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle=90, hjust=1),
        axis.text = element_text(color="black")) +
  labs(x='K', y='5-Fold Cross-Validation Error')

show(p)
save_plot('5-fold_cv.png', p, base_height=8, base_width=13)
save_plot('5-fold_cv.pdf', p, base_height=8, base_width=13)



###############################################################################

tmp_fam = fam_file[2]
tmp_meta_s = global_region_id['s']


write.table(tmp_fam, file = 'tmp_fam.txt', col.names = FALSE,
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(tmp_meta_s, file = 'tmp_meta_s.txt', col.names = FALSE,
            row.names = FALSE, sep = "\t", quote = FALSE)