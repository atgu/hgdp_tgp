library(tidyverse)
library(plotly)
library(RColorBrewer)
library(cowplot)
library(GGally)
library(gtable)
library(maptools)
library(mapdata)
library(raster)
library(ggmap)
library(ggrepel)
library(data.table)


# Make map plot of HGDP+1kG from metadata info ----------------------------

setwd('/Users/alicia/martin_lab/projects/hgdp_tgp')

hgdp_tgp <- read.delim('gnomad_meta_hgdp_tgp_v1.txt', header=T, sep='\t') %>%
  dplyr::select(s, project_meta.title, starts_with('hgdp_tgp_meta'), starts_with('bergstrom')) 

outliers <- read.table('scripts/hgdp_tgp/pca_outliers.txt')$V1 # determined iteratively by subcontinental vs global PCA

data(wrld_simpl)
world <- fortify(wrld_simpl)

region_vec <- setNames(unique(hgdp_tgp$hgdp_tgp_meta.Continent.colors), 
                       unique(hgdp_tgp$hgdp_tgp_meta.Genetic.region))

simple_meta <- read.delim('hgdp_tgp_pop_summary.txt', sep='\t') %>%
  left_join(hgdp_tgp) %>%
  dplyr::select(hgdp_tgp_meta.Population, hgdp_tgp_meta.Longitude, hgdp_tgp_meta.Latitude, hgdp_tgp_meta.Genetic.region, hgdp_tgp_meta.Project, count) %>%
  unique()

p_world <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group=group), fill='lightyellow', color='lightgrey') +
  geom_point(data = simple_meta, aes(hgdp_tgp_meta.Longitude, hgdp_tgp_meta.Latitude,
                                  color=hgdp_tgp_meta.Genetic.region, fill=hgdp_tgp_meta.Genetic.region,
                                  shape=hgdp_tgp_meta.Project, size=count)) +
  coord_fixed(xlim = c(-165,165), ylim = c(-55,85)) +
  labs(x='Longitude', y='Latitude') +
  theme_classic() +
  scale_fill_manual(name = "Region", values = region_vec) +
  scale_color_manual(name = "Region", values = region_vec) +
  scale_shape(name='Project') +
  scale_size(name = 'N', range=c(0,10)) +
  guides(color=F, fill=F, shape=F) +
  theme(panel.background = element_rect(fill = "lightblue"),
        plot.background = element_rect(fill = "transparent", color = NA),
        #legend.justification = c(0,1),
        text = element_text(size=14),
        axis.text = element_text(color='black'))


# Make global PCA plots ---------------------------------------------------

setwd("/Users/alicia/martin_lab/projects/hgdp_tgp/pca/global/")

ref <- read.table(gzfile('hgdp_tgp_pca_snps_scores.txt.bgz'), header=T)
proj <- read.table(gzfile('projected_scores.txt.bgz'), header=T)

# combine the ref and proj data into one dataset and then join that to the metadata to add additional info for those samples 
pca <- ref %>%
  left_join(hgdp_tgp) %>%
  mutate(continent_name=ifelse(s %in% outliers, 'outlier', hgdp_tgp_meta.Genetic.region),
         continent_color=ifelse(s %in% outliers, 'black', hgdp_tgp_meta.Continent.colors))
color_vec <- pca %>%
  dplyr::select(continent_name, continent_color) %>%
  unique() %>%
  deframe()

global_pca_plot <- function(pca, first_PC, second_PC, guides=F) {
  p_global <- ggplot(pca, aes_string(x=first_PC, y=second_PC, color='continent_name', shape='hgdp_tgp_meta.Project')) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values=color_vec, name = "Region") +
    scale_shape(name='Project') +
    theme(text = element_text(size=14),
          axis.text = element_text(color='black'))
  if(!guides) {
    p_global <- p_global + guides(shape=F, color=F)
  }
  return(p_global)
}

p_global_pc_legend <- global_pca_plot(pca, 'PC1', 'PC2', guides=T)
p_global_pc1_2 <- global_pca_plot(pca, 'PC1', 'PC2')
p_global_pc3_4 <- global_pca_plot(pca, 'PC3', 'PC4')
p_global_pc5_6 <- global_pca_plot(pca, 'PC5', 'PC6')

global_legend <- get_legend(p_global_pc_legend)

p_global_legend <- plot_grid(p_global_pc1_2, p_global_pc3_4, global_legend, rel_widths = c(1,1,.4), nrow=1, labels=c('B', 'C'))
p_global_legend2 <- plot_grid(p_global_pc1_2, global_legend, rel_widths = c(1,.4), nrow=1)

save_plot('global_pca.pdf', p_global_legend, base_height=4, base_width=10)
save_plot('global_pca2.pdf', p_global_legend2, base_height=4, base_width=5)

p_map_pcs <- plot_grid(p_world, p_global_legend, labels=c('A', ''), ncol=1)
save_plot('global_map_pca.pdf', p_map_pcs, base_height=8, base_width=10)
save_plot('global_map_pca.png', p_map_pcs, base_height=8, base_width=10)

# Make subcontinental PCA plots -------------------------------------------

# set working directory and import metadata

setwd("/Users/alicia/martin_lab/projects/hgdp_tgp/pca/subcont/")

single_pca_plot <- function(pca_subcont, first_PC, second_PC, pop_name, bottom=F) {
  nb.cols <- length(unique(pca_subcont$hgdp_tgp_meta.Population))
  mycolors <- colorRampPalette(brewer.pal(7, "Set1"))(nb.cols)
  is_tgp <- unique(pca_subcont$hgdp_tgp_meta.Population) %in% subset(pca_subcont, hgdp_tgp_meta.Project=='1000 Genomes')$hgdp_tgp_meta.Population
  myshapes <- ifelse(is_tgp, 16, 4)
  
  p_subcont <- ggplot(pca_subcont, aes_string(x=first_PC, y=second_PC, color='hgdp_tgp_meta.Population', shape='hgdp_tgp_meta.Population', text='s')) +
    geom_point(size=2) +
    theme_classic() +
    scale_shape_manual(values=myshapes, labels=unique(pca_subcont$hgdp_tgp_meta.Population), name='Population') +
    scale_color_manual(values=mycolors, labels=unique(pca_subcont$hgdp_tgp_meta.Population), name='Population') +
    labs(title=pop_name) +
    theme(text = element_text(size=14, color='black'),
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 9),
          legend.key.size = unit(0.7, "lines")) +
    guides(shape = guide_legend(override.aes = list(size = 0.85)))
  if(bottom) {
    p_subcont <- p_subcont +
      theme(axis.text = element_text(color='black'),
            legend.title = element_text(size = 12), 
            legend.text = element_text(size = 9),
            legend.key.size = unit(0.7, "lines"),
            legend.position="bottom") +
      guides(fill=guide_legend(ncol=1,byrow=FALSE, title.position="top"),
             shape=guide_legend(ncol=1,byrow=FALSE, title.position="top"))
  }
  return(p_subcont)
}


load_pop_make_pca_plot <- function(pop_name) {
  # Subcont PCA 
  ref_unrel <- read.table(gzfile(paste0('subcont_pca_', pop_name, '_scores.txt.bgz')), header=T)
  proj_rel <- read.table(gzfile(paste0('subcont_pca_', pop_name, '_projected_scores.txt.bgz')), header=T)
  
  # combine the ref and proj data into one dataset and then join that to the metadata to add additional info for those samples 
  pca_subcont <- bind_rows(ref_unrel, proj_rel) %>%
    left_join(hgdp_tgp) %>%
    arrange(hgdp_tgp_meta.Population)
  
  p_1_2 <- single_pca_plot(pca_subcont, 'PC1', 'PC2', pop_name)
  p_1_2_bottom <- single_pca_plot(pca_subcont, 'PC1', 'PC2', pop_name, TRUE) 
  p_3_4 <- single_pca_plot(pca_subcont, 'PC3', 'PC4', pop_name)
  p_5_6 <- single_pca_plot(pca_subcont, 'PC5', 'PC6', pop_name)
  p_7_8 <- single_pca_plot(pca_subcont, 'PC7', 'PC8', pop_name)
  p_subcont <- plot_grid(p_1_2 + theme(legend.position="none"), p_3_4 + theme(legend.position="none"), 
                         p_5_6 + theme(legend.position="none"), p_7_8 + theme(legend.position="none"), nrow=1)
  legend <- get_legend(p_1_2)
  p_subcont_legend <- plot_grid(p_subcont, legend, rel_widths = c(4, 1))
  
  ggsave(paste0('subcont_pca_', pop_name, '.pdf'), p_subcont_legend, width=15, height=4)
  
  ggplotly(p_subcont_legend, tooltip='text')
  return(list(p_subcont_legend, p_1_2_bottom))
}

p_AFR <- load_pop_make_pca_plot('AFR') # outlier: NA20314 (PC1/2), NA20299 (PC3/4), HG01880 (PC6), HG01881 (PC6)
p_AMR <- load_pop_make_pca_plot('AMR')
p_CSA <- load_pop_make_pca_plot('CSA') # outlier: HGDP00130 (PC1/2), HGDP00013 (PC4&5), HGDP00150 (PC5), HGDP00029
p_EAS <- load_pop_make_pca_plot('EAS') # outlier: HGDP01298 (PC3), HGDP01303 (PC4&5), LP6005443-DNA_B02 (PC5), HGDP01300 (PC6)
p_EUR <- load_pop_make_pca_plot('EUR') # outlier: HG01628 (PC5), HG01629 (PC5), HG01630 (PC5), HG01694 (PC6), HG01696 (PC6)
p_MID <- load_pop_make_pca_plot('MID') # outlier: HGDP00621, HGDP01270, HGDP01271
p_OCE <- load_pop_make_pca_plot('OCE') # outlier: HGDP00554 - not actually an outlier when looking at global PCs

p_region_outliers <- plot_grid(p_AFR[[1]], p_AMR[[1]], p_CSA[[1]], p_EAS[[1]], p_EUR[[1]], p_MID[[1]], p_OCE[[1]], 
                               labels=LETTERS[1:7], ncol=1, align='v')
p_region_outliers_1_2 <- plot_grid(p_AFR[[2]], p_AMR[[2]], p_CSA[[2]], p_EAS[[2]], p_EUR[[2]], p_MID[[2]], p_OCE[[2]], nrow=1, align='h', labels=LETTERS[5:11])
ggsave('region_pca_1_10_no_outliers.pdf', p_region_outliers, height=28, width=15)
ggsave('region_pca_1_10_no_outliers.png', p_region_outliers, height=15, width=12)


# ADMIXTURE analysis ------------------------------------------------------

setwd('/Users/alicia/martin_lab/projects/hgdp_tgp/admixture')

fam <- read.table('unrel_pass_no_outliers.fam') %>%
  dplyr::select(V2) %>%
  rename(s=V2)
plot_admixture <- function(K, to_label=T) {
  Q <- read.table(paste0('unrel_pass_no_outliers.', K, '.Q')) %>%
    mutate(k=paste0('K=',K)) %>%
    bind_cols(fam) %>%
    setnames(old=paste0('V', 1:K), new=paste0('component', 1:K)) %>%
    #rename(component1=V1, component2=V2) %>%
    pivot_longer(paste0('component', 1:K)) %>%
    left_join(hgdp_tgp) %>%
    arrange(hgdp_tgp_meta.Genetic.region, hgdp_tgp_meta.Population) %>%
    #group_by(s) %>%
    mutate(pos= ceiling(seq_along(s)/K))
  
  Q$hgdp_tgp_meta.Genetic.region <- factor(Q$hgdp_tgp_meta.Genetic.region, levels=unique(Q$hgdp_tgp_meta.Genetic.region))
  Q$hgdp_tgp_meta.Population <- factor(Q$hgdp_tgp_meta.Population, levels=unique(Q$hgdp_tgp_meta.Population))
  Q$s <- factor(Q$s, levels=unique(Q$s))
  
  pop_mids <- Q %>%
    group_by(hgdp_tgp_meta.Genetic.region, hgdp_tgp_meta.Population) %>%
    mutate(index=row_number()) %>%
    filter(row_number()==ceiling(n()/K)) %>%
    dplyr::select(hgdp_tgp_meta.Genetic.region, hgdp_tgp_meta.Population, s, index) %>%
    ungroup() %>%
    mutate(cum_index = cumsum(index))
  
  pop_maxs <- Q %>%
    group_by(hgdp_tgp_meta.Genetic.region, hgdp_tgp_meta.Population) %>%
    filter(row_number()==ceiling(n())) %>%
    dplyr::select(hgdp_tgp_meta.Genetic.region, hgdp_tgp_meta.Population, s) %>%
    ungroup()
  
  pop_label <- rep("", nrow(Q)/K)
  pop_label[pop_mids$cum_index] <- as.character(pop_mids$hgdp_tgp_meta.Population)
  
  p_Q <- ggplot(Q, aes(x=s, y=value, fill=name)) +
    geom_vline(xintercept=pop_maxs$s) +
    geom_col(size=1) +
    facet_grid(k~fct_inorder(hgdp_tgp_meta.Genetic.region), scales='free', space='free') +
    theme_classic() + 
    #labs(x='Populations', y='Ancestry') +
    labs(x='', y='Ancestry') +
    scale_fill_brewer(palette = 'Set1') +
    #scale_x_discrete(labels = pop_mids$hgdp_tgp_meta.Population, breaks=pop_mids$s) +
    #scale_x_discrete(labels = rep("", nrow(pop_mids)), breaks=pop_mids$s) +
    guides(fill=F) +
    theme(
      # axis.text.x = element_text(angle=90, hjust=1),
      axis.text.x = element_blank(),
      axis.text = element_text(color='black'),
      text = element_text(size=14))
  # p_label <- ggplot(Q, aes(x=s, y=value, fill=name)) +
  #   geom_text_repel()
  p_label <- ggplot(data.frame(x = 1:nrow(Q)/K,
                            pop = pop_label),
                 aes(x = x, y = 1, label = pop)) +
    geom_text_repel(min.segment.length = grid::unit(0, "pt"),
                    angle=90,
                    #color = "grey30",  ## ggplot2 theme_grey() axis text
                    size = 2  ## ggplot2 theme_grey() axis text
    ) +
    # scale_x_continuous(limits = c(0, nrow(Q)/K), expand = c(0, 0),
    #                    breaks = NULL, labels = NULL, name = NULL) +
    scale_x_continuous(limits = c(0, nrow(Q)/2), expand=c(0,0),
                       breaks = NULL, labels = NULL, name = NULL) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                       breaks = NULL, labels = NULL, name = NULL) +
    theme(panel.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt"))
  
  if(to_label) {
    p_comb <- plot_grid(p_Q, p_label, align='v', axis='lr', ncol=1)
  } else {
    p_comb <- p_Q
  }
  return(p_comb)
}

p_q2 <- plot_admixture(2, to_label=F)
p_q3 <- plot_admixture(3, to_label=F)
p_q4 <- plot_admixture(4, to_label=F)
p_q5 <- plot_admixture(5, to_label=F)
p_q6 <- plot_admixture(6, to_label=F)
p_q7 <- plot_admixture(7, to_label=F)
p_q8 <- plot_admixture(8, to_label=F)
p_q9 <- plot_admixture(9, to_label=F)
save_plot('admixture_k2.pdf', p_q2, base_width=12, base_height=4)
save_plot('admixture_k3.pdf', p_q3, base_width=12, base_height=4)

across_Ks <- plot_grid(p_q2, p_q3, p_q4, p_q5, p_q6, p_q7, p_q8, p_q9, ncol=1)
  
p_map_pcs_admixture <- plot_grid(p_map_pcs, p_q2, labels=c('', 'D'), ncol=1, rel_heights=c(5, 1.5))
p_map_pcs_admixture_subcont <- plot_grid(p_map_pcs_admixture, p_region_outliers_1_2, ncol=1, rel_heights=c(7, 4))

save_plot('global_map_pca_admixture.pdf', p_map_pcs_admixture, base_height=12, base_width=12)
save_plot('global_map_pca_admixture_subcont.pdf', p_map_pcs_admixture_subcont, base_height=14, base_width=10)

save_plot('admixture_K2_9.pdf', across_Ks, base_height=12, base_width=12)
save_plot('admixture_K2_9.png', across_Ks, base_height=12, base_width=12)

