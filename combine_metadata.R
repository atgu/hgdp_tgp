library(tidyverse)
library(gtable)
library(maptools)
library(mapdata)
library(raster)
library(ggmap)
library(cowplot)

setwd('/Users/alicia/martin_lab/projects/hgdp_tgp')

# gnomAD call set info for 1kG + HGDP subset of samples
# source of gnomAD IDs, sample QC, etc.
gnomad <- read.delim('gnomad.genomes.v3.1.hgdp_1kg_subset.sample_meta.julia_update.tsv', header=T, sep='\t')
  
# HGDP metadata from Bergstrom et al
# this has the HGDP LP info as well as the source and library_type
hgdp_meta <- read.table('hgdp_wgs.20190516.metadata.txt', header=T) %>%
  separate(library, into=c('hgdp', 'lp'), sep='\\.') %>%
  mutate(s=if_else(grepl('LP', lp), lp, sample)) %>%
  dplyr::select(c(sample, hgdp, lp, source, library_type, region, sex, coverage, freemix, capmq, insert_size_average, array_non_reference_discordance))
names(hgdp_meta) <- paste0('bergstrom.', names(hgdp_meta))

# coordinates of 1000 Genomes (from GGV browser)
tgp_coords <- read.csv('tgp_coords.csv', header=T)

# compiled pop info of HGDP + TGP with harmonized genetic regions
simple_meta <- read.csv('tgp_hgdp.csv', header=T) %>%
  left_join(tgp_coords, by='Population') %>%
  mutate(Latitude = ifelse(is.na(Latitude.x), Latitude.y, Latitude.x),
         Longitude = ifelse(is.na(Longitude.x), Longitude.y, Longitude.x)) %>%
  dplyr::select(-c(ends_with('.x'), ends_with('.y')) ) %>%
  left_join(hgdp_meta, by=c('Sample.ID' = 'bergstrom.sample'))
  
simple_meta_coords <- simple_meta %>% 
  dplyr::select(Project, Study.region, Population, Genetic.region, Latitude, Longitude) %>%
  unique()

write.table(simple_meta_coords, 'tgp_hgdp_coords.tsv', sep='\t', row.names=F, quote=F)


# Create style sheet ------------------------------------------------------

simple_meta_coords <- read.delim('tgp_hgdp_coords2.tsv', sep='\t') %>%
  arrange(Genetic.region, Population)

super_pop_colors <- data.frame(
  Genetic.region=unique(simple_meta_coords$Genetic.region),
  # rough matching with gnomAD scheme
  Continent.colors = c('#984EA3', '#E41A1C', '#FF7F00', '#4DAF4A', '#377EB8', '#A65628', '#999999') 
)

simple_meta_coords <- simple_meta_coords %>% 
  left_join(super_pop_colors) %>%
  group_by(Genetic.region) %>%
  mutate(n=n(),
         rownum = row_number()) %>%
  mutate(Pop.colors = colorRampPalette(c('white', Continent.colors[1], 'black'))(n+2)[rownum+1],
         Pop.shapes = rep(15:19, 5)[1:n])

write.table(simple_meta_coords, 'tgp_hgdp_coords3.tsv', sep='\t', row.names=F, quote=F)
simple_meta_coords <- read.delim('tgp_hgdp_coords3.tsv', sep='\t')

color_vec <- setNames(simple_meta_coords$Pop.colors, simple_meta_coords$Population)
region_vec <- setNames(unique(simple_meta_coords$Continent.colors), unique(simple_meta_coords$Genetic.region))
shape_vec <- setNames(simple_meta_coords$Pop.shapes, simple_meta_coords$Population)

names(simple_meta_coords) <- paste0('hgdp_tgp_meta.', names(simple_meta_coords))

# HGDP+TGP plots ----------------------------------------------------------------

# load in PCA data from the Pan-UKB project (for now)
# TODO: update with PCA from seq data
pca <- read.table(gzfile('/Users/alicia/daly_lab/ukbb_diverse_pops/pca/globalref_ukbb_scores.txt.bgz'), header=T)
pca_meta <- pca %>%
  left_join(simple_meta, by=c('s'='Sample.ID')) %>%
  #left_join(simple_meta_coords, ) %>%
  arrange(Genetic.region, Population)

pca_meta$Population <- factor(pca_meta$Population, levels=unique(pca_meta$Population))

p_pop <- ggplot(pca_meta, aes(x=PC1, y=PC2, color=Population, shape=Population)) +
  geom_point() +
  scale_color_manual(values=color_vec) +
  scale_shape_manual(values=shape_vec) +
  theme_classic()

p_region <- ggplot(pca_meta, aes(x=PC1, y=PC2, color=Genetic.region, shape=Population)) +
  geom_point() +
  scale_color_manual(values=region_vec, name='Genetic region') +
  scale_shape_manual(values=shape_vec) +
  theme_classic() +
  guides(shape=F)


# gnomAD meta combinations ------------------------------------------------

gnomad_hgdp_meta <- gnomad %>%
  inner_join(hgdp_meta, by=c('project_meta.sample_id'='bergstrom.sample'))
gnomad_hgdp_meta2 <- gnomad %>%
  inner_join(hgdp_meta, by=c('project_meta.sample_id'='bergstrom.lp'))
gnomad_hgdp_neither <- gnomad %>%
  filter(subsets.hgdp == 'true') %>%
  subset(!(project_meta.sample_id %in% hgdp_meta$bergstrom.hgdp) & !(project_meta.sample_id %in% hgdp_meta$bergstrom.lp))
gnomad_tgp <- gnomad %>% filter(subsets.hgdp == 'false')
gnomad_meta <- bind_rows(gnomad_hgdp_meta, gnomad_hgdp_meta2, gnomad_hgdp_neither, gnomad_tgp)

gnomad_meta_coords <- gnomad_meta %>%
  inner_join(simple_meta_coords, by=c('project_meta.project_subpop' = 'hgdp_tgp_meta.project_meta.project_subpop'))

write.table(gnomad_meta_coords, 'gnomad_meta_v1.txt', row.names=F, quote=F, sep='\t')


# Plot map of samples -----------------------------------------------------

data(wrld_simpl)
world <- fortify(wrld_simpl)

lims = SpatialPoints(coords = data_frame(x = xlim, y = ylim), proj4string = CRS("+proj=longlat +datum=WGS84"))%>%
  spTransform(CRS("+init=epsg:3857"))

# Read the latitude/longitdue data for populations in AGVP/TGP
pop_pos <- read.csv('african_pop_summary3.csv', header=T)#, stringsAsFactors = F)
pop_pos_plot <- subset(pop_pos, Population != 'Ethiopia')
pop_pos_plot$Population <- factor(pop_pos_plot$Population, levels=as.character(pop_pos_plot$Population))
pop_pos_plot <- pop_pos_plot %>% 
  filter(Population %in% pca_tgp_agvp$FID) %>%
  droplevels()
#pop_pos_plot <- subset(pop_pos_plot, Population %in% pca$FID)
color_vec <- as.character(pop_pos$Color)
names(color_vec) <- pop_pos$Population
shape_vec <- pop_pos$Shape
names(shape_vec) <- pop_pos$Population

gnomad_meta_coords_pop <- gnomad_meta_coords %>% 
  dplyr::select(starts_with('hgdp_tgp_meta')) %>% 
  group_by(hgdp_tgp_meta.Population) %>%
  mutate(hgdp_tgp_meta.n=n()) %>%
  unique() %>%
  ungroup()
# plot the map of Africa with data points labeled
p_world <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group=group), fill='lightyellow', color='lightgrey') +
  geom_point(data = gnomad_meta_coords_pop, aes(hgdp_tgp_meta.Longitude, hgdp_tgp_meta.Latitude, 
                                                color=hgdp_tgp_meta.Genetic.region, fill=hgdp_tgp_meta.Genetic.region, 
                                                shape=hgdp_tgp_meta.Project, size=hgdp_tgp_meta.n)) +
  coord_fixed(xlim = c(-165,165), ylim = c(-55,85)) +
  labs(x='Longitude', y='Latitude') +
  theme_classic() +
  scale_fill_manual(name = "Region", values = region_vec) +
  scale_color_manual(name = "Region", values = region_vec) +
  scale_size(name = 'N') +
  scale_shape(name='Project') +
  theme(panel.background = element_rect(fill = "lightblue"),
        plot.background = element_rect(fill = "transparent", color = NA),
        #legend.justification = c(0,1),
        text = element_text(size=14),
        axis.text = element_text(color='black'))

ggsave('hgdp_tgp_world_map.pdf', p_world, width=10, height=5)
