# HGDP + 1kGP

This repository contains information for QC, analyses, and tutorials relating to the combined HGDP (Human Genome Diversity Project) + 1kGP (1000 Genomes Project) data

## Data access
Metadata available here: gs://hgdp-1kg/tutorial_datasets/metadata_and_qc/gnomad_meta_updated.tsv
- The metadata can be downloaded from google cloud as [described here](https://cloud.google.com/storage/docs/downloading-objects#downloading-an-object) 

All data are freely available and described in more detail [here](https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#the-gnomad-hgdp-and-1000-genomes-callset) 

The gnomAD HGDP+1kGP callset (pre-QC mt) can be found [here](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg)
- Note that files ending with `.bgz` can be viewed using `zcat` on the command line

Datasets used in the tutorials are located [here](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg-tutorials)

Phased haplotypes are available as BCFs on Google Cloud at the following path:
- gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/ 
- More details can be found [here](https://docs.google.com/document/d/1LCx74zREJaJwtN0MzonSv1QB3UahVtgTfjkepXaQUxc/edit)

Datasets found on the Downloads page of the gnomAD browser are released on Google Cloud Platform, Amazon Web Services, and Microsoft Azure. Instructions on how to download them can be found [here](https://gnomad.broadinstitute.org/downloads).  

## Additional information 
PCA plotting and projection scripts available here (used for the [COVID-19 Host Genetics Initiative](https://www.covid19hg.org/), [Global Biobank Meta-analysis Initiative](https://www.globalbiobankmeta.org/), and related projects to align external cohorts to this resource): https://github.com/atgu/pca_projection/blob/master/hgdp_tgp_reference/hgdp_tgp_pca_intersection.py


