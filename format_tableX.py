# The purpose of this script is to format and write out a matrix table which will be used to create 'table_x'
# for our resource manuscript
# author: Zan Koenig

import hail as hl
hl.init()

# reading in the post QC version of the merged dataset (with metadata)
mt = hl.read_matrix_table('gs://african-seq-data/hgdp_tgp/hgdp_tgp_postQC.mt')

# Running sample_qc to get the n_snp and n_singleton counts
mt = hl.sample_qc(mt, name = "new_sample_qc")

# Grabbing only the columns from the matrix table (outputs table of just columns)
col_table = mt.cols()


# writing out a col table with only the columns needed for table x
col_table = col_table.select(col_table.hgdp_tgp_meta.Study.region,
                             col_table.hgdp_tgp_meta.Population,
                             col_table.new_sample_qc.n_snp,
                             col_table.new_sample_qc.n_singleton,
                             col_table.bam_metrics.mean_coverage)


# writing out col_table as a checkpoint to make the downstream steps run faster
col_table.checkpoint('gs://african-seq-data/hgdp_tgp/table_x_checkpoint.ht')


# this is a table of only the columns with only information
col_table = hl.read_table('gs://african-seq-data/hgdp_tgp/table_x_checkpoint.ht')


# Need to get number of unrelateds annotated to the table
# Reading in the unrelated and related matrix tables
unrelated = hl.read_matrix_table('gs://african-seq-data/hgdp_tgp/unrel.mt')
related = hl.read_matrix_table('gs://african-seq-data/hgdp_tgp/rel.mt')

unrelated = unrelated.annotate_cols(unrelated = True)
related = related.annotate_cols(unrelated = False)


unrelated_cols = unrelated.cols()
related_cols = related.cols()

# Annotating the related/unrelated mts with counts per population
related_cols.aggregate(hl.agg.counter(related_cols.hgdp_tgp_meta.Population))
unrelated_cols.aggregate(hl.agg.counter(unrelated_cols.hgdp_tgp_meta.Population))


mt_rel = unrelated.union_cols(related)


mt_rel.aggregate_cols(hl.agg.counter(mt_rel.unrelated))


rel_table = mt_rel.cols()


col_table = col_table.annotate(unrel = rel_table[col_table.s].unrelated)


# Creating three separate tables for each column
# Calculating stats for that column grouped by genetic region (pop)
n_snp = col_table.group_by(
    col_table.region, col_table.Population).aggregate(
    n_snp_stats = hl.agg.stats(col_table.n_snp),
    n_unrelated = hl.agg.count_where(col_table.unrel == True))

n_singleton = col_table.group_by(
    col_table.region, col_table.Population).aggregate(
    n_singleton_stats = hl.agg.stats(col_table.n_singleton),
    n_unrelated = hl.agg.count_where(col_table.unrel == True))

mean_coverage = col_table.group_by(
    col_table.region, col_table.Population).aggregate(
    cov_stats = hl.agg.stats(col_table.mean_coverage),
    n_unrelated = hl.agg.count_where(col_table.unrel == True))



mean_coverage = mean_coverage.flatten()

n_snp = n_snp.flatten()

n_singleton = n_singleton.flatten()


mean_coverage = mean_coverage.select('region',
                                     'Population',
                                     'n_unrelated',
                                     'cov_stats.n',
                                     'cov_stats.mean',
                                     'cov_stats.stdev')

n_snp = n_snp.select('region',
                     'Population',
                     'n_unrelated',
                     'n_snp_stats.n',
                     'n_snp_stats.mean',
                     'n_snp_stats.stdev')

n_singleton = n_singleton.select('region',
                                 'Population',
                                 'n_unrelated',
                                 'n_singleton_stats.n',
                                 'n_singleton_stats.mean',
                                 'n_singleton_stats.stdev')



mean_coverage = mean_coverage.key_by('region','Population')

n_snp = n_snp.key_by('region', 'Population')

n_singleton = n_singleton.key_by('region', 'Population')


table = mean_coverage.annotate(n_snp = n_snp[mean_coverage.region, mean_coverage.Population])
table = table.annotate(n_singleton = n_singleton[table.region, table.Population])

# Flattening out the structs created from annotating the tables
table = table.flatten()

table = table.key_by(table.region, table.Population)


table.export('gs://african-seq-data/hgdp_tgp/table_x.tsv', header=True)
