import hail as hl

def read_qc(
        raw: bool = False,
        post_qc:bool = False,
        sample_qc: bool = False,
        variant_qc: bool = False,
        outlier_removal: bool = False,
        ld_pruning: bool = False,
        rel_unrel: str = 'default',
        n_partitions: int = 0) -> hl.MatrixTable:
    """
    Wrapper function to get HGDP+1kGP data as Matrix Table at different stages of QC/filtering.
    By raw, returns pre QC MatrixTable with qc filters annotated but not filtered.

    :param bool raw: if True will return a preQC version of the dataset
    :param bool post_qc: if True will return a post QC matrix table that has gone through:
        - sample QC
        - variant QC
        - duplicate removal
        - outlier removal
    :param bool sample_qc: if True will return a post sample QC matrix table
    :param bool variant_qc: if True will return a post variant QC matrix table
    :param bool outlier_removal: if True will return a matrix table with PCA outliers removed
    :param bool ld_pruning: if True will return a matrix table that has gone through:
        - sample QC
        - variant QC
        - duplicate removal
        - LD pruning
        - additional variant filtering
    :param bool rel_unrel: default will return same mt as ld pruned above
        if 'all' will return the same matrix table as if ld_pruning is True
        if 'related_pre_outlier' will return a matrix table with only related samples pre pca outlier removal
        if 'unrelated_pre_outlier' will return a matrix table with only unrelated samples pre pca outlier removal
        if 'related_post_outlier' will return a matrix table with only related samples post pca outlier removal
        if 'unrelated_post_outlier' wil return a matrix table with only unrelated samples post pca outlier removal
    :param int n_partitions: if specified, will read in dataset with given number of partitions for the following arguments:
        - ld_pruning
        - rel_unrel
    """
    # Reading in all the tables and matrix tables needed to generate the pre_qc matrix table
    sample_meta = hl.import_table('gs://hgdp-1kg/hgdp_tgp/qc_and_figure_generation/gnomad_meta_v1.tsv')
    sample_qc_meta = hl.read_table('gs://hgdp_tgp/output/gnomad_v3.1_sample_qc_metadata_hgdp_tgp_subset.ht')
    dense_mt = hl.read_matrix_table(
        'gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt')
    
    dense_mt = dense_mt.naive_coalesce(5000)


    # Takes a list of dicts and converts it to a struct format (works with nested structs too)
    def dict_to_struct(d):
        fields = {}
        for k, v in d.items():
            if isinstance(v, dict):
                v = dict_to_struct(v)
            fields[k] = v
        return hl.struct(**fields)

    # un-flattening a hail table with nested structure
    # dict to hold struct names as well as nested field names
    d = {}

    # Getting the row field names
    row = sample_meta.row_value

    # returns a dict with the struct names as keys and their inner field names as values
    for name in row:
        def recur(dict_ref, split_name):
            if len(split_name) == 1:
                dict_ref[split_name[0]] = row[name]
                return
            existing = dict_ref.get(split_name[0])
            if existing is not None:
                assert isinstance(existing, dict), existing
                recur(existing, split_name[1:])
            else:
                existing = {}
                dict_ref[split_name[0]] = existing
                recur(existing, split_name[1:])
        recur(d, name.split('.'))

    # using the dict created from flattened struct, creating new structs now un-flattened
    sample_meta = sample_meta.select(**dict_to_struct(d))
    sample_meta = sample_meta.key_by('s')

    # grabbing the columns needed from HGDP metadata
    new_meta = sample_meta.select(sample_meta.hgdp_tgp_meta, sample_meta.bergstrom)

    # creating a table with gnomAD sample metadata and HGDP metadata
    ht = sample_qc_meta.annotate(**new_meta[sample_qc_meta.s])

    # stripping 'v3.1::' from the names to match with the densified MT
    ht = ht.key_by(s=ht.s.replace("v3.1::", ""))

    # Using hl.annotate_cols() method to annotate the gnomAD variant QC metadata onto the matrix table
    mt = dense_mt.annotate_cols(**ht[dense_mt.s])
    
    if raw:
        print("Returning default preQC matrix table")
        # returns preQC dataset
        return mt
    
    if post_qc:
        print("Returning post sample and variant QC matrix table with duplicates and PCA outliers removed")
        sample_qc = True
        variant_qc = True
        duplicate = True
        outlier_removal = True
    
    if sample_qc:
        print("Applying sample QC")
        # Apply sample QC filters to dataset
        # filtering samples to those who should pass gnomADs sample QC
        # this filters to only samples that passed gnomad sample QC hard filters
        mt = mt.filter_cols(~mt.sample_filters.hard_filtered)

    if variant_qc:
        print("Applying variant QC")
        # Apply variant QC filters to dataset
        # Subsetting the variants in the dataset to only PASS variants (those which passed gnomAD's variant QC)
        # PASS variants are variants which have an entry in the filters field.
        # This field contains an array which contains a bool if any variant qc filter was failed
        # This is the last step in the QC process
        mt = mt.filter_rows(hl.len(mt.filters) != 0, keep=False)

    if outlier_removal:
        print("Removing PCA outliers")
        # remove PCA outliers
        # reading in the PCA outlier list
        # To read in the PCA outlier list, first need to read the file in as a list
        # using hl.hadoop_open here which allows one to read in files into hail from Google cloud storage
        pca_outlier_path = 'gs://hgdp-1kg/hgdp_tgp/pca_outliers_v2.txt'
        with hl.utils.hadoop_open(pca_outlier_path) as file:
            outliers = [line.rstrip('\n') for line in file]

        # Using hl.literal here to convert the list from a python object to a hail expression so that it can be used
        # to filter out samples
        outliers_list = hl.literal(outliers)

        # Using the list of PCA outliers, using the ~ operator which is a negation operator and obtains the compliment
        # In this case the compliment is samples which are not contained in the pca outlier list
        mt = mt.filter_cols(~outliers_list.contains(mt['s']))

    if ld_pruning:
        print("Returning ld pruned post variant and sample QC matrix table pre PCA outlier removal ")
        # read in dataset which has additional variant filtering and ld pruning run
        # data has gone through:
        #   - sample QC
        #   - variant QC
        if n_partitions != 0:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/intermediate_files/filtered_n_pruned_output_updated.mt',
            _n_partitions = n_partitions)
        else:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/intermediate_files/filtered_n_pruned_output_updated.mt')

    if rel_unrel == "default":
        # do nothing
        # created a default value because there are multiple options for rel/unrel datasets
        mt = mt

    elif rel_unrel == 'related_pre_outlier':
        print("Returning post sample and variant QC matrix table " \
              "pre PCA outlier removal with only related individuals")
        # data has gone through:
        #   - sample QC
        #   - variant QC
        #   - duplicate removal
        #   - LD pruning
        #   - pc_relate 
        #   - filter to only related individuals   
        if n_partitions != 0:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/rel_updated.mt',
            _n_partitions = n_partitions)
        else:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/rel_updated.mt')
        
    elif rel_unrel == 'unrelated_pre_outlier':
        print("Returning post sample and variant QC matrix table " \
              "pre PCA outlier removal with only unrelated individuals")
        # data has gone through:
        #   - sample QC
        #   - variant QC
        #   - duplicate removal
        #   - LD pruning
        #   - pc_relate 
        #   - filter to only unrelated individuals
        if n_partitions != 0:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/unrel_updated.mt',
            _n_partitions = n_partitions)
        else:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/unrel_updated.mt')

    elif rel_unrel == 'related_post_outlier':
        print("Returning post sample and variant QC matrix table " \
              "post PCA outlier removal with only related individuals")
        # data has gone through:
        #   - sample QC
        #   - variant QC
        #   - duplicate removal
        #   - LD pruning
        #   - pc_relate 
        #   - filter to only related individuals
        #   - PCA outlier removal
        if n_partitions != 0:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/related.mt',
            _n_partitions = n_partitions)
        else:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/related.mt')

    elif rel_unrel == 'unrelated_post_outlier':
        print("Returning post sample and variant QC matrix table " \
              "post PCA outlier removal with only unrelated individuals")
        # data has gone through:
        #   - sample QC
        #   - variant QC
        #   - duplicate removal
        #   - LD pruning
        #   - pc_relate 
        #   - filter to only unrelated individuals
        #   - PCA outlier removal
        if n_partitions != 0:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/unrelated.mt',
            _n_partitions = n_partitions)
        else:
            mt = hl.read_matrix_table('gs://hgdp-1kg/hgdp_tgp/datasets_for_others/lindo/ds_without_outliers/unrelated.mt')
        
    # Calculating both variant and sample_qc metrics on the mt before returning
    # so the stats are up to date with the version being written out
    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)
    
    return mt