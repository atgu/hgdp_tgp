# Creating a filemap for all runs (1-10) for all values of k (2-10)
# The filemap is used to run ADMIXTURE
with open('pong_filemap', 'w') as f:
    for k in range(2, 11):
        for run in range(1, 11):
            file_path = f'/Users/zkoenig/Documents/hgdp_tgp_local/ADMIXTURE/admixture_results/' \
                        f'run{run}/post_qc_ld_unrel_only.{k}.Q'
            if k == 10 and run == 10:
                f.write(f'k{k}r{run}\t{k}\t{file_path}')
            else: 
                f.write(f'k{k}r{run}\t{k}\t{file_path}\n')
