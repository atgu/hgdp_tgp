import argparse

def main(args):
    # Creating a filemap for all runs (1-10) for all values of k (2-10)
    with open(args.output_path, 'w') as f:
        for k in range(2, 11):
            for run in range(1, 11):
                file_path = f'{args.run_folder}{run}/{args.file_name}.{k}.Q'
                if k == 10 and run == 10:
                    f.write(f'k{k}r{run}\t{k}\t{file_path}')
                else:
                    f.write(f'k{k}r{run}\t{k}\t{file_path}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # file path to the run folders (do not include numbers)
    parser.add_argument('--run_folder',
                        default='/Users/zkoenig/Documents/hgdp_tgp_local/ADMIXTURE/admixture_results_v1/run')
    # name of the .P and .Q files (will be the same name as the plink files used to run ADMIXTURE
    parser.add_argument('--file_name',
                        default='hgdp_tgp')
    # output file path, default will write filepath in current directory
    parser.add_argument('--output_path', default='pong_filemap')

    args = parser.parse_args()

    main(args)
