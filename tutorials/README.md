## Overview

The HGDP (Human Genome Diversity Project) and 1kGP (1000 Genomes Project) are the two largest public genetic datasets with individual level data. In this project, we have joint called these datasets integrated with gnomAD and provide a harmonized resource ready for downstream analysis or integration with custom datasets. We provide a set of tutorials of our analyses implemented primarily in Hail as educational resources.

## Software Requirements

The tutorials have the following software requirements if you would like to run them as written using Google Cloud:
- [Hail](https://hail.is/#install)
- [Google Cloud sdk](https://cloud.google.com/sdk/docs/install)
- [Python 3.7 or later](https://www.python.org/downloads/)
- [Plotly](https://plotly.com/python/getting-started/)
- [Pandas](https://pandas.pydata.org/getting_started.html)
- [Pickle](https://docs.python.org/3/library/pickle.html#module-pickle)
- [Scikit-learn Ensemble](https://scikit-learn.org/stable/modules/ensemble.html)
- [Typing](https://docs.python.org/3/library/typing.html)

To run the tutorials, the user will need to start up a Google Cloud cluster on their own Google Cloud project. More information on creating Google Cloud projects can be found [here](https://cloud.google.com/resource-manager/docs/creating-managing-projects).

### How to start a Google Cloud cluster:
- The tutorials use hail, so we recommend starting clusters using `hailctl`. 
    - More information on using hail on the cloud can be found [here](https://hail.is/docs/0.2/hail_on_the_cloud.html).
    -  More information on `hailctl` can be found [here](https://hail.is/docs/0.2/cloud/google_cloud.html#hailctl-dataproc). 
- To run the tutorials with the gnomAD components you will need to add the following arguments to your `hailctl dataproc` command:
    - `--packages gnomad` OR `--packages "git+https://github.com/broadinstitute/gnomad_methods.git@main"`
- **If you plan to run Notebook 5**, we suggest you start your cluster with the following commmand:
    >`hailctl dataproc start qc-notebook5 --project [YOUR_PROJECT_NAME] --num-secondary-workers 50 --region=us-central1 --zone=us-central1-b --packages git+https://github.com/broadinstitute/gnomad_methods.git --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --big-executors`

## Links to view notebooks
To view the tutorials rendered so that jupyter notebook links function properly click on the links below for each notebook:
<br>[Notebook 1](https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb1.ipynb) - Metadata and QC
<br>[Notebook 2](https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb2.ipynb) - PCA and Ancestry Analyses
<br>[Notebook 3](https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb3.ipynb) - Summarizing Data Post QC
<br>[Notebook 4](https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb4.ipynb) - Computing Population Genetics Statistics (F2 and F<sub>st</sub>)
<br>[Notebook 5](https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb5.ipynb) - Assigning Ancestry Labels Using a Random Forest Model