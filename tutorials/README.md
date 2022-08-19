To view the tutorials rendered so that included links function properly click on the links below for each notebook:
<br>Notebook 1: [link](https://nbviewer.org/github/atgu/hgdp_tgp/blob/tutorial_reformat/tutorials/nb1.ipynb)
<br>Notebook 2: [link](https://nbviewer.org/github/atgu/hgdp_tgp/blob/tutorial_reformat/tutorials/nb2.ipynb)
<br>Notebook 3: [link](https://nbviewer.org/github/atgu/hgdp_tgp/blob/tutorial_reformat/tutorials/nb3.ipynb)
<br>Notebook 4: [link](https://nbviewer.org/github/atgu/hgdp_tgp/blob/tutorial_reformat/tutorials/nb4.ipynb)
<br>Notebook 5: [link_TBD]()

The tutorials have the following software requirements if you would like to run them as written using google cloud:
- [Hail](https://hail.is/#install)
- [Google Cloud sdk](https://cloud.google.com/sdk/docs/install)
- [Python 3.7 or later](https://www.python.org/downloads/)
- 

To run the tutorials, the user will need to start up a google cloud cluster on their own google cloud project. More information on creating google cloud projects can be found [here](https://cloud.google.com/resource-manager/docs/creating-managing-projects).

How to start a google cloud cluster:
- Since the tutorials use hail, a cluster will be started using `hailctl`. More information on hail on the cloud can be found [here](https://hail.is/docs/0.2/hail_on_the_cloud.html).
- To start a cluster you will need to use the command line tool [hailctl](https://hail.is/docs/0.2/cloud/google_cloud.html#hailctl-dataproc). 
- To run the tutorials with the gnomAD components you will need to add the following arguments to your `hailctl dataproc` command:
    - `--packages gnomad`
    - `--requester-pays-allow-all`