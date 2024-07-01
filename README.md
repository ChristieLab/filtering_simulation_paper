# Associated code for: [Next-generation data filtering in the genomics era](www.doi.org/10.1038/s41576-024-00738-6)

This repository contains the supplementary workflow notebooks for the paper [Next-generation data filtering in the genomics era](www.doi.org/10.1038/s41576-024-00738-6), DOI: 10.1038/s41576-024-00738-6, and the associated scripts necissary to run them on a unix-based file system.

It also contains the the scripts used to produce the figure in Box 2 and the figures in the Supplementary material.

To view the workflow notebooks, go to [main/workflow_notebooks/](https://github.com/ChristieLab/filtering_simulation_paper/tree/main/workflow_notebooks) and select either:

1. [pre_variant_notebook.rmd](https://github.com/ChristieLab/filtering_simulation_paper/blob/main/workflow_notebooks/pre_variant_notebook.rmd); or
2. [post_variant_notebook.rmd](https://github.com/ChristieLab/filtering_simulation_paper/blob/main/workflow_notebooks/post_variant_notebook.rmd)

Pre-knit versions of these are also available for the [pre-variant](https://github.com/ChristieLab/filtering_simulation_paper/blob/main/workflow_notebooks/pre_variant_notebook.nb.html) and [post-variant](https://github.com/ChristieLab/filtering_simulation_paper/blob/main/workflow_notebooks/post_variant_notebook.nb.html) workflow notebookss as well from these links or in the same directory as the `.rmd` files. These can be downloaded and viewed with most web browsers.

Each of the shell (`.sh`) scripts in [main/workflow_notebooks/](https://github.com/ChristieLab/filtering_simulation_paper/tree/main/workflow_notebooks) are also necissary for pre-variant notebook, and should be downloaded and placed in the working directory. All other dependencies are described in the notebooks.

Lastly, the scripts needed to produce the Box 2 figure are all contained at [/main/scripts](https://github.com/ChristieLab/filtering_simulation_paper/tree/main/scripts), with the steps required to prepare the raw data and produce the figures described in [Box_1_figure_and_SI.Rmd](https://github.com/ChristieLab/filtering_simulation_paper/blob/main/scripts/Box_1_figure_and_SI.Rmd). Some of the shell (`.sh`) scripts are written for our `slurm` for our specific computing cluster and may need to be tweaked for other systems.
