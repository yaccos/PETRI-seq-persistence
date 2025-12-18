# Petrisnake: A secondary analysis pipeline for PETRI-seq data.

This is a Snakemake pipeline for the secondary computational analysis of single cell RNA-seq data from the PETRI-seq protocol (https://www.nature.com/articles/s41564-020-0729-6 and https://www.nature.com/articles/s41586-024-08124-2).

![Workflow rulegraph](./images/rulegraph.svg)

# Dependencies

This pipeline is only tested under Linux running on x86-64 architecture, although it might work on macOS and other CPU architectures if some work is put into dependency management. The workflow will complain, but still work if the system dependency `fuse2fs` is not installed. For running the pipeline the following dependencies are required (tested versions):

* [Snakemake](https://snakemake.github.io/) (9.14.5)
* [conda](https://www.anaconda.com/docs/getting-started/miniconda/main) (25.7.0)
* [Apptainer](https://apptainer.org/) (1.4.5)

The configuration tool GUI is a R Shiny application with the following dependencies:

* R (4.5.2)
* shiny (1.12.1)
* shinyBS (0.61.1)
* purrr (1.2.0)
* yaml (2.3.12)

For convenience, the file `config/environment.yml` describes the conda dependencies for the configuration GUI. In order to install it, run `conda env create --file config/environment.yml`. This should create a conda environment named `petrisnake-config` with the required dependencies.


The `posDemux` R package (https://github.com/yaccos/posDemux) is a keystone dependency of this workflow. This package is scheduled to be released as part of Bioconductor 3.23  For the pipeline itself, the package is available as an image Docker Hub (https://hub.docker.com/repository/docker/yaccos/posdemux/general) and is automatically used by Snakemake in the workflow provided that the argument `--software-deployment-method apptainer` is used.


When running the container, you may get the warning messages:
```
During startup - Warning messages:
1: Setting LC_CTYPE failed, using "C" 
2: Setting LC_COLLATE failed, using "C" 
3: Setting LC_TIME failed, using "C" 
4: Setting LC_MESSAGES failed, using "C" 
5: Setting LC_MONETARY failed, using "C" 
6: Setting LC_PAPER failed, using "C" 
7: Setting LC_MEASUREMENT failed, using "C"
```

These warnings can be safely ignored, but if you still want to silence them, run `export LC_ALL=C` in your shell before running Snakemake.

For using the cutoff selection tool, the `posDemux` package must be run outside the workflow. As of writing Bioconductor release 3.23 is not out yet and the package requires a minimum R version of 4.6.0 as per Bioconductor's policies. This makes the package trickier to install, but there are some options:

* Install the devel version of R which is 4.6. You may use the `rig` version manager (https://github.com/r-lib/rig) for this as it allows you to install multiple versions of R and switch between them as needed.
* Download the sources from GitHub and manually edit the minimum version requirements to your installed R version (everything above 4.4 should work).
* Pull the container image with `posDemux` with Docker and use the bundled RStudio Server instance: `docker run -e PASSWORD=<YOUR_PASS> -p 8787:8787 yaccos/posdemux:0.99.8`. Then open `http://localhost:8787/` in your web browser, type `rstudio` as the username and `<YOUR_PASS>` as the password.

For the rules inside the workflow which do not use R scripts, Snakemakes automatically handles the conda environments for these rules such that no manual intervention is required for the dependency management. These conda environment files are located in `workflow/envs`.

# Configuration



# Running the pipeline

Due to the fact that Snakemake uses Apptainer and Conda to handle the dependencies of the rules, the flag `--software-deployment-method apptainer conda` must be provided when launching the pipeline.

When the barcode threshold to use is unknown, we first run the pipeline to the end of the demultiplexer. In that case, run
`snakemake --software-deployment-method apptainer conda --cores <NUMBER_OF_CORES> determine_bc_cutoff`.

This will perform the demultiplexing of the forward reads and produce the following files in interest:

* The barcode table at `results/<sample>/<sample>_barcode_table.txt` which contains the barcode and UMI assignments to each of the reads
* The frequency table at `results/<sample>/<sample>_frequency_table.txt` which lists the frequency of each barcode combination in decending order (most common ones on top). This table is used by `posDemux::interactive_bc_cutoff()` to determine the barcode cutoff, this is how many of the top barcode combinations are kept.
* The log file of the demultiplexing at `results/<sample>/demultiplex.log`. This file both serves as a progress indicator of the process and provides an informative summary of the results in the end.

After the barcode cutoff is determined and entered into the config file, we are ready to run the rest of the pipeline:

`snakemake --software-deployment-method apptainer conda --cores <NUMBER_OF_CORES> all`

The additional files of interest produced after the entire pipeline has finished are:

* The selected frequency table at `results/<sample>/<sample>_selected_frequency_table.txt`. This is a truncated and renormalized version of the frequency table containing only the reads within the barcode cutoff.
* The figures of the Knee plot and density of barcode frequencies located at `results/<sample>/<sample>_kneePlot.pdf` and `results/<sample>/<sample>_ReadsPerBC.pdf`. The dashed vertical line indicates the selected cutoff. Note that the density plot is scaled by the number of reads on the y-axis in order to make it easier to interpret. In case you are unhappy with of the plots are made, you can recreate the plots from the frequency table using `posDemux::knee_plot()` and `posDemux::freq_plot()`.
* The gene count matrix at `results/<sample>/<sample>_gene_count_matrix.txt`. This is a tab-seperated file (remember to explicitly set `\t` as the delimiter) showing the number of detected UMIs per gene for all cells. The genes form the columns, whereas the cells form the rows.
* The UMI count table at `results/<sample>/<sample>_umi_count_table.txt` which can be seen as a linear form of the gene count matrix. Each row contains the cell barcode, the UMI, its associated gene and how many reads were found of this UMI. 
* The log file from the gene counting at `results/<sample>/count_genes.log` which serves both as a progress indicator and provides a summary of the gene count table.

This repository is already prepared with an example dataset `random20000` and a corresponding pre-populated config file, meaning that the user can try the pipeline right out of the box.