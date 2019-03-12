# Performing quality checks on a test zebrafish compendium

We performed a series of analyses on test zebrafish compendia from refine.bio to assess **quality**. 
You can find an overview figure of compendium processing [here](https://github.com/AlexsLemonade/refinebio-docs/issues/70#issuecomment-438733819).

### Docker

We use 3 separate Docker images for the analyses in this directory.

#### Comparing gene-level estimates between _D. rerio_ genome builds

For the following scripts:

```
01-build_salmon_index.sh
02-salmon_quant.sh
03-run_tximport.sh
```

We use the Docker image for processing RNA-seq data in production for refine.bio, which can be obtained and run like so:

```sh
docker pull ccdl/dr_salmon:v1.4.9
docker run -it --volume $(pwd):/media ccdl/dr_salmon:v1.4.9 bash
```

#### Identifying genes that are "differentially expressed" between technologies

We calculate gene lengths in the notebook `07-technology_diff_exp`.
This required us to roll back the version of R (to `3.4.4`) we were using to get `Rsamtools` to install properly.
The Docker image can be built and run from the top directory of this repository like so:

```sh
docker build -t <name:tag> quality_check/docker/3.4.4
docker run -it -e PASSWORD=<PASSWORD> --volume $(pwd):/home/rstudio/kitematic -p 8787:8787 <name:tag>
```
Once these have been run, you can navigate to `localhost:8787`, log in with username `rstudio` and the password you've set above, and open `07-technology_diff_exp.Rmd`.

#### All other analyses

All other notebooks use R `3.5.1`.
The Docker image can be built and run from the top directory of this repository like so:

```sh
docker build -t <name:tag> quality_check/docker/3.5.1
docker run -it -e PASSWORD=<PASSWORD> --volume $(pwd):/home/rstudio/kitematic -p 8787:8787 <name:tag>
```

The subsequent steps are the same as above.

### Summary of main findings

##### A small test compendium that is half RNA-seq samples, half microarray samples shows a 'technology bias'

<img src='https://github.com/AlexsLemonade/compendium-processing/raw/6826cc448d8bd8605ba73d30e344e7d20438234c/quality_check/plots/small_test_compendium_PCA.png' width='600'>

Here, samples are colored by whether they were measured via microarray or RNA-seq technology.

##### Genes that are often zero in RNA-seq data have lower average expression in microarray data

<img src='https://github.com/AlexsLemonade/compendium-processing/raw/6826cc448d8bd8605ba73d30e344e7d20438234c/quality_check/plots/average_expression_array_zero_count_genes.png' width='600'>

Part of any difference between technologies could be driven by lowly expressed genes (see [this notebook](https://alexslemonade.github.io/compendium-processing/quality_check/06-lowly_expressed_genes.nb.html)). We looked at the average expression values of genes in microarray samples and found that those genes that are often zero in RNA-seq samples (75th percentile) tend to have lower expression values than other genes.

##### Genes that are longer tend to have higher values in RNA-seq data as compared to microarray data

<img src='https://github.com/AlexsLemonade/compendium-processing/raw/6826cc448d8bd8605ba73d30e344e7d20438234c/quality_check/plots/larger_compendium_tstat_length_scatter.png' width='600'>

We tested for genes that are "differentially expressed" between RNA-seq samples and microarray samples (see [this notebook](https://alexslemonade.github.io/compendium-processing/quality_check/07-technology_diff_exp.nb.html)).
Genes with a positive t statistic have higher values in RNA-seq data. 
We can see that longer genes tend to have positive t statistics.
([Here](https://github.com/AlexsLemonade/compendium-processing/blob/6826cc448d8bd8605ba73d30e344e7d20438234c/quality_check/07-technology_diff_exp.Rmd#L120)'s how we calculate gene length.)

##### Shorter genes are less likely to be observed in RNA-seq data

<img src='https://github.com/AlexsLemonade/compendium-processing/raw/6826cc448d8bd8605ba73d30e344e7d20438234c/quality_check/plots/gene_length_vs_zero_count_all_genes.png' width='600'>

We [process RNA-seq data in a manner that _should_ take into account gene length](http://docs.refine.bio/en/latest/main_text.html#tximport), but if short genes have less observations this is likely to be less effective.

##### There appear to be two groups of RNA-seq samples in a larger test compendium

<img src='https://github.com/AlexsLemonade/compendium-processing/raw/6826cc448d8bd8605ba73d30e344e7d20438234c/quality_check/plots/larger_test_compendium_PCA.png' width='600'>

We can still see some difference between RNA-seq and microarray samples, but PC1 separates what appear to be two groups in the RNA-seq data.
It turns out that the samples on the left are all from the Wellcome Sanger Institute Zebrafish mutation project (see [this notebook](https://alexslemonade.github.io/compendium-processing/quality_check/11-rnaseq_bias.nb.html)).

##### Despite obvious technical biases, increasing sample size allows us to capture more biological pathways

<img src='https://github.com/AlexsLemonade/compendium-processing/raw/6826cc448d8bd8605ba73d30e344e7d20438234c/quality_check/plots/larger_compendium_plier_metrics_by_tech.png' width='600'>

Training on the whole compendium (technology=`ALL`) results in PLIER ([Mao et al. _bioRxiv_. 2017.](https://doi.org/10.1101/116061)) models with higher pathway coverage ([Taroni et al. _bioRxiv_. 2018.](https://doi.org/10.1101/395947)). 
The other results suggest that many of the major components of variance in the RNA-seq data are due to technical effects. 
(There are high numbers of latent variables, but a small proportion are associated with pathways when technology=`RNA-seq`.)