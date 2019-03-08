# Selecting an imputation method for refine.bio species compendia

### Background

[refine.bio](https://github.com/AlexsLemonade/refinebio) contains gene expression data measured with multiple technologies and many platforms. 
For instance, there are over 25 refine.bio-supported microarray platforms for human assays alone. 
Microarrays, which currently comprise the bulk of publicly available gene expression data, are limited by the sequences contained on the array itself; 
older arrays measure fewer genes compared to more recent platforms. 

The goal of refine.bio is to provide uniformly processed and normalized gene expression data. 
This includes the capability to combine samples that were measured on different technologies. 
In many cases, removing genes that are not measured for every sample will result in a greatly reduced feature space. 
To avoid this problem, it is desirable to impute missing values. 
Here, we're evaluating imputation approaches for constructing species compendia. Species compendia will be comprised of all samples from a species in our system.

### Experiment

We constructed a gene expression matrix from ~800 _Danio rerio_ samples split approximately evenly between microarray and RNA-seq data and only retained genes measured in all samples.

We masked this data using two methods (see #5):

1. Missing completely at random (MCAR) - 30% of values are randomly selected and replaced with NA
2. What we have termed "missing random rows" (MRR) which is likely more representative of what will occur in a species compendia:
    * 10% of genes are missing 10% of values
    * 10% of genes are missing 20% of values
    * 5% of genes are missing 30% of values

This process was repeated 10 times. We then processed the data five ways (see `04-process_PCL`):

* Splitting by technology (e.g., imputing microarray and RNA-seq data 
separately), performing no transformations
* Splitting by technology, `log2(x + 1)` transformation of RNA-seq data
* Imputing technologies together, no transformation
* Imputing technologies together, `log2(x + 1)` transformation of RNA-seq data
* Imputing technologies together, all samples are quantile normalized

And imputed missing values with 5 methods (`05-impute.sh`):

* [`KNNImputer` (Sleipnir)](https://libsleipnir.bitbucket.io/KNNImputer.html)
* `IterativeSVD` ([fancyimpute](https://github.com/iskandr/fancyimpute))
* `SoftImpute` (fancyimpute)
* `KNN` row orientation (fancyimpute)
* `KNN` column orientation (fancyimpute)

We calculated the RMSE between the imputed quantile normalized values and the true quantile normalized masked values.

### Results

![](https://raw.githubusercontent.com/AlexsLemonade/compendium-processing/master/select_imputation_method/plots/RMSE_all.png)

* In general, we find that imputing both technologies together with some kind of processing (e.g., log2-transformation of RNA-seq data, quantile normalization) is the best strategy.
We attribute this in part to the increase in sample size.
The rationale for the log2 transformation was to reduce the effect of outliers. 
In practice, we aim to have the distributions between samples as similar as possible so we probably want quantile normalization to be the _last_ step of any species compendium processing.

* `KNNImputer` did not perform particularly well. 
Many of the imputed values were 0 and it was unlikely to be appropriate for the number of samples that will have (see [this Bitbucket issue](https://bitbucket.org/libsleipnir/sleipnir/issues/28/segmentation-fault-with-100-000-samples)).

* `KNN` performed well and has been regularly used for imputing gene expression values ([Troyanskaya et al., 2001](https://doi.org/10.1093/bioinformatics/17.6.520)), which is why we elected to explore it further.
Unfortunately, it would likely take too long to run (see [`impute_requirements`](https://github.com/AlexsLemonade/compendium-processing/tree/master/impute_requirements))!

Ultimately, we elected to use `IterativeSVD` for initial species compendia processing.
