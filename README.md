# Exploring [refine.bio](https://github.com/AlexsLemonade/refinebio) species compendia

Once refine.bio reaches production, we will periodically release compendia comprised of all the samples from a species that we were able to process. 
We refer to these as species compendia and envision that these collections will be useful for extracting features from a diverse set of biological conditions.
Creating a species compendia pipeline required us to tackle various problems such as selecting a method for imputing missing values.
This repository holds a series of analyses related to refine.bio species compendia divided up into related modules.
See the `README` files in the individual directories for more information.

### Modules

* [`select_imputation_method`](https://github.com/AlexsLemonade/compendium-processing/tree/master/select_imputation_method) - A series of experiments/evaluations for selecting a method for imputing missing values.
* [`human_missingness`](https://github.com/AlexsLemonade/compendium-processing/tree/master/human_missingness) - Typically genes that are measured in less than 30% of samples are removed before imputing missing values in gene expression data. How many genes would be left in the human compendium using this cutoff?
* [`impute_requirements`](https://github.com/AlexsLemonade/compendium-processing/tree/master/impute_requirements) - How long does it take to run KNN impute? (Too long for our use case.)
* [`quality_check`](https://github.com/AlexsLemonade/compendium-processing/tree/master/quality_check) - Exploring test zebrafish compendia.
