## How long does it take to run KNN impute?

**TL;DR:**
In our performance evaluations, we found that KNN imputation performed well (see [this plot](https://github.com/AlexsLemonade/compendium-processing/blob/3e78f6ed3be7647f9e1846ecbbce692118c6b422/plots/RMSE_all.png)).
We also found that the implementation we evaluated will likely take too long to run for our particular use case (species compendia).
We did this by noting the run time for gene expression matrices of varying scales.

----

This is adapted from [`Miserlou/SampleGenerator`](https://github.com/miserlou/samplegenerator).

#### Set up and Docker

```
mkdir sampled_data
docker build -t jtaroni/impute_requirements:blas docker/.
docker run -it --volume $(pwd):/home jtaroni/impute_requirements:blas /bin/bash -c 'cd .. && cd home; export LC_ALL=C.UTF-8 && export LANG=C.UTF-8; bash'
```
#### Randomly generate data 

This is generating a bunch of PCL files of specified sizes using the MCAR data included in this directory 
(`all_log_1663_MCAR.pcl`).
Note that takes a while and the files are big (up to ~60GB)!

```
./generate_sampled_data.sh
```

#### Run imputation (KNN with k=10)

```
python impute.py sampled_data/all_log_1663_MCAR_genes_<NUM_GENES>_samples_<NUM_SAMPLES>.pcl
```

A quick test of the run time for 50 samples and 50 genes can performed like so:

```
python3 gen_gene_expression_matrix.py --num-samples 50 --num-genes 50 all_log_1663_MCAR.pcl 
python3 impute.py sampled_data/all_log_1663_MCAR_genes_50_samples_50.pcl
```

_Note: Cython warnings that muddy up the logs: https://github.com/numpy/numpy/issues/11788_

### Results

All performed locally, 30% missing completely at random data sampled to be of the correct size 

#### Sample + gene "sweep"

| Sample size | Number of genes | Time to first row (s) | Approx. time per sample (row) |
| :---------: | :-------------: | :-------------------: | :---------------------------: |
|     500     |      3000       |        6.311          |           0.010               |
|    1000     |      3000       |        24.531         |           0.012               |
|    10000    |      3000       |       2455.273        |           0.029               |
|    25000    |      3000       |       15232.713       |           0.061               |
|     500     |      15000      |        30.947         |        0.053                  |
|    1000     |      15000      |      122.254          |      0.059                    |
|    10000    |      15000      |      12271.147        |       0.138                   |
|    25000    |      15000      |      75158.243        |       0.3                     |


#### Sample size held constant at 500

| Number of genes | Time to first row (s) | Approx. time per sample (row) |
| :-------------: | :-------------------: | :---------------------------: |
| 3000 | 6.962 | 0.011 |
| 6000 | 13.917 | 0.021 |
| 9000 | 20.366 | 0.031 |
| 12000 | 27.177 | 0.043 |
| 15000 | 34.242 | 0.054 |
| 18000 | 40.276 | 0.061 |


#### Number of genes held constant at 500

| Sample size | Time to first row (s) | Approx. time per sample (row) |
| :---------: | :-------------------: | :---------------------------: |
| 100 | 0.048 | 0.002 |
| 200 | 0.141 | 0.002 |
| 300 | 0.295 | 0.002 |
| 400 | 0.535 | 0.002 |
| 500 | 0.825 | 0.002 |
| 600 | 1.284 | 0.002 |
| 700 | 1.587 | 0.002 |
| 800 | 2.158 | 0.002 |
| 900 | 2.717 | 0.002 |
| 1000 | 3.371 | 0.002 |
| 5000 | 117.030 | 0.004 |
| 10000 | 472.211 | 0.005 |
| 20000 | 1742.242 | 0.011 |
| 40000 | 6735.273 | 0.012 |

When we plot the same data (with `run_time_graph.R`), we can see that the run time scales linearly with the number of genes, but looks quite a bit "worse" when looking by the number of samples.
Note that we expect many more samples than genes in the more "popular" organisms.

![](https://github.com/AlexsLemonade/compendium-processing/blob/master/impute_requirements/plots/num_samples_constant.png)

![](https://github.com/AlexsLemonade/compendium-processing/blob/master/impute_requirements/plots/num_genes_constant.png)

