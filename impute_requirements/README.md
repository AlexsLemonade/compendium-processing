Adapted from `Miserlou/SampleGenerator`

```
# setup + docker
mkdir sampled_data
docker build -t jtaroni/impute_requirements:blas docker/.
docker run -it --volume $(pwd):/home jtaroni/impute_requirements:blas /bin/bash -c 'cd .. && cd home; export LC_ALL=C.UTF-8 && export LANG=C.UTF-8; bash'

# randomly generate data -- takes a while, files are big (up to ~60GB)!
./01-generate_sampled_data.sh

# run imputation (KNN with k=10)
python impute.py sampled_data/all_log_1663_MCAR_genes_<NUM_GENES>_samples_<NUM_SAMPLES>.pcl
```

A quick test can performed like so:
```
python3 gen_gene_expression_matrix.py --num-samples 50 --num-genes 50 all_log_1663_MCAR.pcl 
python3 impute.py sampled_data/all_log_1663_MCAR_genes_50_samples_50.pcl
```

Cython warnings that muddy up the logs: https://github.com/numpy/numpy/issues/11788