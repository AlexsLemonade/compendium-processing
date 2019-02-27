#! /bin/bash

declare -a sample_arr=(500 1000 10000 25000 50000 75000)
declare -a gene_arr=(3000 10000 15000 25000 50000 75000)

for sample_size in "${sample_arr[@]}"
do
	for gene_number in "${gene_arr[@]}"
	do
		echo "Sample size = $sample_size + Number of genes = $gene_number"
		python3 gen_gene_expression_matrix.py --num-samples $sample_size \
			--num-genes $gene_number all_log_1663_MCAR.pcl
	done
done
