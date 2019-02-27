import click
import pandas as pd
import random

# adapted from Miserlou/SampleGenerator

@click.command()
@click.option("--num-samples", default=10000, help="Number of Samples")
@click.option("--num-genes", default=10000, help="Number of Genes")
@click.option("--seed", default=12345, help="Random seed for reproducibility")
@click.argument('file')
def generate_gene_expression_matrix(num_samples, num_genes, seed, file):
    """ 
    From specified file, generate a random gene expression matrix that
    is num_genes x num_samples
    """
    data = pd.read_csv(file, sep='\t', error_bad_lines=False, index_col=0)
    new_data = data.copy()

    # random set of rows (genes)
    sampled_data = new_data.sample(n=num_genes, replace=True, random_state=seed, axis=0)
    # random set of genes 
    sampled_data = sampled_data.sample(n=num_samples, replace=True, random_state=seed, axis=1)

    # write to PCL file
    output_file = file.replace(".pcl", "")
    output_file = output_file + "_genes_" + str(num_genes) + "_samples_" + str(num_samples) + ".pcl"
    output_file = "sampled_data/" + output_file
    sampled_data.to_csv(output_file, sep='\t')

if __name__ == '__main__':
    generate_gene_expression_matrix()
