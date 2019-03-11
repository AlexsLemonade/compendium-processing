# Performing quality checks on a test zebrafish compendium

### Docker

We use 3 separate Docker images for the analyses in this directory.

#### Comparing gene-level estimates between _D. rerio_ genome builds

For the following scripts

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


