#!/bin/sh

# these two lines are if you have to load the snakemake software
module load conda
mamba activate sm

snakemake --snakefile scripts/GetCoordGraphs.smk -j 10 --use-conda
