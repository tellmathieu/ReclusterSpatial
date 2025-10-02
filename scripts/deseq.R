suppressPackageStartupMessages(library(tidyverse))

data_dir <- "/Users/tell/Desktop/spatialtrans/xenium/data/0055468_wt-d5_r_0.2-0.8_b_g"
transcripts <- read.csv(gzfile(file.path(data_dir,"transcripts.csv.gz"),'rt'))
cells <- read.csv(gzfile(file.path(data_dir,"cells.csv.gz"),'rt'))
cell_boundaries <- read.csv(gzfile(file.path(data_dir,"cell_boundaries.csv.gz"),'rt'))
features <- read.csv(gzfile(file.path(data_dir,"cell_feature_matrix","matrix.mtx.gz"),'rt'))


geneCountFile = 'combined_gene_count.csv'
colDataFile = 'colData.csv'
alpha = 0.05
treatVsCtrl = 'treatVsControl.csv'
resultsPath = 'results'
ctrl = "CL"