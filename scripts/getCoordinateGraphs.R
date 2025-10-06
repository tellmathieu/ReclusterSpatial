suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(spatial))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(arrow))

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop("Four arguments are required", call.=FALSE)
} else if (length(args)==4) {
  wd = args[1]
  file_id = args[2]
  file_details = args[3]
  outputfile = args[4]
}

# setting variables for testing
#wd <- "/Users/tell/Desktop/spatialtrans/xenium/data"
#file_id <- "0055596_wt_r_2-3_a"
#data.dir <- file.path(wd,file_id)
#file_details <- "/Users/tell/Desktop/spatialtrans/xenium/xenium_file_details.csv"
#outputfile = file.path(data.dir,"coordinate_graph_all.pdf")


# set working directory, based on input, and pick xenium folder
setwd(wd)
data.dir <- file.path(wd,file_id)
details <- read.csv(file_details)

# create transcripts.csv if it doesn't already exist. 
# Older version of Xenium output included it, but with this version, 
# we have to create it.
if(file.exists(file.path(data.dir,"transcripts.csv.gz"))){
  print("The 'transcripts.csv.gz' already exists, no need to create it.")
} else {
  transcripts <- read_parquet(file.path(data.dir, "transcripts.parquet"))
  write.csv(transcripts, gzfile(file.path(data.dir, "transcripts.csv.gz")), row.names = FALSE)
}

# loading the whole output folder from xenium
xenium.obj.orig <- LoadXenium(data.dir, fov = "fov")

xenium.obj <- xenium.obj.orig

# remove cells with zero counts
xenium.obj <- subset(xenium.obj.orig, subset = nCount_Xenium > 0)

pdf(outputfile)
ImageDimPlot(xenium.obj, axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = TRUE, nmols = 10000, flip_xy = FALSE)
dev.off()

