suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(spatial))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(arrow))

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 7) {
  stop("Seven arguments are required", call.=FALSE)
} else if (length(args)==7) {
  data.dir = args[1]
  file_details = args[2]
}

# setting variables for testing
data.dir <- "/Users/tell/Desktop/spatialtrans/xenium/data"
file_details <- "/Users/tell/Desktop/spatialtrans/xenium/xenium_file_details.csv"


# set working directory, based on input, and pick xenium folder
setwd(data.dir)
details <- read.csv(file_details)

all.xenium.obj <- CreateSeuratObject(counts = matrix(0, nrow = 0, ncol = 0))

details <- na.omit(details)
for (section in 1:nrow(details)) {
  print(paste0("Getting data for section: ",details$section_id[section])
 
  #setting individual line variables
  section_id <- details$section_id[section]
  file_id <- details$tar_file_id[section]
  data.dir <- file.path(wd,file_id)
  
  #getting coordinates for subsetting
  x_min <- as.numeric(details$x_min[details$section_id==section_id])
  x_max <- as.numeric(details$x_max[details$section_id==section_id])
  y_min <- as.numeric(details$y_min[details$section_id==section_id])
  y_max <- as.numeric(details$y_max[details$section_id==section_id])
  
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
  print(paste0("Loading xenium object using seurat for section: ",section_id))
  xenium.obj.orig <- LoadXenium(data.dir, fov = "fov")
  
  #keeping original for cell names - needed for creating zarr file for .xenium file
  xenium.obj <- xenium.obj.orig
  
  # remove cells with zero counts
  xenium.obj <- subset(xenium.obj.orig, subset = nCount_Xenium > 0)
  
  
  # *************************
  # SELECTING ONLY SECTION
  # *************************
  print(paste0("Starting to crop ", section_id))
  cropped.coords <- Crop(xenium.obj[["fov"]], x = c(x_min,x_max), y = c(y_min, y_max), coords = "plot")
  print(paste0("Finished cropping section: ",section_id))
  
  #DefaultBoundary(xenium.obj.section) <- "segmentation"
  # *************************
  #    CLUSTERING
  # *************************
  print(paste0("Running SCTransform: ",section_id))
  xenium.obj.section <- SCTransform(xenium.obj.section, assay = "Xenium")
  print(paste0("Running PCA: ",section_id))
  xenium.obj.section <- RunPCA(xenium.obj.section, npcs = 30, features = rownames(xenium.obj.section))
  print(paste0("Running UMAP: ",section_id))
  xenium.obj.section <- RunUMAP(xenium.obj.section, dims = 1:30)
  print(paste0("Finding neighbors: ",section_id))
  xenium.obj.section <- FindNeighbors(xenium.obj.section, reduction = "pca", dims = 1:30)
  print(paste0("Finding clusters: ",section_id))
  xenium.obj.section <- FindClusters(xenium.obj.section, resolution = 0.3)
  xenium.obj.section <- SetIdent(xenium.obj.section, value = xenium.obj.section@meta.data$seurat_clusters)
  
  all.xenium.obj <- merge(all.xenium.obj, xenium.obj.section)
}






