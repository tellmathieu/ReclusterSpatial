suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(spatial))
suppressPackageStartupMessages(library(sf))

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 7) {
  stop("Seven arguments are required", call.=FALSE)
} else if (length(args)==7) {
  wd = args[1]
  file_id = args[2]
  file_details = args[3]
  section_id = args[4]
  cluster_csv_orig = args[5]
  cluster_csv_final = args[6]
  post_analysis_path = args[7]
}

# setting variables for testing
#data.dir <- "/Users/tell/Desktop/spatialtrans/xenium/data/0055468_wt-d5_r_0.2-0.8_b_g"
#file_details = "/Users/tell/Desktop/spatialtrans/xenium/xenium_file_details.csv"
#section_id = '19_wt_r_0.8b'
#cluster_folder = file.path(data.dir,"analysis","clustering","gene_expression_custom_clusters")
#cluster_csv_orig = file.path(cluster_folder,paste0("clusters_orig_", section_id, ".csv"))
#cluster_csv_final = file.path(cluster_folder,paste0("clusters_final_", section_id, ".csv"))
#post_analysis_path = "/Users/tell/Desktop/spatialtrans/xenium/data/0055468_wt-d5_r_0.2-0.8_b_g/post_analysis"

# set working directory, based on input, and pick xenium folder
setwd(wd)
data.dir <- file.path(wd,file_id)
details <- read.csv(file_details)
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
xenium.obj.orig <- LoadXenium(data.dir, fov = "fov")

xenium.obj <- xenium.obj.orig

# remove cells with zero counts
xenium.obj <- subset(xenium.obj.orig, subset = nCount_Xenium > 0)

ImageDimPlot(xenium.obj, axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = TRUE, nmols = 10000, flip_xy = FALSE)

# *************************
# SELECTING ONLY SECTION
# *************************
print("About to crop")
cropped.coords <- Crop(xenium.obj[["fov"]], x = c(x_min,x_max), y = c(y_min, y_max), coords = "plot")
print("finished crop")

xenium.obj.section <- subset(xenium.obj, cells = Cells(cropped.coords))
xenium.obj.section[["fov"]] <- cropped.coords

ImageDimPlot(xenium.obj.section, fov = "fov",axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = TRUE, nmols = 10000)

#DefaultBoundary(xenium.obj.section) <- "segmentation"
# *************************
#    CLUSTERING
# *************************
xenium.obj.section <- SCTransform(xenium.obj.section, assay = "Xenium")
xenium.obj.section <- RunPCA(xenium.obj.section, npcs = 30, features = rownames(xenium.obj.section))
xenium.obj.section <- RunUMAP(xenium.obj.section, dims = 1:30)
xenium.obj.section <- FindNeighbors(xenium.obj.section, reduction = "pca", dims = 1:30)
xenium.obj.section <- FindClusters(xenium.obj.section, resolution = 0.3)

xenium.obj.section <- SetIdent(xenium.obj.section, value = xenium.obj.section@meta.data$seurat_clusters)
# *************************
#  BASIC CLUSTER MAPS
# *************************

pdf(file.path(post_analysis_path,paste0(section_id,"_clusters.pdf")))
ImageDimPlot(xenium.obj.section, axes = TRUE, 
             coord.fixed = TRUE,size = 1.5, border.size = NULL, nmols = 10000, flip_xy = FALSE)
dev.off()

pdf(file.path(post_analysis_path,paste0(section_id,"_pca.pdf")))
PCAPlot(xenium.obj.section)
dev.off()

pdf(file.path(post_analysis_path,paste0(section_id,"_umap.pdf")))
UMAPPlot(xenium.obj.section)
dev.off()


# *************************
#    EXPORTING CLUSTERS
# *************************

print("Past Plots")
# here I'm following the tutorial from https://www.10xgenomics.com/analysis-guides/importing-customized-clustering-into-xenium-explorer
# to get a customized zarr file with our chosen clustering
cell_clusters_section <- as.data.frame(Idents(xenium.obj.section))
cell_clusters <- as.data.frame(Idents(xenium.obj.orig))
cell_clusters_all <- merge(cell_clusters,cell_clusters_section, by=0, all.x=TRUE)
cell_clusters_all$Barcode <- as.integer(rownames(cell_clusters_all)) 
cell_clusters_all$Cluster <- cell_clusters_all[,3]
cell_cluster_with_name <- cell_clusters_all[,c(1,4,5)]
cell_cluster_with_name <- rename(cell_cluster_with_name, Cell_id = "Row.names")
write.csv(cell_cluster_with_name, file.path(post_analysis_path,"clusters.csv"), quote=FALSE, row.names = FALSE)

cell_clusters_all <- cell_clusters_all[,4:5]
write.csv(cell_clusters_all, cluster_csv_orig, quote=FALSE, row.names = FALSE)

cluster_cell_char_section <- read.csv(cluster_csv_orig, header=TRUE)
num_clusters <- max(replace_na(cluster_cell_char_section$Cluster,0))
unassigned <- max(replace_na(cluster_cell_char_section$Cluster,0)) + 1
cluster_cell_char_section$Cluster <- replace_na(cluster_cell_char_section$Cluster,unassigned)
write.csv(cluster_cell_char_section, cluster_csv_final, quote=FALSE, row.names = F)

# *************************
#    GETTING MARKERS
# *************************
for (cluster in 0:as.numeric(num_clusters)) {
  marker <- FindMarkers(xenium.obj.section, ident.1 = cluster)
  marker <- marker %>% filter(marker$p_val_adj < 0.05)
  write.csv(marker, file.path(post_analysis_path, paste0("markers_",cluster, ".csv")))
}

