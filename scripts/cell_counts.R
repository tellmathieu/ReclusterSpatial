suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))

if (length(args) != 2) {
  stop("Two arguments are required", call.=FALSE)
} else if (length(args)==2) {
  data.dir = args[1]
  output = args[2]
}

xenium.obj <- LoadXenium(data.dir, fov = "fov")
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

# *****************************
#     CELL COUNTS INTO CSV
# *****************************

all_cells_counts <- GetAssayData(xenium.obj, slot = "counts")
all_cells_counts <- as.data.frame(all_cells_counts)

write.csv(all_cells_counts, output, quote=FALSE, row.names = TRUE)