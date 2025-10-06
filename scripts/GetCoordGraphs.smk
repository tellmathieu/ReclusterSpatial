import os, sys, glob
import pandas as pd
from multiprocessing import cpu_count

#######################
# Getting coordinate graphs to set coords for reclustering from Xenium Data
# By Tell Mathieu
#######################

#######################
#        Input        #
#######################

dataDir = "/Users/tell/Desktop/spatialtrans/xenium/data"
fileDetails = "/Users/tell/Desktop/spatialtrans/xenium/xenium_file_details.csv"

#######################
#     Variables       #
#######################



#######################
#      Pipeline       #
#######################

mainDir = os.getcwd()

filesAvail = glob.glob('*', root_dir = dataDir, recursive=False)

manifest = pd.read_csv(fileDetails)
tar_file_id = []
for num in range(0,len(manifest)-1):
	if manifest["tar_file_id"][num] in filesAvail:
		if pd.isnull(manifest["x_min"][num]) is False:
			tar_file_id.append(manifest["tar_file_id"][num])

tar_file_id = list(set(tar_file_id))

coordScript = os.path.join(mainDir, "scripts", "getCoordinateGraphs.R")
graph = os.path.join(dataDir, "{tar_file_id}", "{tar_file_id}_coordinate_graph_all.pdf")

rule all:
	input:
		expand(os.path.join(dataDir, "{tar_file_id}", "{tar_file_id}_coordinate_graph_all.pdf"), tar_file_id = tar_file_id)

rule reclustering:
	threads: cpu_count()
	conda: "env/getCoords.yaml"
	input:
		fileDetails = fileDetails,
	params: 
		dataDir = dataDir,
		coordScript = coordScript
	output:
		graph = graph
	shell: '''
		#run script for getting coord graphs
		Rscript {params.coordScript} \
			{params.dataDir} \
			{wildcards.tar_file_id}
			{input.fileDetails} \
			{output.graph}
	'''
