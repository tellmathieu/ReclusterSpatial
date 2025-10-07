import os, sys, glob
import pandas as pd
from multiprocessing import cpu_count

#######################
# Reclustering from Xenium Data
# By Tell Mathieu
#######################

#######################
#        Input        #
#######################

fileDetails = "/home/tmathieu/xenium/xenium_file_details.csv"
dataDir = "/home/tmathieu/xenium/data"

#######################
#     Variables       #
#######################

customClusterDirName = "gene_expression_custom_clusters"

#######################
#      Pipeline       #
#######################

mainDir = os.getcwd()

filesAvail = glob.glob('*', root_dir = dataDir, recursive=False)

manifest = pd.read_csv(fileDetails)
section_id = []
tar_file_id = []
for num in range(0,len(manifest)-1):
	if manifest["tar_file_id"][num] in filesAvail:
		if pd.isnull(manifest["x_min"][num]) is False:
			section_id.append(manifest["section_id"][num])
			tar_file_id.append(manifest["tar_file_id"][num])
print(section_id)

reclusterScript = os.path.join(mainDir, "scripts", "recluster.R")
zarrScript = os.path.join(mainDir, "scripts", "custom_zarr.py")

fileIdDir = os.path.join(dataDir, "{tar_file_id}")
customClusterDir = os.path.join(fileIdDir, "analysis", "clustering", customClusterDirName)
postAnalysisDir = os.path.join(fileIdDir ,"post_analysis", "{section_id}")
analysisPath = os.path.join(fileIdDir, "analysis")

cluster_csv_orig = os.path.join(customClusterDir, "clusters_orig_{section_id}.csv")
cluster_csv_final = os.path.join(customClusterDir, "clusters_final_{section_id}.csv")
zarr = os.path.join(fileIdDir, "custom_analysis_{section_id}.zarr.zip")
zarrText1 = zarr.split("/")[-1].split("{")[0]
zarrText2 = zarr.split("/")[-1].split("}")[1]
orig_exp = os.path.join(fileIdDir, "experiment.xenium")
experiment = os.path.join(fileIdDir, "experiment_recluster_{section_id}.xenium")

print(zarrText1)
print(zarrText2)

rule all:
	input:
		expand(os.path.join(fileIdDir, "experiment_recluster_{section_id}.xenium"), zip, tar_file_id = tar_file_id, section_id = section_id)

rule reclustering:
	threads: cpu_count()
	conda: "env/recluster.yaml"
	input:
		fileDetails = fileDetails
	params: 
		dataDir = dataDir,
		reclusterScript = reclusterScript,
		customClusterDir = customClusterDir,
		postAnalysisDir = postAnalysisDir
	output:
		cluster_csv_orig = cluster_csv_orig,
		cluster_csv_final = cluster_csv_final
	shell: '''
		mkdir -p {params.customClusterDir}
		mkdir -p {params.postAnalysisDir}

		#run script for reclustering and mfinding markers
		Rscript {params.reclusterScript} \
			{params.dataDir} \
			{wildcards.tar_file_id} \
			{input.fileDetails} \
			{wildcards.section_id} \
			{output.cluster_csv_orig} \
			{output.cluster_csv_final} \
			{params.postAnalysisDir}
	'''

rule getZarr:
	threads: cpu_count()
	conda: "env/zarr.yaml"
	input:
		cluster_csv_final = cluster_csv_final
	params:
		zarrScript = zarrScript,
		customClusterDir = customClusterDir,
		analysisPath = analysisPath
	output:
		zarr = zarr
	shell: '''
		python3 {params.zarrScript} \
			{output.zarr} \
			{params.customClusterDir} \
			{input.cluster_csv_final} \
			{params.analysisPath}
	'''

rule duplicateExperiment:
	threads: cpu_count()
	input:
		zarr = zarr,
		orig_exp = orig_exp
	params:
		zarrText1 = zarrText1,
		zarrText2 = zarrText2
	output:
		experiment = experiment
	shell: '''
		sed "s/analysis.zarr.zip/{params.zarrText1}{wildcards.section_id}{params.zarrText2}/" {input.orig_exp} > {output.experiment}
	'''

rule getCellCounts:
	threads: cpu_count()
	input:
	output:
	shell: '''

	'''


