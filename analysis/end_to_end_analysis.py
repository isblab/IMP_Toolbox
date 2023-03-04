################ Perform end-to-end analysis for IMP models ################
############## ------>"May the Force serve u well..." <------###############
############################################################################

############# One above all #############
##-------------------------------------##
import numpy as np
import sys
import os
import subprocess

############ This is the way ############
##-------------------------------------##
# Check all the path variables before proceeding.
# Following scripts must be in the same directory --
	# run_analysis_trajectories.py
	# run_extract_models.py
	# variable_filter_v1.py
	# get_xl_viol_validation_set.py
	# contact_maps_all_pairs_surface.py
# Define the protein pairs and protein lengths for contact map script.
imp_path = "/home/kartik/imp-clean/build/setup_environment.sh"
modeling_dir_path = "../modeling/"
output_dir = "Output_"
density_file = "./density_Sin3A_hdac1.txt"

def return_major_cluster():
	file = "./model_analysis/summary_hdbscan_clustering.dat"

	with open (file, "r") as f:
		cluster_summary = f.readlines()

	models_count = {}
	#Get the index for the N_models header
	index = cluster_summary[0].split(",").index("N_models") 

	# Get cluster no. and model counts for all clusters.
	for i in range(1,len(cluster_summary)):
		line = cluster_summary[i].split(",")
		if line[0] != "-1":
			models_count[int(line[0])] = int(line[index])
		else:
			continue                              #"Jaisa chal raha h, chalne do!!!"

	clust = [(i,models_count[i]) for i in range(len(models_count)) if models_count[i] == max(models_count.values())][0]
	print("Max models in cluster: ", clust[0], " no. of models = ", clust[1])
	return clust


############# PMI Analysis ##############
##-------------------------------------##
print("\n<-----------PMI Analysis----------->")
subprocess.call( [f"{imp_path}", "python", "run_analysis_trajectories.py", f"{modeling_dir_path}", f"{output_dir}"] )
cluster = return_major_cluster()


############ Model Extraction ###########
##-------------------------------------##
print("\n<-----------Model Extraction----------->")
if int( cluster[1] ) >= 30000:	
	print(">30000 models....\n using Variable filter \n")
	var_filter_path = "./Variable_filter"
	os.mkdir( var_filter_path )

	subprocess.call( ["python", "variable_filter_v1.py", "-c", f"{cluster[0]}"] )
	subprocess.call( [f"{imp_path}", "python", "run_extract_models.py", f"{modeling_dir_path}", f"{output_dir}", f"{cluster[0]}", "True"] )

	os.rename( "ignore.KS_Test.txt", f"{var_filter_path}/ignore.KS_Test.txt" )
	os.rename( "ignore.Score_Hist_A.txt", f"{var_filter_path}/ignore.Score_Hist_A.txt" )
	os.rename( "ignore.Score_Hist_B.txt", f"{var_filter_path}/ignore.Score_Hist_B.txt" )
	os.rename( "var_filt_out.log", f"{var_filter_path}/var_filt_out.log" )
	os.rename( "var_filt_out.png", f"{var_filter_path}/var_filt_out.png" )
	os.rename( "variable_filter_v1.py", f"{var_filter_path}/variable_filter_v1.py" )
else:
	subprocess.call( [f"{imp_path}", "python", "run_extract_models.py", f"{modeling_dir_path}", f"{output_dir}", f"{cluster[0]}", "False"] )


########### Sampcon analysis ############
##-------------------------------------##
print("\n<-----------Sampcon----------->")
sampcon_path = "/home/kartik/imp-clean/imp/modules/sampcon/pyext/src/exhaust.py"
subprocess.call([
	f"{imp_path}",
	"python",
	f"{sampcon_path}",
	"-n", "sin3a",
	"-m", "cpu_omp",
	"-c", "10",
	"-d", f"{density_file}",
	"-gp",
	"-g", "5",
	"-sa", f"model_analysis/A_models_clust{cluster[0]}.txt",
	"-sb", f"model_analysis/B_models_clust{cluster[0]}.txt",
	"-ra", f"model_analysis/A_models_clust{cluster[0]}.rmf3",
	"-rb", f"model_analysis/B_models_clust{cluster[0]}.rmf3",
	"--align",
	"--prism"
	])


########### Fit to Input XLs ############
##-------------------------------------##
print("\n<-----------Fit to input XLs----------->")
xl_viol_path = "./XL_Violations"
os.mkdir( xl_viol_path )
os.rename( "get_xl_viol_validation_set.py", f"{xl_viol_path}/get_xl_viol_validation_set.py" )
os.chdir( xl_viol_path )

XLs_path = "../../Data/XL/Sin3A_Mapped_XLs_HDAC1.csv"
XL_cutoff = 35
os.system(
	f"{imp_path} python get_xl_viol_validation_set.py \
	-ia ../cluster.0.sample_A.txt \
	-ib ../cluster.0.sample_B.txt \
	-ra ../model_analysis/A_models_clust{cluster[0]}.rmf3 \
	-rb ../model_analysis/B_models_clust{cluster[0]}.rmf3 \
	-c ../cluster.0/cluster_center_model.rmf3 \
	-ta ../model_analysis/A_models_clust{cluster[0]}.txt \
	-x {XLs_path} \
	-t {XL_cutoff}"
	)

# subprocess.call([
# 	f"{imp_path}",
# 	"python",
# 	"get_xl_viol_validation_set.py",
# 	"-ia", "../cluster.0.sample_A.txt",
# 	"-ib", "../cluster.0.sample_B.txt",
# 	"-ra", f"../model_analysis/A_models_clust{cluster[0]}.rmf3",
# 	"-rb", f"../model_analysis/B_models_clust{cluster[0]}.rmf3",
# 	"-c", "../cluster.0/cluster_center_model.rmf3",
# 	"-ta", f"../model_analysis/A_models_clust{cluster[0]}.txt",
# 	"-x", f"{XLs_path}",
# 	"-t", f"{XL_cutoff}"
# 	])


############# Contact Maps ##############
##-------------------------------------##
print("\n<-----------Creating Contact maps----------->")
cmap_path = "../contact_maps"
os.mkdir( cmap_path )
os.rename( "../contact_maps_all_pairs_surface.py", f"{cmap_path}/contact_maps_all_pairs_surface.py" )
os.chdir( cmap_path )

protein1 = ["HDAC1.0", "SAP30.0", "SUDS3.0", "SIN3A.0"] 
protein2 = ["HDAC1.0", "SAP30.0", "SUDS3.0", "SIN3A.0"] 

prots = []

for index1 in range(len(protein1)):
	for index2 in range(len(protein2)):
		prot1 = protein1[index1]
		prot2 = protein2[index2]
		if prot1 == prot2:
			continue
		elif f"{prot2},{prot1}" in prots:
			continue
		else:
			prots.append(f"{prot1},{prot2}")


print( f"Contact maps will be created for the following protein pairs -- \n{prots}" )
print( f"No. of Contact maps created = {len( prots )}" )
for prot in prots:
	os.system(
		f"~/imp-clean/build/setup_environment.sh \
		python contact_maps_all_pairs_surface.py \
		-ia ../cluster.0.sample_A.txt \
		-ib ../cluster.0.sample_B.txt \
		-ra ../model_analysis/A_models_clust{cluster[0]}.rmf3 \
		-rb ../model_analysis/B_models_clust{cluster[0]}.rmf3 \
		-c ../cluster.0/cluster_center_model.rmf3 \
		-ta ../model_analysis/A_models_clust{cluster[0]}.txt \
		-p {prot} &" 
		)

################ PrISM ##################
##-------------------------------------##
print("\n<-----------Precision Analysis: PrISM----------->")
prism_path = "../prism"
prism_script = "/home/kartik/PrISM/prism/src/main.py"
os.mkdir( prism_path )
os.chdir( prism_path )

subprocess.call( [
	f"{imp_path}",
	"python", "{prism_script}",
	"--input", "../cluster.0.prism.npz" 
	"--input_type", "npz",
	"--output", "output",
	"--voxel_size", "2",
	"--return_spread",
	"--classes", "2",
	"--cores", "20",
	"--models", "1.0",
	"--n_breaks", "50",
	"--resolution", "1"
	] )
subprocess.call( ["echo", "Calculating precision completed..."] )

prism_colour_script = "/home/kartik/PrISM/prism/src/color_precision.py"
subprocess.call( 
	[f"{imp_path}",
	"python", "{prism_colour_script}",
	"--resolution", "1",
	"--annotations_file", "./output/annotations_cl2.txt",
	"--input", "../cluster.0/cluster_center_model.rmf3",
	"--output", "./patch_colored_cluster_center_model.rmf3"]
	)
os.rename( "./patch_colored_cluster_center_model.rmf3", "./output/patch_colored_cluster_center_model.rmf3" )

print( "Alas, the journey comes to an end...\n" )
print( "May the Force be with you..." )
