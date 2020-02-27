#!/usr/bin/env python3
import os
import re
import argparse
from SigProfilerMatrixGenerator.scripts import convert_input_to_simple_files as convertIn
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from . import hotspot
import SigProfilerMatrixGenerator as sig
import matplotlib as plt
import sigProfilerPlotting as sigPlt
import datetime
import shutil
import platform
import numpy as np
import statistics
import scipy
import time
import sys
from SigProfilerExtractor import sigpro as sigs


def distance_multiple_files_sims (output_path, output_path_chrom, simulation_path, interdistance, file_context2, inter, genome, sim, centromeres):
	sample_path = simulation_path
	output_path_interdistance = output_path 
	output_path_chromosome_mutations = output_path_chrom

	folders = os.listdir(sample_path)


	if os.path.exists(output_path_interdistance):
		shutil.rmtree(output_path_interdistance)
		os.makedirs(output_path_interdistance)
	else:
		os.makedirs(output_path_interdistance)

	for file in os.listdir(sample_path):
		with open(sample_path + file) as f:
			lines = [line.strip().split() for line in f]
		original_lines = sorted(lines, key = lambda x: (x[15], ['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'].index(x[4]), int(x[5])))


		if len(original_lines)>1:
			distances = [[int(y[5])-int(x[5]), x[15]] + ['c'] if x[15] == y[15] and x[4] == y[4] and int(y[5]) < centromeres[genome][x[4]][0] else [int(y[5])-int(x[5]), x[15]] + ['c'] if x[15] == y[15] and x[4] == y[4] and centromeres[genome][x[4]][1] < int(x[5]) else ['a',x[15]] if x[15] != y[15] or x[4] != y[4] else [int(x[5]) - int(original_lines[original_lines.index(x)-1][5]), x[15]] + ['d'] for x,y in zip(original_lines, original_lines[1:])]			
			distances_new = distances[:] + [[distances[-1][0]] + [original_lines[-1][15]] + [distances[-1][-1]]]
			final_distances = ['bb' if x[1] != y[1] else [min(x[0],y[0])] + y[1:-1] if x[0] != 'a' and y[0] != 'a' and x[-1] != 'd' else [y[0]] + y[1:-1] if x[-1] == 'd' and x[0] != 'a' and y[0] != 'a' else [x[0],x[1]] if y[0] == 'a' and x[0] != 'a' else [y[0],y[1]] if x[0] == 'a' and y[0] != 'a' else 'bb' for x,y in zip(distances_new, distances_new[1:])]
			final_distances = [distances_new[0]] + final_distances
			final_distances_filtered = [x for x in final_distances if x[0] != 'b' and x[0] != 'a']


			sample = final_distances_filtered[0][1]
			
			file_sample = sample.split("_")[:-1]
			file_sample = ("_").join([x for x in file_sample])
			if not os.path.exists(output_path_interdistance + file_sample + "/"):
				os.makedirs(output_path_interdistance + file_sample + "/")
			out_file = open(output_path_interdistance  + file_sample + "/"+ sample + "_" + file_context2 + "_intradistance.txt", 'a')
			for x in final_distances_filtered:
				if x[1] == sample:
					print(str(x[0]), file=out_file)
				else:
					out_file.close()
					sample = x[1]
					file_sample = sample.split("_")[:-1]
					file_sample = ("_").join([x for x in file_sample])
					if sample.split("_")[0] != file_sample:
						file_sample = sample.split("_")[:-1]
						file_sample = ("_").join([x for x in file_sample])
						if not os.path.exists(output_path_interdistance + file_sample + "/"):
							os.makedirs(output_path_interdistance + file_sample + "/")
						out_file = open(output_path_interdistance  + file_sample + "/"+ sample + "_" + file_context2 + "_intradistance.txt", 'a')
						print(str(x[0]), file=out_file)
					else:
						if not os.path.exists(output_path_interdistance + file_sample + "/"):
							os.makedirs(output_path_interdistance + file_sample + "/")
						out_file = open(output_path_interdistance + file_sample + "/" + sample + "_" + file_context2 + "_intradistance.txt", 'a')
						print(str(x[0]), file=out_file)

			out_file.close()




def distance_one_file (original_samples, output_path_original, output_path_chrom_original, file_context2, interdistance, inter, project, genome, centromeres):
	sample_path = original_samples
	output_path_interdistance = output_path_original
	output_path_chromosome_mutations = output_path_chrom_original


	if os.path.exists(output_path_interdistance):
		shutil.rmtree(output_path_interdistance)
		os.makedirs(output_path_interdistance)
	else:
		os.makedirs(output_path_interdistance)
	for file in os.listdir(sample_path):
		with open(sample_path + file) as f:
			lines = [line.strip().split() for line in f]
		original_lines = sorted(lines, key = lambda x: (x[0], int(x[2])))


		if len(original_lines)>1:
			centro_start = centromeres[genome][original_lines[0][1]][0]
			centro_end = centromeres[genome][original_lines[0][1]][1]
			distances = [[int(y[2])-int(x[2])] + x + ['c'] if x[0] == y[0] and int(y[2])<centro_start else [int(y[2])-int(x[2])] + x +['c']if x[0] == y[0] and centro_end<int(x[2]) else 'aa' if x[0] != y[0] else [int(x[2]) - int(original_lines[original_lines.index(x)-1][2])] + x +['d'] for x,y in zip(original_lines, original_lines[1:])]
			distances_new = distances[:] + [[distances[-1][0]] + original_lines[-1] + [distances[-1][-1]]]
			final_distances = ['bb' if x[1] != y[1] else [min(x[0],y[0])]+y[1:-1] if x[0] != 'a' and y[0] != 'a' and x[-1] != 'd' else [y[0]] + y[1:-1] if x[-1] == 'd' and x[0] != 'a' and y[0] != 'a' else x if y[0] == 'a' and x[0] != 'a' else y if x[0] == 'a' and y[0] != 'a' else 'bb' for x,y in zip(distances_new[:], distances_new[1:])]		
			final_distances = [distances_new[0]] + final_distances
			final_distances_filtered = [x for x in final_distances if x[0] != 'b' and x[0] != 'a']



			sample = final_distances_filtered[0][1]
			out_file = open(output_path_interdistance + sample + "_" + file_context2 + "_intradistance.txt", 'a')
			for x in final_distances_filtered:
				if x[1] == sample:
					print("\t".join([str(y) for y in x]), file=out_file)
				else:
					out_file.close()
					sample = x[1]
					out_file = open(output_path_interdistance + sample + "_" + file_context2 + "_intradistance.txt", 'a')
					print("\t".join([str(y) for y in x]), file=out_file)
			out_file.close()



def analysis (project, genome, contexts, simContext, input_path, output_type='all', analysis='all', interdistance='96', clustering_vaf=False, extraction=False, startProcess=1, endProcess=25, totalIterations=1000):
	ncbi_chrom = {'NC_000067.6':'1', 'NC_000068.7':'2', 'NC_000069.6':'3', 'NC_000070.6':'4', 
				  'NC_000071.6':'5', 'NC_000072.6':'6', 'NC_000073.6':'7', 'NC_000074.6':'8',
				  'NC_000075.6':'9', 'NC_000076.6':'10', 'NC_000077.6':'11', 'NC_000078.6':'12',
				  'NC_000079.6':'13', 'NC_000080.6':'14', 'NC_000081.6':'15', 'NC_000082.6':'16', 
				  'NC_000083.6':'17', 'NC_000084.6':'18', 'NC_000085.6':'19', 'NC_000086.7':'X', 
				  'NC_000087.7':'Y'}

	centromeres = {'GRCh38':{'1': [122026460,125184587],'10': [39686683,41593521],'11':[51078349,54425074],'12':[34769408,37185252],'13': [16000001,18051248],
				'14': [16000001,18173523],'15': [17000001,19725254],'16': [36311159,38280682],'17': [22813680,26885980],'18': [15460900,20861206],
				'19': [24498981,27190874],'2': [92188146,94090557],'20': [26436233,30038348],'21': [10864561,12915808],'22': [12954789,15054318],
				'3': [90772459,93655574],'4': [49708101,51743951],'5': [46485901,50059807],'6': [58553889,59829934],'7': [58169654,60828234],
				'8': [44033745,45877265],'9': [43236168,45518558],'X': [58605580,62412542],'Y': [10316945,10544039]}, 

				'GRCh37':{'1': [121535434,124535434],'10': [39254935,42254935],'11':[51644205,54644205],'12':[34856694,37856694],'13': [16000000,19000000],
				'14': [16000000,19000000],'15': [17000000,20000000],'16': [35335801,38335801],'17': [22263006,25263006],'18': [15460898,18460898],
				'19': [24681782,27681782],'2': [92326171,95326171],'20': [26369569,29369569],'21': [11288129,14288129],'22': [13000000,16000000],
				'3': [90504854,93504854],'4': [49660117,52660117],'5': [46405641,49405641],'6': [58830166,61830166],'7': [58054331,61054331],
				'8': [43838887,46838887],'9': [47367679,50367679],'X': [58632012,61632012],'Y': [10316945,10544039]},

				'mm10':{'1': [110000,3000000],'10': [110000,3000000],'11':[110000,3000000],'12':[110000,3000000],'13': [110000,3000000],
				'14': [110000,3000000],'15': [110000,3000000],'16': [110000,3000000],'17': [110000,3000000],'18': [110000,3000000],
				'19': [110000,3000000],'2': [110000,3000000],'3': [110000,3000000],'4': [110000,3000000],'5': [110000,3000000],
				'6': [110000,3000000],'7': [110000,3000000],'8': [110000,3000000],'9': [110000,3000000],'X': [110000,3000000],'Y': [110000,3000000]},

				'mm9':{'1': [0,3000000],'10': [0,3000000],'11':[0,3000000],'12':[0,3000000],'13': [0,3000000],
				'14': [0,3000000],'15': [0,3000000],'16': [0,3000000],'17': [0,3000000],'18': [0,3000000],
				'19': [0,3000000],'2': [0,3000000],'3': [0,3000000],'4': [0,3000000],'5': [0,3000000],
				'6': [0,3000000],'7': [0,3000000],'8': [0,3000000],'9': [0,3000000],'X': [0,3000000],'Y': [0,3000000]}}

	inter = ''
	inter_label = interdistance
	interdistance = False
	ref_dir = input_path

	if analysis == 'hotspot':
		output_type = None

	if inter_label == 'SI':
		interdistance = True

	else:
		if inter_label == 'INDEL' or inter_label == 'ID':
			inter = 'ID'
		else:
			inter = 'SNP'

	simContext = sorted(simContext, reverse=True)
	file_context = "_".join(simContext)

	if ref_dir[-1] != "/":
		ref_dir += "/"
		
	output_path = ref_dir + "output/vcf_files/single/"
	if os.path.exists(output_path) or output_type != None:
		os.system("rm -r " + output_path)
		os.makedirs(output_path)
	output_log_path = ref_dir + "logs/"
	if not os.path.exists(output_log_path):
		os.makedirs(output_log_path)

	time_stamp = datetime.date.today()


	log_file = output_log_path + 'SigProfilerHotSpots_' + project + "_" + genome + "_" + str(time_stamp) + ".out"
	error_file = output_log_path + 'SigProfilerHotSpots_' + project + "_" + genome + "_" + str(time_stamp) + ".err"
	if os.path.exists(error_file):
		os.remove(error_file)
	if os.path.exists(log_file):
		 os.remove(log_file)
	sys.stderr = open(error_file, 'w')
	log_out = open(log_file, 'w')
	log_out.write("THIS FILE CONTAINS THE METADATA ABOUT SYSTEM AND RUNTIME\n\n\n")
	log_out.write("-------System Info-------\n")
	log_out.write("Operating System Name: "+ platform.uname()[0]+"\n"+"Nodename: "+ platform.uname()[1]+"\n"+"Release: "+ platform.uname()[2]+"\n"+"Version: "+ platform.uname()[3]+"\n")
	log_out.write("\n-------Python and Package Versions------- \n")
	log_out.write("Python Version: "+str(platform.sys.version_info.major)+"."+str(platform.sys.version_info.minor)+"."+str(platform.sys.version_info.micro)+"\n")
	log_out.write("SigProfilerMatrixGenerator Version: "+sig.__version__+"\n")
	log_out.write("SigProfilerPlotting version: "+sigPlt.__version__+"\n")
	log_out.write("matplotlib version: "+plt.__version__+"\n")
	log_out.write("scipy version: "+scipy.__version__+"\n")
	log_out.write("numpy version: "+np.__version__+"\n")

	log_out.write("\n-------Vital Parameters Used for the execution -------\n")
	log_out.write("Project: {}\nGenome: {}\nContext: {}\ninterdistance: {}\ninput_path: {}\noutput_type: {}\n".format(project, genome, simContext, interdistance, ref_dir, output_type))
	log_out.write("\n-------Date and Time Data------- \n")
	tic = datetime.datetime.now()
	log_out.write("Date and Clock time when the execution started: "+str(tic)+"\n\n\n")
	log_out.write("-------Runtime Checkpoints------- \n")
	log_out.close()

	print("\n\n=====================================", flush=True)
	print("Beginning SigProfilerHotSpot Analysis", flush=True)
	print("=====================================\n\n", flush=True)
	start = time.time()
	context_ref = None
	if len(simContext) == 1:
		if simContext[0] == 'INDEL' or simContext[0] == 'ID':
			context_ref = 'INDEL'

		else:
			context_ref = 'SNV' 
		original_vcf_path = ref_dir + "input/"

		vcf_files = os.listdir(original_vcf_path)
		vcf_path = original_vcf_path

		if '.DS_Store' in vcf_files:
			vcf_files.remove('.DS_Store')
		file_name = vcf_files[0].split('.')
		file_extension = file_name[-1]

		if file_extension == 'genome':
				convertIn.convertTxt(project, vcf_path, genome, output_path)
		else:
			if file_extension == 'txt':
				snv, indel, skipped, samples = convertIn.convertTxt(project, vcf_path,  genome,  output_path, ncbi_chrom, log_file)
			elif file_extension == 'vcf':
				snv, indel, skipped, samples = convertIn.convertVCF(project, vcf_path,  genome, output_path, ncbi_chrom, log_file)
			elif file_extension == 'maf':
				snv, indel, skipped, samples = convertIn.convertMAF(project, vcf_path,  genome, output_path, ncbi_chrom, log_file)
			elif file_extension == 'tsv':
				snv, indel, skipped, samples = convertIn.convertICGC(project, vcf_path,  genome, output_path, ncbi_chrom, log_file)
			else:
				print("File format not supported")

	else:
		pass
		# original_vcf_path_INDEL = ref_dir + "/references/vcf_files/" + project + "/INDEL/"
		# original_vcf_path_SNV = ref_dir + "/input/"
		# INDEL_files = os.listdir(original_vcf_path_INDEL)
		# SNV_files = os.listdir(original_vcf_path_SNV)

		# if '.DS_Store' in INDEL_files:
		# 	INDEL_files.remove('.DS_Store')
		# if '.DS_Store' in SNV_files:
		# 	SNV_files.remove('.DS_Store')
		# file_name_INDEL = INDEL_files[0].split('.')
		# file_extension_INDEL = file_name_INDEL[-1] 
		# if file_extension_INDEL == 'genome':
		# 	os.system("bash convert_txt_files_to_simple_files.sh " + project + " " + original_vcf_path_INDEL + " INDEL" )
		# else:
		# 	os.system("bash convert_" + file_extension_INDEL + "_files_to_simple_files.sh " + project + " " + original_vcf_path_INDEL + " INDEL")
		
		# os.system("mv " + ref_dir + '/references/vcf_files/single/' + project + "_indels.genome " + ref_dir + '/references/vcf_files/single/' + project + "_indels2.genome")

		# file_name_SNV = SNV_files[0].split('.')
		# file_extension_SNV = file_name_SNV[-1]
		# if file_extension_SNV == 'genome':
		# 	os.system("bash convert_txt_files_to_simple_files.sh " + project + " " + original_vcf_path_SNV + " SNP")
		# else:
		# 	os.system("bash convert_" + file_extension_SNV + "_files_to_simple_files.sh " + project + " " + original_vcf_path_SNV + " SNP")
		# os.system("mv " + ref_dir + '/references/vcf_files/single/' + project + "_indels.genome " + ref_dir + '/references/vcf_files/single/' + project + "_indels3.genome")


		# os.system("cat " + ref_dir + '/references/vcf_files/single/'+project+"_indels3.genome " + ref_dir + '/references/vcf_files/single/'+project+"_indels2.genome >> " + ref_dir + '/references/vcf_files/single/' +project + "_indels.genome")
		# os.system("rm "+ ref_dir + '/references/vcf_files/single/'+project+"_indels3.genome")
		# os.system("rm " + ref_dir + '/references/vcf_files/single/'+project+"_indels2.genome")

	file_context2 = inter_label

	output_path = ref_dir + "output/simulations/" + project + "_intradistance_" + genome + "_" + file_context2 + "/"
	output_path_chrom = ref_dir + "output/simulations/" + project + "_chromosome_mutation_count_" + genome + "_" + file_context2 + "/"
	simulation_path = ref_dir + "output/simulations/" + project + "_simulations_" + genome + "_" + file_context + "/"	
	original_samples = ref_dir + "output/vcf_files/single/" + context_ref + "/"
	output_path_original = ref_dir + "output/simulations/" + project + "_intradistance_original_" + genome + "_" + file_context2 + "/"
	output_path_chrom_original = ref_dir + "output/simulations/" + project + "_chromosome_mutation_count_original_" + genome + "_" + file_context2 + "/"
	
	# Checks for simulated data. If there are not simulations, then inform the user to generate simulations first before
	# performing the hotspot analysis:
	if not os.path.exists(simulation_path):
		print("There are no simulated data present for this project. Please generate simulations before running SigProfilerHotSpots.\n"\
				"\tThe package can be installed via pip:\n\t\t\t$ pip install SigProfilerSimulator\n"\
				"\n\tand used within a python3 sessions as follows:\n\t\t\t$ python3\n\t\t\t>> from SigProfilerSimulator import SigProfilerSimulator as sigSim\n"\
				"\t\t\t>> sigSim.SigProfilerSimulator(project, project_path, genome, contexts=['6144'], simulations=100)\n\n"\
				"\tFor a complete list of parameters, visit the github repo (https://github.com/AlexandrovLab/SigProfilerSimulator) or the documentation page (https://osf.io/usxjz/wiki/home/)")
		sys.quit()

	if len(os.listdir(simulation_path)) < 100:
		print("Please simulate a minimum of 100 simulations per sample to successfully run SigProfilerHotSpots.")
		sys.quit()


	if output_type != None:
		print("Calculating mutational distances...", end='', flush=True)
	if output_type == 'original':
		distance_one_file (original_samples, output_path_original, output_path_chrom_original, file_context2, interdistance, inter, project, genome, centromeres)
	elif output_type == 'all':
		distance_one_file (original_samples, output_path_original, output_path_chrom_original, file_context2, interdistance, inter, project, genome, centromeres)
		distance_multiple_files_sims (output_path, output_path_chrom, simulation_path, interdistance, file_context2, inter, genome, file_context, centromeres)
	elif output_type == 'simulations':
		distance_multiple_files_sims (output_path, output_path_chrom, simulation_path, interdistance, file_context2, inter, genome, file_context, centromeres)
	else:
		pass
	if output_type != None:
		print("Completed!", flush=True)


	if analysis == 'hotspot' or analysis == 'all':
		original= True
		firstRun=True
		signature = False
		percentage = False

		hotspot.hotSpotAnalysis(project, genome, contexts, simContext, ref_dir, original, signature, percentage, firstRun, clustering_vaf, centromeres)

	if extraction:
		print("Beginning signature extraction...")
		print(ref_dir+"output/extraction_clustered/")
		print(ref_dir+"output/vcf_files/"+project+"_clustered/SNV/output/SBS/"+project+"_clustered.SBS6.all")
		print(genome)
		print(startProcess)
		print(endProcess)
		print(totalIterations)
		sigs.sigProfilerExtractor("table", ref_dir+"output/extraction_clustered/", ref_dir+"output/vcf_files/"+project+"_clustered/SNV/output/SBS/"+project+"_clustered.SBS96.all", genome, startProcess=startProcess, endProcess=endProcess, totalIterations=totalIterations)#, totalIterations=totalIterations)
	# sys.stderr.close()
	end = time.time() - start
	print("SigProfilerHotSpots successfully finished! Elapsed time: " + str(round(end, 2)) + " seconds.")

