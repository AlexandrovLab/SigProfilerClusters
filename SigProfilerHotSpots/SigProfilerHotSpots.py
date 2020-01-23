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


def distance_multiple_files_sims (output_path, output_path_chrom, simulation_path, interdistance, file_context2, inter, genome, sim):
	sample_path = simulation_path
	output_path_interdistance = output_path 
	output_path_chromosome_mutations = output_path_chrom

	folders = os.listdir(sample_path)
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


	if os.path.exists(output_path_interdistance):
		shutil.rmtree(output_path_interdistance)
		os.makedirs(output_path_interdistance)
	else:
		os.makedirs(output_path_interdistance)

	for file in folders:
		if file == '.DS_Store':
			continue

		sim = file.split(".")[0]
		with open(sample_path + file) as f:
			lines = [line.strip().split() for line in f]

		output = open(sample_path + file, 'w')

		for line in sorted(lines, key = lambda x: (x[15], ['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'].index(x[4]), int(x[5]))):
			print('\t'.join(line), file=output)

		output.close()


		with open(sample_path + file) as f:

			initial_line = True
			initial_mutation = True
			first_sample_occurence = True
			previous_line = None
			prev_dist = 0

			for lines in f:
				line = lines.strip().split()
				sample = line[15].split("_")
				# sim = sample[-1]
				sample = "_".join([x for x in sample[:-1]])
				chrom = line[4]
				start = int(line[5])
				ref = line[10]
				mut = line[12]
				mutType = line[9]

				if mutType == 'INDEL':
					if len(ref) > 1 and len(mut) > 1:
						continue
				else:
					if initial_line:
						init_samp = sample
						init_chrom = chrom
						init_start = start
						init_ref = ref
						init_mut = mut
						initial_line = False
						previous_line = line
						record_prev = True
						centro_start = centromeres[genome][init_chrom][0]
						centro_end = centromeres[genome][init_chrom][1]

						out_folder_dist = output_path_interdistance + init_samp + "/"
						out_folder_chrom = output_path_chromosome_mutations + init_samp + "/"
						if not os.path.exists(out_folder_dist):
							os.makedirs(out_folder_dist)
						if not os.path.exists(out_folder_chrom):
							os.makedirs(out_folder_chrom)
						if first_sample_occurence:
							out_file_intra = output_path_interdistance + init_samp + "/" + init_samp + "_" + file_context2  + "_" + sim + "_intradistance.txt"
							out_file_chrom = output_path_chromosome_mutations + init_samp + "/" + init_samp + "_" + file_context2  + "_" + sim + "_chrom.txt"
							out_int = open(out_file_intra, 'w')
							out_chr = open (out_file_chrom, 'w')
							first_sample_occurence = False
						continue
					
					else:
						if sample != init_samp:
							# print(str(prev_dist), file=out_int)
							print("\t".join([str(prev_dist), "-1"]), file=out_int)

							out_int.close()
							out_chr.close()
							init_samp = sample
							init_chrom = chrom
							init_start = start
							init_ref = ref
							init_mut = mut
							previous_line = line
							record_prev = True
							centro_start = centromeres[genome][init_chrom][0]
							centro_end = centromeres[genome][init_chrom][1]
							out_folder_dist = output_path_interdistance + init_samp + "/"
							out_folder_chrom = output_path_chromosome_mutations + init_samp + "/"
							if not os.path.exists(out_folder_dist):
								os.makedirs(out_folder_dist)

							if not os.path.exists(out_folder_chrom):
								os.makedirs(out_folder_chrom)
							out_file_intra = output_path_interdistance + init_samp + "/" + init_samp + "_" + file_context2  + "_" + sim + "_intradistance.txt"
							out_file_chrom = output_path_chromosome_mutations + init_samp + "/" + init_samp + "_" + file_context2  + "_" + sim + "_chrom.txt"
							out_int = open(out_file_intra, 'w')
							out_chr = open (out_file_chrom, 'w')
							continue
						else:
							if chrom != init_chrom:
								# print(str(prev_dist), file=out_int)
								print("\t".join([str(prev_dist), "-1"]), file=out_int)
								init_chrom = chrom
								init_start = start
								init_ref = ref
								init_mut = mut
								previous_line = line
								record_prev = True
								centro_start = centromeres[genome][init_chrom][0]
								centro_end = centromeres[genome][init_chrom][1]
								continue

							else:
								if init_start < centro_start and start > centro_start:
									# print(str(prev_dist), file=out_int)
									print("\t".join([str(prev_dist), "-1"]), file=out_int)
									initial_line = True
									initial_mutation = True
									continue
								if init_start < centro_end and start > centro_end:
									initial_line = True
									initial_mutation = True
									continue
								if init_start > centro_start and start < centro_end:
									initial_line = True
									initial_mutation = True
									continue
								if ref == init_ref and mut == init_mut and init_start == start:
									continue
								dist = start - init_start
								if initial_mutation:
									# print(str(dist), file=out_int)
									print("\t".join([str(dist), "+1"]), file=out_int)
									initial_mutation = False
								else:
									dist_min = min(dist, prev_dist)
									# print(str(dist_min), file=out_int)
									print("\t".join([str(x) for x in sorted([dist, prev_dist])]), file=out_int)

								previous_line = line
								prev_dist = dist
				init_start = start
		out_int.close()
		out_chr.close()




def distance_one_file (original_samples, output_path_original, output_path_chrom_original, file_context2, interdistance, inter, project, genome):
	sample_path = original_samples
	output_path_interdistance = output_path_original
	output_path_chromosome_mutations = output_path_chrom_original

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

	if os.path.exists(output_path_interdistance):
		shutil.rmtree(output_path_interdistance)
		os.makedirs(output_path_interdistance)
	else:
		os.makedirs(output_path_interdistance)


	context = file_context2

	#out_chr = open (out_file_chrom, 'w')
	for chromosome in centromeres[genome]:
		count = 0
		with open (sample_path + chromosome + "_" + project + ".genome") as f:
			initial_line = True
			initial_mutation = True
			first_sample_occurence = True
			previous_line = None

			prev_dist = 0
			print_file = False

			for lines in f:
				print_file=True
				line = lines.strip().split()
				sample2 = line[0]
				chrom = line[1]
				start = int(line[2])
				ref = line[3]
				mut = line[4]
				try:
					mutType = line[4]
				except:
					if len(ref) > 1 or len(mut) > 1:
						if len(ref) == 2 and len(mut) == 2:
							mutType = "DINUC"
						else:
							mutType = "INDEL"
					else:
						mutType = 'SNV' 

				if mutType == 'INDEL':
					if len(ref) > 1 and len(mut) > 1:
						continue

				first_sample_occurence = False
				if initial_line:
					sample = sample2
					out_file_intra = output_path_interdistance + sample + "_" + context + "_intradistance.txt"
					# out_file_chrom = output_path_chromosome_mutations + sample + "_" + context + "_chrom.txt"

					out_int = open (out_file_intra, 'a')
					init_chrom = chrom
					init_start = start
					init_ref = ref
					init_mut = mut
					initial_line = False
					previous_line = line
					record_prev = True
					centro_start = centromeres[genome][init_chrom][0]
					centro_end = centromeres[genome][init_chrom][1]
					continue
				
				else:
					if sample2 != sample:
						# print(str(prev_dist), file=out_int)
						print("\t".join([str(prev_dist), "-1"]), file=out_int)
						count += 1
						out_int.close()
						#out_chr.close()
						sample = sample2
						out_file_intra = output_path_interdistance + sample + "_" + context + "_intradistance.txt"
						#out_file_chrom = output_path_chromosome_mutations + sample + "_" + context + "_chrom.txt"

						out_int = open (out_file_intra, 'a')
						#out_chr = open (out_file_chrom, 'w')

						init_chrom = chrom
						init_start = start
						init_ref = ref
						init_mut = mut
						initial_line = True
						previous_line = line
						record_prev = True
						initial_mutation = True
						centro_start = centromeres[genome][init_chrom][0]
						centro_end = centromeres[genome][init_chrom][1]
						

					else:
						if chrom != init_chrom:
							# print(str(prev_dist), file=out_int)
							print("\t".join([str(prev_dist), "-1"]), file=out_int)

							count += 1
							init_chrom = chrom
							init_start = start
							init_ref = ref
							init_mut = mut
							previous_line = line
							initial_line = True
							record_prev = True
							centro_start = centromeres[genome][init_chrom][0]
							centro_end = centromeres[genome][init_chrom][1]
							initial_mutation = True
							continue

						else:
							if init_start < centro_start and start > centro_start:
								# print(str(prev_dist), file=out_int)
								print("\t".join([str(prev_dist), "-1"]), file=out_int)

								count += 1
								initial_line = True
								initial_mutation = True
								continue
							if init_start < centro_end and start > centro_end:
								initial_line = True
								initial_mutation = True
								continue
							if init_start > centro_start and start < centro_end:
								initial_line = True
								initial_mutation = True
								continue

							if ref == init_ref and mut == init_mut and init_start == start:
								continue
							dist = start - init_start
							if initial_mutation:
								# print(str(dist), file=out_int)
								print("\t".join([str(dist), "+1"]), file=out_int)

								count += 1
								initial_mutation = False
							else:
								dist_min = min(dist, prev_dist)
								# print(str(dist_min), file=out_int)
								print("\t".join([str(x) for x in sorted([dist, prev_dist])]), file=out_int)

								count += 1

							previous_line = line
							prev_dist = dist

				init_start = start

			if print_file:
				# print(str(prev_dist), file=out_int)
				print("\t".join([str(prev_dist), "-1"]), file=out_int)

				count += 1
		out_int.close()


def analysis (project, genome, contexts, simContext, input_path, output_type='all', analysis='all', interdistance='96', clustering_vaf=False):
	ncbi_chrom = {'NC_000067.6':'1', 'NC_000068.7':'2', 'NC_000069.6':'3', 'NC_000070.6':'4', 
				  'NC_000071.6':'5', 'NC_000072.6':'6', 'NC_000073.6':'7', 'NC_000074.6':'8',
				  'NC_000075.6':'9', 'NC_000076.6':'10', 'NC_000077.6':'11', 'NC_000078.6':'12',
				  'NC_000079.6':'13', 'NC_000080.6':'14', 'NC_000081.6':'15', 'NC_000082.6':'16', 
				  'NC_000083.6':'17', 'NC_000084.6':'18', 'NC_000085.6':'19', 'NC_000086.7':'X', 
				  'NC_000087.7':'Y'}


	inter = ''
	inter_label = interdistance
	interdistance = False
	ref_dir = input_path

	if analysis == 'hotspot':
		output_type = None

	if inter_label == 'SI':
		interdistance = True

	else:
		if inter_label == 'INDEL':
			inter = 'INDEL'
		else:
			inter = 'SNP'
			inter2 = 'SNV'

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
		if simContext[0] == 'INDEL':
			context_ref = 'INDEL'
			original_vcf_path = ref_dir + "references/vcf_files/" + project + "/INDEL/"

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


		vcf_files_path = os.listdir(ref_dir + "output/vcf_files/single/")

	else:
		original_vcf_path_INDEL = ref_dir + "/references/vcf_files/" + project + "/INDEL/"
		original_vcf_path_SNV = ref_dir + "/input/"
		INDEL_files = os.listdir(original_vcf_path_INDEL)
		SNV_files = os.listdir(original_vcf_path_SNV)

		if '.DS_Store' in INDEL_files:
			INDEL_files.remove('.DS_Store')
		if '.DS_Store' in SNV_files:
			SNV_files.remove('.DS_Store')
		file_name_INDEL = INDEL_files[0].split('.')
		file_extension_INDEL = file_name_INDEL[-1] 
		if file_extension_INDEL == 'genome':
			os.system("bash convert_txt_files_to_simple_files.sh " + project + " " + original_vcf_path_INDEL + " INDEL" )
		else:
			os.system("bash convert_" + file_extension_INDEL + "_files_to_simple_files.sh " + project + " " + original_vcf_path_INDEL + " INDEL")
		
		os.system("mv " + ref_dir + '/references/vcf_files/single/' + project + "_indels.genome " + ref_dir + '/references/vcf_files/single/' + project + "_indels2.genome")

		file_name_SNV = SNV_files[0].split('.')
		file_extension_SNV = file_name_SNV[-1]
		if file_extension_SNV == 'genome':
			os.system("bash convert_txt_files_to_simple_files.sh " + project + " " + original_vcf_path_SNV + " SNP")
		else:
			os.system("bash convert_" + file_extension_SNV + "_files_to_simple_files.sh " + project + " " + original_vcf_path_SNV + " SNP")
		os.system("mv " + ref_dir + '/references/vcf_files/single/' + project + "_indels.genome " + ref_dir + '/references/vcf_files/single/' + project + "_indels3.genome")


		os.system("cat " + ref_dir + '/references/vcf_files/single/'+project+"_indels3.genome " + ref_dir + '/references/vcf_files/single/'+project+"_indels2.genome >> " + ref_dir + '/references/vcf_files/single/' +project + "_indels.genome")
		os.system("rm "+ ref_dir + '/references/vcf_files/single/'+project+"_indels3.genome")
		os.system("rm " + ref_dir + '/references/vcf_files/single/'+project+"_indels2.genome")

	vcf_files_path = os.listdir(ref_dir + "output/vcf_files/single/SNV")
	if interdistance:
		file_context2 = inter_label
	else:

		file_context2 = inter_label

	output_path = ref_dir + "output/simulations/" + project + "_intradistance_" + genome + "_" + file_context2 + "/"
	output_path_chrom = ref_dir + "output/simulations/" + project + "_chromosome_mutation_count_" + genome + "_" + file_context2 + "/"
	simulation_path = ref_dir + "output/simulations/" + project + "_simulations_" + genome + "_" + file_context + "/"	
	original_samples = ref_dir + "output/vcf_files/single/SNV/"
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


	# Sorts and organizes the original input files:
	for files in os.listdir(original_samples):
		with open(original_samples + files) as f:
			lines = [line.strip().split() for line in f]
		lines = sorted(lines, key = lambda x: (x[0], int(x[2])))

		with open(original_samples + files, "w") as f:
			for line in lines:
				print("\t".join([x for x in line]), file=f)

	if output_type != None:
		print("Calculating mutational distances...", end='', flush=True)
	if output_type == 'original':
		distance_one_file (original_samples, output_path_original, output_path_chrom_original, file_context2, interdistance, inter, project, genome)
	elif output_type == 'all':
		distance_one_file (original_samples, output_path_original, output_path_chrom_original, file_context2, interdistance, inter, project, genome)
		distance_multiple_files_sims (output_path, output_path_chrom, simulation_path, interdistance, file_context2, inter, genome, file_context)
	elif output_type == 'simulations':
		distance_multiple_files_sims (output_path, output_path_chrom, simulation_path, interdistance, file_context2, inter, genome, file_context)
	else:
		pass
	if output_type != None:
		print("Completed!", flush=True)


	if analysis == 'hotspot' or analysis == 'all':
		original= True
		firstRun=True
		signature = False
		percentage = False

		hotspot.hotSpotAnalysis(project, genome, contexts, simContext, ref_dir, original, signature, percentage, firstRun, clustering_vaf)

	
	end = time.time() - start
	print("SigProfilerHotSpots successfully finished! Elapsed time: " + str(round(end, 2)) + " seconds.")
