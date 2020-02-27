import matplotlib.pyplot as plt
import os
import numpy as np
import statistics
from statistics import stdev as stdev
from statistics import mean as mean
from scipy.stats import norm
import statsmodels.stats.multitest as sm
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import shutil
import re
import matplotlib
from collections import OrderedDict
from matplotlib.patches import Rectangle
import matplotlib.image as mpimg
from decimal import Decimal
import warnings
from statistics import median
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import time
from . import plottingFunctions

warnings.filterwarnings("ignore", message="invalid value encountered in long_scalars")
warnings.filterwarnings("ignore", message="Data has no positive values, and therefore cannot be log-scaled.")

def z_test (x, mu, sigma):
	z = (x-mu)/sigma
	p = 2*(1-norm.cdf(z))
	return(z, p)

def plot_clustered (orig_bins, sim_bins, bincenters2, sample, panel4, lower_CI, upper_CI):
	sim, = panel4.plot(bincenters2[:len(sim_bins)], sim_bins, '-', marker='o', markersize = 2, color = 'red')
	panel4.fill_between(bincenters2[:len(sim_bins)], upper_CI[:len(sim_bins)], lower_CI[:len(sim_bins)], alpha=0.5, color='red', zorder=1000)
	orig, = panel4.plot(bincenters2[:len(orig_bins)], orig_bins, '-', marker='o', markersize = 2, color = 'green')

	panel4.set_yscale('log')
	panel4.set_xscale('log')
	
def plot_non_clustered (orig_bins, sim_bins, bincenters2, sample, panel6, lower_CI, upper_CI):
	sim, = panel6.plot(bincenters2[-len(sim_bins):], sim_bins, '-', marker='o', markersize = 2, color = 'red')
	panel6.fill_between(bincenters2[-len(sim_bins):], upper_CI[-len(sim_bins):], lower_CI[-len(sim_bins):], alpha=0.5, color='red', zorder=1000)
	orig, = panel6.plot(bincenters2[-len(orig_bins):], orig_bins, '-', marker='o', markersize = 2, color = 'green')

	panel6.set_yscale('log')
	panel6.set_xscale('log')
	panel6.set_xlabel("Inter-Mutational Distance (IMD)", fontweight ='bold', fontsize=12)

def plot_hist (y2, bincenters2, q_vals, interval_line, orig_mutations, lower_CI, upper_CI, avg_bin_counts, sample, original, panel2, panel3, panel4, panel5, panel6):
	axes = plt.gca()
	distance_cutoff = bincenters2[interval_line]
	if original:
		panel2.axvline(x=distance_cutoff, linestyle='--', linewidth=0.5, color='darkred')
		extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
		panel3.text(0.08, 0.87, "Total=" + str(sum(y2)),transform=plt.gcf().transFigure)
		panel3.text(0.08, 0.575, "Clustered=" + str(orig_mutations) + "(IMD<" + str(distance_cutoff) + ")",transform=plt.gcf().transFigure)
		panel3.text(0.08, 0.56, "Expected=" + str(statistics.mean([sum(lower_CI[:interval_line]), sum(upper_CI[:interval_line])])) + "(" +str(sum(lower_CI[:interval_line])) + "-" + str(sum(upper_CI[:interval_line])) + ")",transform=plt.gcf().transFigure)
		panel5.text(0.08, .28, "Non-Clustered=" + str(sum(y2)-orig_mutations), transform=plt.gcf().transFigure)
	panel2.set_yscale('log')
	panel2.set_xscale('log')
	l = 0

	y2 = y2[l:]
	avg_bin_counts = avg_bin_counts[l:]
	bincenters2 = bincenters2[l:]
	upper_CI = upper_CI[l:]
	lower_CI = lower_CI[l:]

	sim, = panel2.plot(bincenters2, avg_bin_counts, '-', marker='o', markersize = 2, color = 'red')
	panel2.fill_between(bincenters2, upper_CI, lower_CI, alpha=0.5, color='red', zorder=1000)
	if original:
		orig, = panel2.plot(bincenters2, y2, '-', marker='o', markersize = 2, color = 'green')

	if original:
		panel2.legend([sim, orig, extra], ['simulated', 'real samples', "q_value = " + "{:.2E}".format(Decimal(q_vals))])
	else:
		panel2.legend([sim], ['simulated'])

	plot_clustered (y2[:interval_line-l+1], avg_bin_counts[:interval_line-l+1], bincenters2, sample, panel4, lower_CI, upper_CI)
	plot_non_clustered (y2[interval_line-l+1:], avg_bin_counts[interval_line-l+1:], bincenters2, sample, panel6, lower_CI, upper_CI)
	return(True)



def first_run (distances, distances_orig_all, distances_orig, original_vcf, vcf_path_clust, vcf_path_nonClust, sample, original, sim_count, matrix_path_clustered, matrix_path_nonClustered, cutoff, vcf_path_all_clust, vcf_path_all_nonClust, project, distance_path, genome, clustering_vaf, centromeres):
	maximum_sim = 0
	try:
		maximum_orig = max(distances_orig)
	except:
		maximum_orig = 0
	max_dist = max(maximum_orig, maximum_sim)


	bin1 = [0]
	i = 0
	while 2**i <= max_dist:
		bin1.append(2**i)
		i += 1
	bin1.append(max_dist)

	total_bin_counts = []
	for dist in distances:
		total_bin_counts.append(np.histogram(dist, bins = bin1)[0])

	avg_bin_counts = []
	upper_CI = []
	lower_CI = []
	mean_counts = []
	std_dev = []
	CI = int(round(.95 * sim_count, 0))

	for i in range (0, len(total_bin_counts[0]), 1):
		count = 0
		bin_std_dev = []
		for binned in total_bin_counts:
			count += binned[i]
			bin_std_dev.append(binned[i])

		bin_std_dev.sort()
		std_dev.append(stdev(bin_std_dev))
		mean_counts.append(mean(bin_std_dev))
		avg_bin_counts.append(int(round(median(bin_std_dev),0)))
		upper_CI.append(bin_std_dev[CI-1])
		lower_CI.append(bin_std_dev[-CI])


	y2, binEdges2 = np.histogram(distances_orig, bins = bin1)


	bincenters2 = binEdges2[1:-1]
	y2 = y2[1:]
	avg_bin_counts = avg_bin_counts[1:]
	mean_counts = mean_counts[1:]
	std_dev = std_dev[1:]
	upper_CI=upper_CI[1:]
	lower_CI=lower_CI[1:]


	total_mutations = 0
	sim_mutations = 0
	orig_mutations = 0

	interval_line = 0
	previous_interval = 0
	interval_percent = 0

	p_vals = []
	percent_clustered = 0.90

	for i in range (0, len(avg_bin_counts), 1):
		current_orig_mutations = y2[i]
		current_mean = mean_counts[i]
		current_stdev = std_dev[i]
		if current_stdev == 0:
			p = 0
		else:
			z, p = z_test (current_orig_mutations, current_mean, current_stdev)
		p_vals.append(p)
	q_vals = sm.fdrcorrection(p_vals)[1]

	# interval_line = 8 -> imd<256
	# interval_line=10
	# orig_mutations = sum(y2[:interval_line+1])

	for i in range (0, len(avg_bin_counts), 1):
		sim_mutations += avg_bin_counts[i]
		orig_mutations += y2[i]
		total_mutations = total_mutations + avg_bin_counts[i] + y2[i]

		current_orig_mutations = y2[i]
		current_mean = mean_counts[i]
		current_stdev = std_dev[i]
		current_q = q_vals[i]
		if orig_mutations/total_mutations < percent_clustered or current_q > 0.01:
			if abs((orig_mutations/total_mutations) - percent_clustered) < abs(previous_interval - percent_clustered):
				if current_q > 0.01:
					interval_line = i-1
					orig_mutations -= y2[i]
					interval_percent = previous_interval*100
				else:
					interval_line = i
					interval_percent = orig_mutations/total_mutations*100
			else:
				interval_line = i-1
				orig_mutations -= y2[i]
				interval_percent = previous_interval*100
			break
		previous_interval = orig_mutations/total_mutations		

	distance_cut = bincenters2[interval_line+1]

	clustered_muts = [x[1:] for x in distances_orig_all if int(x[0]) <= distance_cut]
	nonClustered_muts = [x[1:] for x in distances_orig_all if int(x[0]) > distance_cut]

	with open(vcf_path_clust + project + "_clustered.txt", 'a') as clust:
		# print("HEADER", file=clust)
		for muts in clustered_muts:
			sample = muts[0]
			chrom = muts[1]
			pos = muts[2]
			ref = muts[3]
			alt = muts[4]
			if clustering_vaf:
				vaf = line[5]
				print("\t".join([project,sample,".",genome,"SNP",chrom,pos,pos,ref,alt,"SOMATIC",vaf]),file=f)
			else:
				print("\t".join([project,sample,".",genome,"SNP",chrom,pos,pos,ref,alt,"SOMATIC"]),file=clust)

	with open(vcf_path_nonClust + project + "_nonClustered.txt", 'a') as nonclust:
		# print("HEADER", file=nonclust)
		for muts in nonClustered_muts:
			sample = muts[0]
			chrom = muts[1]
			pos = muts[2]
			ref = muts[3]
			alt = muts[4]
			print("\t".join([project,sample,".",genome,"SNP",chrom,pos,pos,ref,alt,"SOMATIC"]),file=nonclust)


	return(y2, bincenters2, q_vals[i-1], interval_line, orig_mutations, lower_CI, upper_CI, avg_bin_counts)

def hotSpotAnalysis (project, genome, contexts, simContext, ref_dir, original=False, signature=False, percentage=False, firstRun=False, clustering_vaf=False, centromeres=None):
	height = 8
	width = 13
	scaled_width = (width/1.75 *.95)/width
	scaled_height = (height/4.5)/height

	if ref_dir[-1] != "/":
		ref_dir += "/"

	simContext = sorted(simContext, reverse=True)
	simContext = "_".join(simContext)
	file_context = contexts
	if file_context == 'INDEL' or file_context == 'ID':
		cutoff = 25
	else:
		cutoff = 50

	if contexts == '96':
		matrix_file_suffix = '.SBS96.'
	elif contexts ==  'INDEL' or contexts ==  'ID':
		matrix_file_suffix = '.ID83.'


	original_vcf = ref_dir + "output/vcf_files/single/" + project + "_all.txt"
	directory = ref_dir + 'output/simulations/' + project + "_intradistance_" + genome + "_" + contexts + "/"
	directory_out = ref_dir + "output/simulations/" + project + "_simulations_" + genome + "_" + contexts + "_intradistance_plots/"
	directory_orig = ref_dir + "output/simulations/" + project + "_intradistance_original_" + genome + "_" + contexts + "/"
	distance_path = ref_dir + "output/vcf_files/" + project + "_distances/"

	if os.path.exists(distance_path):
		shutil.rmtree(distance_path)
	os.mkdir(distance_path)


	if file_context == '96':
		vcf_path_clust = ref_dir + "output/vcf_files/" + project + "_clustered/SNV/"
		vcf_path_nonClust = ref_dir + "output/vcf_files/" + project + "_nonClustered/SNV/"
		vcf_path_all_clust = ref_dir + "output/vcf_files/" + project + "_all_clustered/SNV/"
		vcf_path_all_nonClust = ref_dir + "output/vcf_files/" + project + "_all_nonClustered/SNV/"
		matrix_path = ref_dir + "output/SBS/" + project + matrix_file_suffix + "all"
		matrix_path_clustered = ref_dir + "output/vcf_files/" + project + "_clustered/SNV/output/SBS/" + project + "_clustered" + matrix_file_suffix + "all" 
		matrix_path_nonClustered = ref_dir + "output/vcf_files/" + project + "_nonClustered/SNV/output/SBS/" + project + "_nonClustered" + matrix_file_suffix + "all" 
		matrix_path_all_clustered = ref_dir + "output/vcf_files/" + project + "_all_clustered/SNV/output/SBS" + project + "_all_clustered" + matrix_file_suffix + "all"
		matrix_path_all_nonClustered = ref_dir + "/references/matrix/" + project + "_all_nonClustered/" + project + "_all_nonClustered" + matrix_file_suffix + "all"
		output_path = directory_out + project + '_intradistance_plots_' + contexts + '.pdf'

	else:
		vcf_path_clust = ref_dir + "output/vcf_files/" + project + "_clustered/INDEL/"
		vcf_path_nonClust = ref_dir + "output/vcf_files/" + project + "_nonClustered/INDEL/"
		vcf_path_all_clust = ref_dir + "output/vcf_files/" + project + "_all_clustered/INDEL/"
		vcf_path_all_nonClust = ref_dir + "output/vcf_files/" + project + "_all_nonClustered/INDEL/"	
		matrix_path = ref_dir + "output/ID/" + project + matrix_file_suffix + "all"
		matrix_path_clustered = ref_dir + "output/vcf_files/" + project + "_clustered/INDEL/output/ID/" + project + "_clustered" + matrix_file_suffix + "all" 
		matrix_path_nonClustered = ref_dir + "output/vcf_files/" + project + "_nonClustered/INDEL/output/ID/" + project + "_nonClustered" + matrix_file_suffix + "all" 
		matrix_path_all_clustered = ref_dir + "output/vcf_files/" + project + "_all_clustered/INDEL/output/ID" + project + "_all_clustered" + matrix_file_suffix + "all"
		matrix_path_all_nonClustered = ref_dir + "/references/matrix/" + project + "_all_nonClustered/" + project + "_all_nonClustered" + matrix_file_suffix + "all"
		output_path = directory_out + project + '_intradistance_plots_' + contexts + '.pdf'

	with open(vcf_path_clust + project + "_clustered.txt", 'a') as clust:
		print("HEADER", file=clust)

	with open(vcf_path_nonClust + project + "_nonClustered.txt", 'a') as nonclust:
		print("HEADER", file=nonclust)

	if os.path.exists(directory_out) == False:
		os.mkdir(directory_out)

	if firstRun:
		if os.path.exists(vcf_path_clust):
			shutil.rmtree(vcf_path_clust)

		if os.path.exists(vcf_path_nonClust):
			shutil.rmtree(vcf_path_nonClust)

	if not os.path.exists(vcf_path_clust):
		os.makedirs(vcf_path_clust)

	if not os.path.exists(vcf_path_nonClust):
		os.makedirs(vcf_path_nonClust)

	y2s = {}
	bincenters2s= {} 
	q_values = {}
	interval_lines = {}
	orig_mutations_samps = {}
	lower_CIs = {}
	upper_CIs = {}
	avg_bin_counts_samp = {}


	folders = os.listdir(directory)
	first_sort = True
	if firstRun:
		if os.path.exists(vcf_path_clust + project + "_clustered.txt"):
			os.remove(vcf_path_clust + project + "_clustered.txt")
		if os.path.exists(vcf_path_nonClust + project + "_nonClustered.txt"):
			os.remove(vcf_path_nonClust + project + "_nonClustered.txt")
		if os.path.exists(vcf_path_all_clust + project + "_all_clustered.txt"):
			os.remove(vcf_path_all_clust + project + "_all_clustered.txt")
		if os.path.exists(vcf_path_all_nonClust + project + "_all_nonClustered.txt"):
			os.remove(vcf_path_all_nonClust + project + "_all_nonClustered.txt")
	
		print("Determining sample-dependent intermutational distance (IMD) cutoff...", end='', flush=True)
		for folder in folders:
			if folder == '.DS_Store_intradistance.txt' or folder == '.DS_Store':
				continue
			sample = folder
			files = os.listdir(directory + sample + "/")
			overall_distances_all = []
			distances_orig_all = []
			distances_orig_all_samps = []



			if original:
				try:
					with open (directory_orig + sample + "_" + contexts + "_intradistance.txt") as f2:
						for lines in f2:
							line = lines.strip().split()
							if int(line[0]) >= 1:
								distances_orig_all_samps.append(line)
								distances_orig_all.append(int(line[0]))
				except:
					print(sample + " does not have nearby IDs to one another. Skipping this sample.")
					continue

			sim_count = len(files)
			for file in files:
				if file == '.DS_Store':
					continue
				distances = []
				with open(directory + sample + "/" + file) as f:
					for lines in f:
						line = lines.strip().split()
						if int(line[0]) >= 1:
							distances.append(int(line[0]))
				overall_distances_all.append(distances)
			y2s[sample], bincenters2s[sample], q_values[sample], interval_lines[sample], orig_mutations_samps[sample], lower_CIs[sample], upper_CIs[sample], avg_bin_counts_samp[sample] = first_run(overall_distances_all, distances_orig_all_samps, distances_orig_all, original_vcf, vcf_path_clust, vcf_path_nonClust, sample, original, sim_count, matrix_path_clustered, matrix_path_nonClustered, cutoff, vcf_path_all_clust, vcf_path_all_nonClust, project, distance_path, genome, clustering_vaf, centromeres)
		print("Completed!", flush=True)

		print("\nAnalyzing clustered mutations...", flush=True)
		matGen.SigProfilerMatrixGeneratorFunc(project + "_clustered", genome, vcf_path_clust,plot=False)
		print("\nAnalyzing non-clustered mutations...", flush=True)
		matGen.SigProfilerMatrixGeneratorFunc(project + "_nonClustered", genome, vcf_path_nonClust,plot=False)
		print(flush=True)

	if contexts == '96':
		with open(ref_dir + "output/vcf_files/" + project + "_clustered/SNV/output/SBS/" + project + "_clustered.SBS6.all" ) as f:
			first_line = f.readline()
			samples = first_line.strip().split()
			samples = samples[1:]
	elif contexts ==  'INDEL' or contexts ==  'ID':
		with open(ref_dir + "output/vcf_files/" + project + "_clustered/INDEL/output/ID/" + project + "_clustered.ID83.all" ) as f:
			first_line = f.readline()
			samples = first_line.strip().split()
			samples = samples[1:]

	pp = PdfPages(directory_out + project + '_intradistance_plots_' + simContext + '.pdf')

	histo = True

	print("Plotting SigProfilerHotSpot Results...", end='', flush=True)

	for folder in folders:
		if folder == '.DS_Store_intradistance.txt' or folder == '.DS_Store':
			continue
		if folder not in samples:
			histo = False
		sample = folder
		files = os.listdir(directory + sample + "/")

		fig = plt.figure(figsize = (width, height))
		panel1=plt.axes([0.075, 0.225 + scaled_height*2, scaled_width, scaled_height])
		panel2=plt.axes([0.125 + scaled_width, 0.225 + scaled_height*2, 0.3, scaled_height])
		panel3=plt.axes([0.075, 0.15 + scaled_height, scaled_width, scaled_height])
		panel4=plt.axes([0.125 + scaled_width, 0.15 + scaled_height, 0.3, scaled_height])
		panel5=plt.axes([0.075, 0.075, scaled_width, scaled_height])
		panel6=plt.axes([0.125 + scaled_width, 0.075, 0.3, scaled_height])


		if histo:
			clustered = plot_hist(y2s[sample], bincenters2s[sample], q_values[sample], interval_lines[sample], orig_mutations_samps[sample], lower_CIs[sample], upper_CIs[sample], avg_bin_counts_samp[sample], sample, original, panel2, panel3, panel4, panel5, panel6)
			if clustered:
				if file_context == '96':
					plottingFunctions.plot96_same (matrix_path, matrix_path_clustered, matrix_path_nonClustered, sample, percentage, signature, panel1, panel3, panel5, fig)
				else:
					plottingFunctions.plotINDEL_same (matrix_path, matrix_path_clustered, matrix_path_nonClustered, sample, percentage, signature, panel1, panel3, panel5, fig)
				pp.savefig()
			plt.close()
		histo = True

	plt.close()
	pp.close()
	print("Completed!", flush=True)

