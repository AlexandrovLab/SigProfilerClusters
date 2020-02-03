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

warnings.filterwarnings("ignore", message="invalid value encountered in long_scalars")
warnings.filterwarnings("ignore", message="Data has no positive values, and therefore cannot be log-scaled.")

def z_test (x, mu, sigma):
	z = (x-mu)/sigma
	p = 2*(1-norm.cdf(z))
	return(z, p)

def plot_hist2d(orig_bins, sim_bins, bincenters2, sample, panel4, lower_CI, upper_CI):
	print(bincenters2[:len(orig_bins)])
	print(orig_bins)
	# plt.hist2d(bincenters2[:len(orig_bins)], orig_bins)
	# plt.show()

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

def plot_hist (distances, distances_orig, distances_orig_adj, original_vcf, vcf_path_clust, vcf_path_nonClust, sample, pp, original, sim_count, matrix_path_clustered, matrix_path_nonClustered, panel2, panel3, panel4, panel5, panel6, cutoff):
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
	upper_conf = []
	lower_conf = []
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

		upper_conf.append(bin_std_dev[CI-1]-int(round(median(bin_std_dev),0)))
		lower_conf.append(int(round(median(bin_std_dev),0)) - bin_std_dev[-CI])
		upper_CI.append(bin_std_dev[CI-1])
		lower_CI.append(bin_std_dev[-CI])


	y, binEdges = np.histogram(avg_bin_counts, bins = bin1)
	y2, binEdges2 = np.histogram(distances_orig, bins = bin1)

	print(distances_orig)
	print(distances_orig_adj)

	bincenters2 = binEdges2[1:-1]
	y2 = y2[1:]
	avg_bin_counts = avg_bin_counts[1:]
	mean_counts = mean_counts[1:]
	std_dev = std_dev[1:]

	upper_conf=upper_conf[1:]
	lower_conf=lower_conf[1:]
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


	axes = plt.gca()
	interval_line = 10
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
	# while y2[l] == 0:
	# 	l += 1

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
		panel2.legend([sim, orig, extra], ['simulated', 'real samples', "q_value = " + "{:.2E}".format(Decimal(q_vals[i-1]))])
	else:
		panel2.legend([sim], ['simulated'])

	if True:
		plot_clustered (y2[:interval_line-l+1], avg_bin_counts[:interval_line-l+1], bincenters2, sample, panel4, lower_CI, upper_CI)
		plot_non_clustered (y2[interval_line-l+1:], avg_bin_counts[interval_line-l+1:], bincenters2, sample, panel6, lower_CI, upper_CI)
		plot_hist2d(y2[:interval_line-l+1], avg_bin_counts[:interval_line-l+1], bincenters2, sample, panel4, lower_CI, upper_CI)
		return(True)
	else:
		return (False)

def plot_hist_final (distances, distances_orig, original_vcf, vcf_path_clust, vcf_path_nonClust, sample, pp, original, sim_count, matrix_path_clustered, matrix_path_nonClustered, panel2, panel3, panel4, panel5, panel6, cutoff):
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
	upper_conf = []
	lower_conf = []
	upper_CI = []
	lower_CI = []
	CI = int(round(.95 * sim_count, 0))

	tracker = 0
	first = True
	for l in range (0, int(len(total_bin_counts)/100), 1):
		for i in range (0, len(total_bin_counts[0]), 1):
			count = 0
			bin_std_dev = []
			for binned in total_bin_counts[tracker:tracker+100]:
				bin_std_dev.append(binned[i])

			bin_std_dev.sort()
			if first:
				avg_bin_counts.append(int(round(median(bin_std_dev),0)))
				upper_conf.append(bin_std_dev[CI-1]-int(round(median(bin_std_dev),0)))
				lower_conf.append(int(round(median(bin_std_dev),0)) - bin_std_dev[-CI])
				upper_CI.append(bin_std_dev[CI-1])
				lower_CI.append(bin_std_dev[-CI])
			else:
				try:
					avg_bin_counts[i] += int(round(median(bin_std_dev),0))
				except:
					print(bin_std_dev)
				upper_conf[i] += bin_std_dev[CI-1]-int(round(median(bin_std_dev),0))
				lower_conf[i] += int(round(median(bin_std_dev),0)) - bin_std_dev[-CI]
				upper_CI[i] += bin_std_dev[CI-1]
				lower_CI[i] += bin_std_dev[-CI]
		if l == 0:
			first = False

		tracker += 100			

	y, binEdges = np.histogram(avg_bin_counts, bins = bin1)
	y2, binEdges2 = np.histogram(distances_orig, bins = bin1)

	bincenters2 = binEdges2[1:-1]
	y2 = y2[1:]
	avg_bin_counts = avg_bin_counts[1:]

	upper_conf=upper_conf[1:]
	lower_conf=lower_conf[1:]
	upper_CI=upper_CI[1:]
	lower_CI=lower_CI[1:]


	total_mutations = 0
	sim_mutations = 0
	orig_mutations = 0

	interval_line = 0
	previous_interval = 0
	interval_percent = 0

	percent_clustered = 0.90
	for i in range (0, len(avg_bin_counts), 1):
		sim_mutations += avg_bin_counts[i]
		orig_mutations += y2[i]
		total_mutations = total_mutations + avg_bin_counts[i] + y2[i]

		current_orig_mutations = y2[i]
		current_mean = mean_counts[i]
		current_stdev = std_dev[i]
		z, p = z_test (current_orig_mutations, current_mean, current_stdev)
		if orig_mutations/total_mutations < percent_clustered or p > 0.01:

			if abs((orig_mutations/total_mutations) - percent_clustered) < abs(previous_interval - percent_clustered):
				interval_line = i
				interval_percent = orig_mutations/total_mutations*100
			else:
				interval_line = i-1
				orig_mutations -= y2[i]
				interval_percent = previous_interval*100
			break
		previous_interval = orig_mutations/total_mutations

	axes = plt.gca()
	

	if original:
		p = None
		t, p = mannwhitneyu(avg_bin_counts[:interval_line+1], y2[:interval_line+1])
		panel2.axvline(x=bincenters2[interval_line], linestyle='--', linewidth=0.5, color='darkred')
		extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)	
		panel3.text(0.08, 0.87, "Total=" + str(sum(y2)),transform=plt.gcf().transFigure)
		panel3.text(0.08, 0.575, "Clustered=" + str(orig_mutations) + "(IMD<" + str(bincenters2[interval_line]) + ")",transform=plt.gcf().transFigure)
		panel3.text(0.08, 0.56, "Expected=" + str(statistics.mean([sum(lower_CI[:interval_line+1]), sum(upper_CI[:interval_line+1])])) + "(" +str(sum(lower_CI[:interval_line+1])) + "-" + str(sum(upper_CI[:interval_line+1])) + ")",transform=plt.gcf().transFigure)
		panel5.text(0.08, .28, "Non-Clustered=" + str(sum(y2)-orig_mutations), transform=plt.gcf().transFigure)
	panel2.set_yscale('log')
	panel2.set_xscale('log')
	i = 0
	while y2[i] == 0:
		i += 1

	y2 = y2[i:]
	avg_bin_counts = avg_bin_counts[i:]
	bincenters2 = bincenters2[i:]
	upper_CI = upper_CI[i:]
	lower_CI = lower_CI[i:]

	sim, = panel2.plot(bincenters2, avg_bin_counts, '-', marker='o', markersize = 2, color = 'red')
	panel2.fill_between(bincenters2, upper_CI, lower_CI, alpha=0.5, color='red', zorder=1000)
	if original:
		orig, = panel2.plot(bincenters2, y2, '-', marker='o', markersize = 2, color = 'green')


	if original:
		panel2.legend([sim, orig, extra], ['simulated', 'real samples', "p_value = " + "{:.2E}".format(Decimal(p))])
	else:
		panel2.legend([sim], ['simulated'])

	plot_clustered (y2[:interval_line+1-i], avg_bin_counts[:interval_line+1-i], bincenters2, sample, panel4, lower_CI, upper_CI)
	plot_non_clustered (y2[interval_line+1-i:], avg_bin_counts[interval_line+1-i:], bincenters2, sample, panel6, lower_CI, upper_CI)
	return(True)



def plot96_same (matrix_path, matrix_path_clustered, matrix_path_nonClustered, sample, percentage, signature, panel1, panel3, panel5, fig):
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	mutations = dict()
	total_count = []
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		sample_index = samples.index(sample) + 1
		mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
							 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0]
			mut_type = line[0][2:5]

			if percentage:
				mutCount = float(line[sample_index])
			else:
				mutCount = int(line[sample_index])
			mutations[sample][mut_type][nuc] = mutCount

	total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
	xlabels = []

	x = 0.4
	ymax = 0
	colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
	i = 0
	for key in mutations[sample]:
		for seq in mutations[sample][key]:
			xlabels.append(seq[0]+seq[2]+seq[6])
			if signature:
				if percentage:
					panel1.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]*100 > ymax:
						ymax = mutations[sample][key][seq]*100
				else:	
					panel1.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]/total_count*100 > ymax:
						ymax = mutations[sample][key][seq]/total_count*100
			else:
				panel1.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
				if mutations[sample][key][seq] > ymax:
						ymax = mutations[sample][key][seq]
			x += 1
		i += 1

	x = .077
	y3 = .895
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 6, 1):
		panel1.add_patch(plt.Rectangle((x,y3), .088, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .09

	yText = y3 + .015
	plt.text(.11, yText, 'C>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.2, yText, 'C>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.288, yText, 'C>T', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.376, yText, 'T>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.466, yText, 'T>C', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.554, yText, 'T>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

	while y%4 != 0:
		y += 1
	ytick_offest = int(y/4)
	if signature:
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
	else:
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]		

	labs = np.arange(0.375,96.375,1)

	panel1.set_xlim([0, 96])
	panel1.set_ylim([0, y])
	panel1.set_xticks(labs)
	panel1.set_yticks(ylabs)
	count = 0
	m = 0
	for i in range (0, 96, 1):
		plt.text(i/177 + .075, .052, xlabels[i][0], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .06, xlabels[i][1], fontsize=6, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .068, xlabels[i][2], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		count += 1
		if count == 16:
			count = 0
			m += 1	


	panel1.set_yticklabels(ylabels, fontsize=8, color='b')
	panel1.set_xlabel('')
	panel1.set_ylabel('')
	panel1.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
	panel1.set_ylabel("All Mutations", fontname='Times New Roman', fontsize=15, fontweight='bold')
	panel1.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=True, labelleft=True,\
					   right=True, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=10, colors='lightgray', width=1)
	[i.set_color("black") for i in panel1.get_yticklabels()]



	#############################################################################################
	# Plot 96 bar plot (clustered)
	#############################################################################################
	mutations = dict()
	total_count = []
	with open (matrix_path_clustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		sample_index = samples.index(sample) + 1
		mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
							 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0]
			mut_type = line[0][2:5]

			if percentage:
				mutCount = float(line[sample_index])
			else:
				mutCount = int(line[sample_index])
			mutations[sample][mut_type][nuc] = mutCount

	total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
	xlabels = []

	x = 0.4
	ymax = 0
	colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
	i = 0
	for key in mutations[sample]:
		for seq in mutations[sample][key]:
			xlabels.append(seq[0]+seq[2]+seq[6])
			if signature:
				if percentage:
					panel3.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]*100 > ymax:
						ymax = mutations[sample][key][seq]*100
				else:	
					panel3.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]/total_count*100 > ymax:
						ymax = mutations[sample][key][seq]/total_count*100
			else:
				panel3.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
				if mutations[sample][key][seq] > ymax:
						ymax = mutations[sample][key][seq]
			x += 1
		i += 1

	x = .077
	y3 = .598
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 6, 1):
		panel3.add_patch(plt.Rectangle((x,y3), .088, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .09

	yText = y3 + .015
	panel3.text(.11, yText, 'C>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.2, yText, 'C>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.288, yText, 'C>T', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.376, yText, 'T>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.466, yText, 'T>C', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.554, yText, 'T>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

	panel3.set_ylabel("Clustered", fontname='Times New Roman', fontsize=15, fontweight='bold')


	while y%4 != 0:
		y += 1
	ytick_offest = int(y/4)

	if signature:
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
	else:
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]		

	labs = np.arange(0.375,96.375,1)

	panel3.set_xlim([0, 96])
	panel3.set_yticklabels(ylabels, fontsize=8, color='b')
	panel3.set_ylim([0, y])
	panel3.set_xticks(labs)
	panel3.set_yticks(ylabs)
	count = 0
	m = 0
	for i in range (0, 96, 1):
		plt.text(i/177 + .075, .3495, xlabels[i][0], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .3575, xlabels[i][1], fontsize=6, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .3655, xlabels[i][2], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		count += 1
		if count == 16:
			count = 0
			m += 1

	panel3.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=True, labelleft=True,\
					   right=True, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=10, colors='lightgray', width=1)

	panel3.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
	[i.set_color("black") for i in panel3.get_yticklabels()]


	#############################################################################################
	# Plot 96 bar plot (non-clustered)
	#############################################################################################
	mutations = dict()
	total_count = []
	with open (matrix_path_nonClustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		sample_index = samples.index(sample) + 1
		mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
							 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0]
			mut_type = line[0][2:5]

			if percentage:
				mutCount = float(line[sample_index])
			else:
				mutCount = int(line[sample_index])
			mutations[sample][mut_type][nuc] = mutCount

	total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
	xlabels = []

	x = 0.4
	ymax = 0
	colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
	i = 0
	for key in mutations[sample]:
		for seq in mutations[sample][key]:
			xlabels.append(seq[0]+seq[2]+seq[6])
			if signature:
				if percentage:
					panel5.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]*100 > ymax:
						ymax = mutations[sample][key][seq]*100
				else:	
					panel5.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]/total_count*100 > ymax:
						ymax = mutations[sample][key][seq]/total_count*100
			else:
				panel5.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
				if mutations[sample][key][seq] > ymax:
						ymax = mutations[sample][key][seq]
			x += 1
		i += 1

	x = .077
	y3 = .301
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 6, 1):
		panel5.add_patch(plt.Rectangle((x,y3), .088, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .09

	yText = y3 + .015
	panel5.text(.11, yText, 'C>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.2, yText, 'C>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.288, yText, 'C>T', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.376, yText, 'T>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.466, yText, 'T>C', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.554, yText, 'T>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

	panel5.set_ylabel("Non-Clustered", fontname='Times New Roman', fontsize=15, fontweight='bold')

	while y%4 != 0:
		y += 1
	ytick_offest = int(y/4)

	if signature:
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
	else:
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]		

	labs = np.arange(0.375,96.375,1)

	panel5.set_xlim([0, 96])
	panel5.set_yticklabels(ylabels, fontsize=8, color='b')
	panel5.set_ylim([0, y])
	panel5.set_xticks(labs)
	panel5.set_yticks(ylabs)
	panel5.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)

	count = 0
	m = 0
	for i in range (0, 96, 1):
		plt.text(i/177 + .075, .647, xlabels[i][0], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .655, xlabels[i][1], fontsize=6, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .663, xlabels[i][2], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		count += 1
		if count == 16:
			count = 0
			m += 1

	panel5.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=True, labelleft=True,\
					   right=True, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=10, colors='lightgray', width=1)

	[i.set_color("black") for i in panel5.get_yticklabels()]

	fig.suptitle(sample, fontsize=30, fontname="Times New Roman")
	plt.text(.008, .55, 'Mutation Counts', rotation='vertical', fontsize=20, fontweight='bold', fontname="Times New Roman", transform=plt.gcf().transFigure)


def plotINDEL_same (matrix_path, matrix_path_clustered, matrix_path_nonClustered, sample, percentage, signature, panel1, panel3, panel5, fig):
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	indel_types = ['1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5', '1:Del:C:6'
				   '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5', '1:Del:T:6'
				   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
				   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', 
						# >1bp INDELS
				   '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
				   '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
				   '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
				   '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
				   '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', 
				   '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', 
				   '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
				   '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
						#MicroHomology INDELS
				   '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
				   '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5', '2:Ins:M:1', 
				   '3:Ins:M:1', '3:Ins:M:2', '4:Ins:M:1', '4:Ins:M:2', '4:Ins:M:3', '5:Ins:M:1', 
				   '5:Ins:M:2', '5:Ins:M:3', '5:Ins:M:4', '5:Ins:M:5']

	mutations = dict()
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		sample_index = samples.index(sample) + 1
		mutations[sample] = {'1DelC':[0,0,0,0,0,0], '1DelT':[0,0,0,0,0,0], '1InsC':[0,0,0,0,0,0], '1InsT':[0,0,0,0,0,0], 
							 '2DelR':[0,0,0,0,0,0], '3DelR':[0,0,0,0,0,0], '4DelR':[0,0,0,0,0,0], '5DelR':[0,0,0,0,0,0],
							 '2InsR':[0,0,0,0,0,0], '3InsR':[0,0,0,0,0,0], '4InsR':[0,0,0,0,0,0], '5InsR':[0,0,0,0,0,0], 
							 '2DelM':[0], '3DelM':[0,0], '4DelM':[0,0,0], '5DelM':[0,0,0,0,0]}

		for lines in f:
			line = lines.strip().split()
			categories = line[0].split(":")
			mut_type = categories[0] + categories[1] + categories[2]
			repeat_size = int(categories[3])
			if categories[2] == 'M':
				repeat_size -= 1

			if mut_type in mutations[sample].keys():
				if percentage:
					mutCount = float(line[sample_index])
				else:
					mutCount = int(line[sample_index])
				mutations[sample][mut_type][repeat_size] = mutCount
			else:
				continue
	total_count = sum(sum(nuc) for nuc in mutations[sample].values())
	xlabels = []
	
	x = .5
	ymax = 0
	colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256], 
			  [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
			  [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
			  [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]

	i = 0
	for key in mutations[sample]:
		l = 1
		for seq in mutations[sample][key]:
			xlabels.append(l)
			if signature:
				if percentage:
					panel1.bar(x, seq*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq*100 > ymax:
						ymax = seq*100

				else:
					panel1.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq/total_count*100 > ymax:
						ymax = seq/total_count*100
			else:
				panel1.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
				if seq > ymax:
						ymax = seq
			x += 1
			l += 1
		i += 1

	x = .0757
	y_top = .895
	y_bottom = .656
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 12, 1):
		panel1.add_patch(plt.Rectangle((x,y_top), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((x,y_bottom), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .0393

	panel1.add_patch(plt.Rectangle((x,y_top), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	panel1.add_patch(plt.Rectangle((x,y_bottom), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	x +=.0058
	panel1.add_patch(plt.Rectangle((x,y_top), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	panel1.add_patch(plt.Rectangle((x,y_bottom), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	x += .0133
	panel1.add_patch(plt.Rectangle((x,y_top), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
	panel1.add_patch(plt.Rectangle((x,y_bottom), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
	x += .0203
	panel1.add_patch(plt.Rectangle((x,y_top), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
	panel1.add_patch(plt.Rectangle((x,y_bottom), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))




	yText = y_top + .0015
	plt.text(.092, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.1313, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.1706, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.2099, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.2492, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.2885, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3278, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3671, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.4064, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.4457, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.485, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5243, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.5467, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5565, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.573, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5977, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

	yText_labels_top = yText + .015
	yText_labels_bottom = y_bottom + .002
	yText_labels_bottom_sec = yText_labels_bottom - .015

	plt.text(.09, yText_labels_top, '1bp Deletion', fontsize=7, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.167, yText_labels_top, '1bp Insertion', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.266, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.42, yText_labels_top, '>1bp Insertions at Repeats\n       (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.55, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	plt.text(.079, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.156, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.275, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.43, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5475, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	x = .0767
	for i in range (0, 8, 1):
		if i != 2 and i != 3:
			if i == 1 or i == 7:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		else:
			if i == 3:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x += .03925

	for i in range (0, 4, 1):
		if i == 3:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		else:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x += .03925

	plt.text(x, yText_labels_bottom, '1', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .006
	plt.text(x, yText_labels_bottom, '1  2', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .0135
	plt.text(x, yText_labels_bottom, '1  2  3', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .02
	plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

	while y%4 != 0:
		y += 1
	if signature:
		ytick_offest = int(y/4)
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]

	else:
		ytick_offest = int(y/4)
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]


	labs = np.arange(0.375,83.375,1)
	panel1.set_xlim([0, 83])
	panel1.set_yticklabels(ylabels, fontsize=8, color='b')
	panel1.set_ylim([0, y])
	panel1.set_xticks(labs)
	panel1.set_yticks(ylabs)
	panel1.set_yticklabels(ylabels, fontsize=8)
	panel1.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
	panel1.set_xlabel('')
	panel1.set_ylabel('')
	panel1.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=False, labelleft=True,\
					   right=False, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=25, colors='gray', width=2)

	[i.set_color("black") for i in panel1.get_yticklabels()]
	panel1.set_ylabel("All Mutations", fontname='Times New Roman', fontsize=15, fontweight='bold')

	########################################################
	# Non-clustered
	########################################################
	mutations = dict()
	with open (matrix_path_clustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		sample_index = samples.index(sample) + 1
		mutations[sample] = {'1DelC':[0,0,0,0,0,0], '1DelT':[0,0,0,0,0,0], '1InsC':[0,0,0,0,0,0], '1InsT':[0,0,0,0,0,0], 
							 '2DelR':[0,0,0,0,0,0], '3DelR':[0,0,0,0,0,0], '4DelR':[0,0,0,0,0,0], '5DelR':[0,0,0,0,0,0],
							 '2InsR':[0,0,0,0,0,0], '3InsR':[0,0,0,0,0,0], '4InsR':[0,0,0,0,0,0], '5InsR':[0,0,0,0,0,0], 
							 '2DelM':[0], '3DelM':[0,0], '4DelM':[0,0,0], '5DelM':[0,0,0,0,0]}

		for lines in f:
			line = lines.strip().split()
			categories = line[0].split(":")
			mut_type = categories[0] + categories[1] + categories[2]
			repeat_size = int(categories[3])
			if categories[2] == 'M':
				repeat_size -= 1

			if mut_type in mutations[sample].keys():
				if percentage:
					mutCount = float(line[sample_index])
				else:
					mutCount = int(line[sample_index])
				mutations[sample][mut_type][repeat_size] = mutCount
			else:
				continue

	total_count = sum(sum(nuc) for nuc in mutations[sample].values())
	xlabels = []
	
	x = 0.5
	ymax = 0

	i = 0
	for key in mutations[sample]:
		l = 1
		for seq in mutations[sample][key]:
			xlabels.append(l)
			if signature:
				if percentage:
					panel3.bar(x, seq*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq*100 > ymax:
						ymax = seq*100

				else:
					panel3.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq/total_count*100 > ymax:
						ymax = seq/total_count*100
			else:
				panel3.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
				if seq > ymax:
						ymax = seq
			x += 1
			l += 1
		i += 1

	x = .0757
	y_top = .597
	y_bottom = .36
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 12, 1):
		panel1.add_patch(plt.Rectangle((x,y_top), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((x,y_bottom), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .0393

	panel3.add_patch(plt.Rectangle((x,y_top), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	panel3.add_patch(plt.Rectangle((x,y_bottom), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	x +=.0058
	panel3.add_patch(plt.Rectangle((x,y_top), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	panel3.add_patch(plt.Rectangle((x,y_bottom), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	x += .0133
	panel3.add_patch(plt.Rectangle((x,y_top), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
	panel3.add_patch(plt.Rectangle((x,y_bottom), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
	x += .0203
	panel3.add_patch(plt.Rectangle((x,y_top), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
	panel3.add_patch(plt.Rectangle((x,y_bottom), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))

	yText = y_top + .00138
	plt.text(.092, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.1313, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.1706, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.2099, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.2492, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.2885, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3278, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3671, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.4064, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.4457, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.485, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5243, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.5467, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5565, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.573, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5977, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

	yText_labels_top = yText + .015
	yText_labels_bottom = y_bottom + .002
	yText_labels_bottom_sec = yText_labels_bottom - .015

	plt.text(.09, yText_labels_top, '1bp Deletion', fontsize=7, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.167, yText_labels_top, '1bp Insertion', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.266, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.42, yText_labels_top, '>1bp Insertions at Repeats\n       (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.55, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	plt.text(.079, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.156, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.275, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.43, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5475, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	x = .0767
	for i in range (0, 8, 1):
		if i != 2 and i != 3:
			if i == 1 or i == 7:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		else:
			if i == 3:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x += .03925

	for i in range (0, 4, 1):
		if i == 3:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		else:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .03925
	plt.text(x, yText_labels_bottom, '1', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .006
	plt.text(x, yText_labels_bottom, '1  2', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .0135
	plt.text(x, yText_labels_bottom, '1  2  3', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .02
	plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	
	while y%4 != 0:
		y += 1
	if signature:
		ytick_offest = int(y/4)
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]

	else:
		ytick_offest = int(y/4)
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]

	panel3.set_xlim([0, 83])
	panel3.set_yticklabels(ylabels, fontsize=8, color='b')
	panel3.set_ylim([0, y])
	panel3.set_xticks(labs)
	panel3.set_yticks(ylabs)
	panel3.set_yticklabels(ylabels, fontsize=8)
	panel3.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
	panel3.set_xlabel('')
	panel3.set_ylabel('')
	panel3.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=False, labelleft=True,\
					   right=False, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=25, colors='gray', width=2)

	[i.set_color("black") for i in panel3.get_yticklabels()]
	panel3.set_ylabel("Clustered", fontname='Times New Roman', fontsize=15, fontweight='bold')

	########################################################
	# Non-clustered
	########################################################
	mutations = dict()
	with open (matrix_path_nonClustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		try:
			sample_index = samples.index(sample) + 1
		except:
			pass
		mutations[sample] = {'1DelC':[0,0,0,0,0,0], '1DelT':[0,0,0,0,0,0], '1InsC':[0,0,0,0,0,0], '1InsT':[0,0,0,0,0,0], 
							 '2DelR':[0,0,0,0,0,0], '3DelR':[0,0,0,0,0,0], '4DelR':[0,0,0,0,0,0], '5DelR':[0,0,0,0,0,0],
							 '2InsR':[0,0,0,0,0,0], '3InsR':[0,0,0,0,0,0], '4InsR':[0,0,0,0,0,0], '5InsR':[0,0,0,0,0,0], 
							 '2DelM':[0], '3DelM':[0,0], '4DelM':[0,0,0], '5DelM':[0,0,0,0,0]}

		for lines in f:
			line = lines.strip().split()
			categories = line[0].split(":")
			mut_type = categories[0] + categories[1] + categories[2]
			repeat_size = int(categories[3])
			if categories[2] == 'M':
				repeat_size -= 1

			if mut_type in mutations[sample].keys():
				if percentage:
					mutCount = float(line[sample_index])
				else:
					mutCount = int(line[sample_index])
				mutations[sample][mut_type][repeat_size] = mutCount
			else:
				continue

	total_count = sum(sum(nuc) for nuc in mutations[sample].values())
	xlabels = []
	x = 0.5
	ymax = 0
	i = 0
	for key in mutations[sample]:
		l = 1
		for seq in mutations[sample][key]:
			xlabels.append(l)
			if signature:
				if percentage:
					panel5.bar(x, seq*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq*100 > ymax:
						ymax = seq*100

				else:
					panel5.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq/total_count*100 > ymax:
						ymax = seq/total_count*100
			else:
				panel5.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
				if seq > ymax:
						ymax = seq
			x += 1
			l += 1
		i += 1

	x = .0757
	y_top = .3
	y_bottom = .0625
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 12, 1):
		panel5.add_patch(plt.Rectangle((x,y_top), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		panel5.add_patch(plt.Rectangle((x,y_bottom), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .0393

	panel5.add_patch(plt.Rectangle((x,y_top), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	panel5.add_patch(plt.Rectangle((x,y_bottom), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	x +=.0058
	panel5.add_patch(plt.Rectangle((x,y_top), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	panel5.add_patch(plt.Rectangle((x,y_bottom), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	x += .0133
	panel5.add_patch(plt.Rectangle((x,y_top), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
	panel5.add_patch(plt.Rectangle((x,y_bottom), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
	x += .0203
	panel5.add_patch(plt.Rectangle((x,y_top), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
	panel5.add_patch(plt.Rectangle((x,y_bottom), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))

	yText = y_top + .00138
	plt.text(.092, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.1313, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.1706, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.2099, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.2492, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.2885, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3278, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3671, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.4064, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.4457, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.485, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5243, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.5467, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5565, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.573, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5977, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

	yText_labels_top = yText + .015
	yText_labels_bottom = y_bottom + .002
	yText_labels_bottom_sec = yText_labels_bottom - .015

	plt.text(.09, yText_labels_top, '1bp Deletion', fontsize=7, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.167, yText_labels_top, '1bp Insertion', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.266, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.42, yText_labels_top, '>1bp Insertions at Repeats\n       (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.55, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	plt.text(.079, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.156, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.275, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.43, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5475, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	x = .0767
	for i in range (0, 8, 1):
		if i != 2 and i != 3:
			if i == 1 or i == 7:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		else:
			if i == 3:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x += .03925

	for i in range (0, 4, 1):
		if i == 3:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		else:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .03925
	plt.text(x, yText_labels_bottom, '1', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .006
	plt.text(x, yText_labels_bottom, '1  2', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .0135
	plt.text(x, yText_labels_bottom, '1  2  3', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .02
	plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	
	while y%4 != 0:
		y += 1
	if signature:
		ytick_offest = int(y/4)
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]

	else:
		ytick_offest = int(y/4)
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]


	panel5.set_xlim([0, 83])
	panel5.set_yticklabels(ylabels, fontsize=8, color='b')
	panel5.set_ylim([0, y])
	panel5.set_xticks(labs)
	panel5.set_yticks(ylabs)
	panel5.set_yticklabels(ylabels, fontsize=8)
	panel5.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)

	fig.suptitle(sample, fontsize=30, fontname="Times New Roman")
	plt.text(.008, .55, 'Mutation Counts', rotation='vertical', fontsize=20, fontweight='bold', fontname="Times New Roman", transform=plt.gcf().transFigure)
	panel5.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=False, labelleft=True,\
					   right=False, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=25, colors='gray', width=2)
	panel5.set_ylabel("Non-Clustered", fontname='Times New Roman', fontsize=15, fontweight='bold')
	[i.set_color("black") for i in panel5.get_yticklabels()]

def plotINDEL_same_final (matrix_path, matrix_path_clustered, matrix_path_nonClustered, sample, percentage, signature, panel1, panel3, panel5, fig):
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	indel_types = ['1:Del:C:1', '1:Del:C:2', '1:Del:C:3', '1:Del:C:4', '1:Del:C:5', '1:Del:C:6'
				   '1:Del:T:1', '1:Del:T:2', '1:Del:T:3', '1:Del:T:4', '1:Del:T:5', '1:Del:T:6'
				   '1:Ins:C:0', '1:Ins:C:1', '1:Ins:C:2', '1:Ins:C:3', '1:Ins:C:4', '1:Ins:C:5',
				   '1:Ins:T:0', '1:Ins:T:1', '1:Ins:T:2', '1:Ins:T:3', '1:Ins:T:4', '1:Ins:T:5', 
						# >1bp INDELS
				   '2:Del:R:0', '2:Del:R:1', '2:Del:R:2', '2:Del:R:3', '2:Del:R:4', '2:Del:R:5',
				   '3:Del:R:0', '3:Del:R:1', '3:Del:R:2', '3:Del:R:3', '3:Del:R:4', '3:Del:R:5',
				   '4:Del:R:0', '4:Del:R:1', '4:Del:R:2', '4:Del:R:3', '4:Del:R:4', '4:Del:R:5',
				   '5:Del:R:0', '5:Del:R:1', '5:Del:R:2', '5:Del:R:3', '5:Del:R:4', '5:Del:R:5',
				   '2:Ins:R:0', '2:Ins:R:1', '2:Ins:R:2', '2:Ins:R:3', '2:Ins:R:4', '2:Ins:R:5', 
				   '3:Ins:R:0', '3:Ins:R:1', '3:Ins:R:2', '3:Ins:R:3', '3:Ins:R:4', '3:Ins:R:5', 
				   '4:Ins:R:0', '4:Ins:R:1', '4:Ins:R:2', '4:Ins:R:3', '4:Ins:R:4', '4:Ins:R:5',
				   '5:Ins:R:0', '5:Ins:R:1', '5:Ins:R:2', '5:Ins:R:3', '5:Ins:R:4', '5:Ins:R:5',
						#MicroHomology INDELS
				   '2:Del:M:1', '3:Del:M:1', '3:Del:M:2', '4:Del:M:1', '4:Del:M:2', '4:Del:M:3',
				   '5:Del:M:1', '5:Del:M:2', '5:Del:M:3', '5:Del:M:4', '5:Del:M:5', '2:Ins:M:1', 
				   '3:Ins:M:1', '3:Ins:M:2', '4:Ins:M:1', '4:Ins:M:2', '4:Ins:M:3', '5:Ins:M:1', 
				   '5:Ins:M:2', '5:Ins:M:3', '5:Ins:M:4', '5:Ins:M:5']

	mutations = dict()
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		sample = 'all'
		mutations[sample] = {'1DelC':[0,0,0,0,0,0], '1DelT':[0,0,0,0,0,0], '1InsC':[0,0,0,0,0,0], '1InsT':[0,0,0,0,0,0], 
							 '2DelR':[0,0,0,0,0,0], '3DelR':[0,0,0,0,0,0], '4DelR':[0,0,0,0,0,0], '5DelR':[0,0,0,0,0,0],
							 '2InsR':[0,0,0,0,0,0], '3InsR':[0,0,0,0,0,0], '4InsR':[0,0,0,0,0,0], '5InsR':[0,0,0,0,0,0], 
							 '2DelM':[0], '3DelM':[0,0], '4DelM':[0,0,0], '5DelM':[0,0,0,0,0]}

		for lines in f:
			line = lines.strip().split()
			categories = line[0].split(":")
			mut_type = categories[0] + categories[1] + categories[2]
			repeat_size = int(categories[3])
			if categories[2] == 'M':
				repeat_size -= 1

			if mut_type in mutations[sample].keys():
				if percentage:
					mutCount = float(sum([int(x) for x in line[1:]]))
				else:
					mutCount = sum([int(x) for x in line[1:]])
				mutations[sample][mut_type][repeat_size] += mutCount
			else:
				continue

	total_count = sum(sum(nuc) for nuc in mutations[sample].values())
	xlabels = []
	
	x = .5
	ymax = 0
	colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256], 
			  [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
			  [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
			  [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]

	i = 0
	for key in mutations[sample]:
		l = 1
		for seq in mutations[sample][key]:
			xlabels.append(l)
			if signature:
				if percentage:
					panel1.bar(x, seq*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq*100 > ymax:
						ymax = seq*100

				else:
					panel1.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq/total_count*100 > ymax:
						ymax = seq/total_count*100
			else:
				panel1.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
				if seq > ymax:
						ymax = seq
			x += 1
			l += 1
		i += 1

	x = .0757
	y_top = .895
	y_bottom = .656
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 12, 1):
		panel1.add_patch(plt.Rectangle((x,y_top), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((x,y_bottom), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .0393

	panel1.add_patch(plt.Rectangle((x,y_top), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	panel1.add_patch(plt.Rectangle((x,y_bottom), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	x +=.0058
	panel1.add_patch(plt.Rectangle((x,y_top), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	panel1.add_patch(plt.Rectangle((x,y_bottom), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	x += .0133
	panel1.add_patch(plt.Rectangle((x,y_top), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
	panel1.add_patch(plt.Rectangle((x,y_bottom), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
	x += .0203
	panel1.add_patch(plt.Rectangle((x,y_top), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
	panel1.add_patch(plt.Rectangle((x,y_bottom), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))




	yText = y_top + .0015
	plt.text(.092, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.1313, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.1706, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.2099, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.2492, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.2885, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3278, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3671, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.4064, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.4457, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.485, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5243, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.5467, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5565, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.573, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5977, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

	yText_labels_top = yText + .015
	yText_labels_bottom = y_bottom + .002
	yText_labels_bottom_sec = yText_labels_bottom - .015

	plt.text(.09, yText_labels_top, '1bp Deletion', fontsize=7, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.167, yText_labels_top, '1bp Insertion', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.266, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.42, yText_labels_top, '>1bp Insertions at Repeats\n       (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.55, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	plt.text(.079, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.156, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.275, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.43, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5475, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	x = .0767
	for i in range (0, 8, 1):
		if i != 2 and i != 3:
			if i == 1 or i == 7:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		else:
			if i == 3:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x += .03925

	for i in range (0, 4, 1):
		if i == 3:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		else:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x += .03925

	plt.text(x, yText_labels_bottom, '1', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .006
	plt.text(x, yText_labels_bottom, '1  2', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .0135
	plt.text(x, yText_labels_bottom, '1  2  3', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .02
	plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

	while y%4 != 0:
		y += 1
	if signature:
		ytick_offest = int(y/4)
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]

	else:
		ytick_offest = int(y/4)
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]


	labs = np.arange(0.375,83.375,1)
	panel1.set_xlim([0, 83])
	panel1.set_yticklabels(ylabels, fontsize=8, color='b')
	panel1.set_ylim([0, y])
	panel1.set_xticks(labs)
	panel1.set_yticks(ylabs)
	panel1.set_yticklabels(ylabels, fontsize=8)
	panel1.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
	panel1.set_xlabel('')
	panel1.set_ylabel('')
	panel1.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=False, labelleft=True,\
					   right=False, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=25, colors='gray', width=2)

	[i.set_color("black") for i in panel1.get_yticklabels()]
	panel1.set_ylabel("All Mutations", fontname='Times New Roman', fontsize=15, fontweight='bold')

	########################################################
	# Non-clustered
	########################################################
	mutations = dict()
	with open (matrix_path_clustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		mutations[sample] = {'1DelC':[0,0,0,0,0,0], '1DelT':[0,0,0,0,0,0], '1InsC':[0,0,0,0,0,0], '1InsT':[0,0,0,0,0,0], 
							 '2DelR':[0,0,0,0,0,0], '3DelR':[0,0,0,0,0,0], '4DelR':[0,0,0,0,0,0], '5DelR':[0,0,0,0,0,0],
							 '2InsR':[0,0,0,0,0,0], '3InsR':[0,0,0,0,0,0], '4InsR':[0,0,0,0,0,0], '5InsR':[0,0,0,0,0,0], 
							 '2DelM':[0], '3DelM':[0,0], '4DelM':[0,0,0], '5DelM':[0,0,0,0,0]}

		for lines in f:
			line = lines.strip().split()
			categories = line[0].split(":")
			mut_type = categories[0] + categories[1] + categories[2]
			repeat_size = int(categories[3])
			if categories[2] == 'M':
				repeat_size -= 1

			if mut_type in mutations[sample].keys():
				if percentage:
					mutCount = float(sum([int(x) for x in line[1:]]))
				else:
					mutCount = sum([int(x) for x in line[1:]])
				mutations[sample][mut_type][repeat_size] = mutCount
			else:
				continue

	total_count = sum(sum(nuc) for nuc in mutations[sample].values())
	xlabels = []
	
	x = 0.5
	ymax = 0

	i = 0
	for key in mutations[sample]:
		l = 1
		for seq in mutations[sample][key]:
			xlabels.append(l)
			if signature:
				if percentage:
					panel3.bar(x, seq*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq*100 > ymax:
						ymax = seq*100

				else:
					panel3.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq/total_count*100 > ymax:
						ymax = seq/total_count*100
			else:
				panel3.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
				if seq > ymax:
						ymax = seq
			x += 1
			l += 1
		i += 1

	x = .0757
	y_top = .597
	y_bottom = .36
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 12, 1):
		panel1.add_patch(plt.Rectangle((x,y_top), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		panel1.add_patch(plt.Rectangle((x,y_bottom), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .0393

	panel3.add_patch(plt.Rectangle((x,y_top), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	panel3.add_patch(plt.Rectangle((x,y_bottom), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	x +=.0058
	panel3.add_patch(plt.Rectangle((x,y_top), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	panel3.add_patch(plt.Rectangle((x,y_bottom), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	x += .0133
	panel3.add_patch(plt.Rectangle((x,y_top), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
	panel3.add_patch(plt.Rectangle((x,y_bottom), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
	x += .0203
	panel3.add_patch(plt.Rectangle((x,y_top), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
	panel3.add_patch(plt.Rectangle((x,y_bottom), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))

	yText = y_top + .00138
	plt.text(.092, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.1313, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.1706, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.2099, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.2492, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.2885, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3278, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3671, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.4064, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.4457, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.485, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5243, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.5467, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5565, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.573, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5977, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

	yText_labels_top = yText + .015
	yText_labels_bottom = y_bottom + .002
	yText_labels_bottom_sec = yText_labels_bottom - .015

	plt.text(.09, yText_labels_top, '1bp Deletion', fontsize=7, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.167, yText_labels_top, '1bp Insertion', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.266, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.42, yText_labels_top, '>1bp Insertions at Repeats\n       (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.55, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	plt.text(.079, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.156, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.275, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.43, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5475, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	x = .0767
	for i in range (0, 8, 1):
		if i != 2 and i != 3:
			if i == 1 or i == 7:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		else:
			if i == 3:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x += .03925

	for i in range (0, 4, 1):
		if i == 3:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		else:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .03925
	plt.text(x, yText_labels_bottom, '1', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .006
	plt.text(x, yText_labels_bottom, '1  2', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .0135
	plt.text(x, yText_labels_bottom, '1  2  3', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .02
	plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	
	while y%4 != 0:
		y += 1
	if signature:
		ytick_offest = int(y/4)
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]

	else:
		ytick_offest = int(y/4)
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]

	panel3.set_xlim([0, 83])
	panel3.set_yticklabels(ylabels, fontsize=8, color='b')
	panel3.set_ylim([0, y])
	panel3.set_xticks(labs)
	panel3.set_yticks(ylabs)
	panel3.set_yticklabels(ylabels, fontsize=8)
	panel3.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)
	panel3.set_xlabel('')
	panel3.set_ylabel('')
	panel3.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=False, labelleft=True,\
					   right=False, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=25, colors='gray', width=2)

	[i.set_color("black") for i in panel3.get_yticklabels()]
	panel3.set_ylabel("Clustered", fontname='Times New Roman', fontsize=15, fontweight='bold')

	########################################################
	# Non-clustered
	########################################################
	mutations = dict()
	with open (matrix_path_nonClustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		mutations[sample] = {'1DelC':[0,0,0,0,0,0], '1DelT':[0,0,0,0,0,0], '1InsC':[0,0,0,0,0,0], '1InsT':[0,0,0,0,0,0], 
							 '2DelR':[0,0,0,0,0,0], '3DelR':[0,0,0,0,0,0], '4DelR':[0,0,0,0,0,0], '5DelR':[0,0,0,0,0,0],
							 '2InsR':[0,0,0,0,0,0], '3InsR':[0,0,0,0,0,0], '4InsR':[0,0,0,0,0,0], '5InsR':[0,0,0,0,0,0], 
							 '2DelM':[0], '3DelM':[0,0], '4DelM':[0,0,0], '5DelM':[0,0,0,0,0]}

		for lines in f:
			line = lines.strip().split()
			categories = line[0].split(":")
			mut_type = categories[0] + categories[1] + categories[2]
			repeat_size = int(categories[3])
			if categories[2] == 'M':
				repeat_size -= 1

			if mut_type in mutations[sample].keys():
				if percentage:
					mutCount = float(sum([int(x) for x in line[1:]]))
				else:
					mutCount = sum([int(x) for x in line[1:]])
				mutations[sample][mut_type][repeat_size] = mutCount
			else:
				continue

	total_count = sum(sum(nuc) for nuc in mutations[sample].values())
	xlabels = []
	x = 0.5
	ymax = 0
	i = 0
	for key in mutations[sample]:
		l = 1
		for seq in mutations[sample][key]:
			xlabels.append(l)
			if signature:
				if percentage:
					panel5.bar(x, seq*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq*100 > ymax:
						ymax = seq*100

				else:
					panel5.bar(x, seq/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if seq/total_count*100 > ymax:
						ymax = seq/total_count*100
			else:
				panel5.bar(x, seq,width=0.4,color=colors[i],align='center', zorder=1000)
				if seq > ymax:
						ymax = seq
			x += 1
			l += 1
		i += 1

	x = .0757
	y_top = .3
	y_bottom = .0625
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 12, 1):
		panel5.add_patch(plt.Rectangle((x,y_top), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		panel5.add_patch(plt.Rectangle((x,y_bottom), .037, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .0393

	panel5.add_patch(plt.Rectangle((x,y_top), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	panel5.add_patch(plt.Rectangle((x,y_bottom), .0035, .01, facecolor=colors[12], clip_on=False, transform=plt.gcf().transFigure))
	x +=.0058
	panel5.add_patch(plt.Rectangle((x,y_top), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	panel5.add_patch(plt.Rectangle((x,y_bottom), .011, .01, facecolor=colors[13], clip_on=False, transform=plt.gcf().transFigure))
	x += .0133
	panel5.add_patch(plt.Rectangle((x,y_top), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))
	panel5.add_patch(plt.Rectangle((x,y_bottom), .018, .01, facecolor=colors[14], clip_on=False, transform=plt.gcf().transFigure))		
	x += .0203
	panel5.add_patch(plt.Rectangle((x,y_top), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))
	panel5.add_patch(plt.Rectangle((x,y_bottom), .03, .01, facecolor=colors[15], clip_on=False, transform=plt.gcf().transFigure))

	yText = y_top + .00138
	plt.text(.092, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.1313, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.1706, yText, 'C', fontsize=7, fontname='Times New Roman', fontweight='bold', transform=plt.gcf().transFigure)
	plt.text(.2099, yText, 'T', fontsize=7, fontname='Times New Roman', fontweight='bold', color='white', transform=plt.gcf().transFigure)
	plt.text(.2492, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.2885, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3278, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.3671, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.4064, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.4457, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.485, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5243, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	plt.text(.5467, yText, '2', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5565, yText, '3', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.573, yText, '4', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5977, yText, '5+', fontsize=7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)

	yText_labels_top = yText + .015
	yText_labels_bottom = y_bottom + .002
	yText_labels_bottom_sec = yText_labels_bottom - .015

	plt.text(.09, yText_labels_top, '1bp Deletion', fontsize=7, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.167, yText_labels_top, '1bp Insertion', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.266, yText_labels_top, '>1bp Deletion at Repeats\n      (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.42, yText_labels_top, '>1bp Insertions at Repeats\n       (Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.55, yText_labels_top, ' Mircohomology\n(Deletion Length)', fontsize=7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	plt.text(.079, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontname='Times New Roman', weight='bold', color='black', transform=plt.gcf().transFigure)
	plt.text(.156, yText_labels_bottom_sec, 'Homopolymer Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.275, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.43, yText_labels_bottom_sec, 'Number of Repeat Units', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	plt.text(.5475, yText_labels_bottom_sec, 'Mircohomology Length', fontsize=6.5, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

	x = .0767
	for i in range (0, 8, 1):
		if i != 2 and i != 3:
			if i == 1 or i == 7:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '1  2  3  4  5  6+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		else:
			if i == 3:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
			else:
				plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)

		x += .03925

	for i in range (0, 4, 1):
		if i == 3:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
		else:
			plt.text(x, yText_labels_bottom, '0  1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
		x += .03925
	plt.text(x, yText_labels_bottom, '1', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .006
	plt.text(x, yText_labels_bottom, '1  2', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .0135
	plt.text(x, yText_labels_bottom, '1  2  3', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='black', transform=plt.gcf().transFigure)
	x += .02
	plt.text(x, yText_labels_bottom, '1  2  3  4  5+', fontsize=5.7, fontweight='bold', fontname='Times New Roman', color='white', transform=plt.gcf().transFigure)
	
	while y%4 != 0:
		y += 1
	if signature:
		ytick_offest = int(y/4)
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]

	else:
		ytick_offest = int(y/4)
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]


	panel5.set_xlim([0, 83])
	panel5.set_yticklabels(ylabels, fontsize=8, color='b')
	panel5.set_ylim([0, y])
	panel5.set_xticks(labs)
	panel5.set_yticks(ylabs)
	panel5.set_yticklabels(ylabels, fontsize=8)
	panel5.grid(which='major', axis='y', color=[0.6,0.6,0.6], zorder=1)

	fig.suptitle("Summary", fontsize=30, fontname="Times New Roman")
	plt.text(.008, .55, 'Mutation Counts', rotation='vertical', fontsize=20, fontweight='bold', fontname="Times New Roman", transform=plt.gcf().transFigure)
	panel5.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=False, labelleft=True,\
					   right=False, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=25, colors='gray', width=2)
	panel5.set_ylabel("Non-Clustered", fontname='Times New Roman', fontsize=15, fontweight='bold')
	[i.set_color("black") for i in panel5.get_yticklabels()]


def plot96_same_final (matrix_path, matrix_path_clustered, matrix_path_nonClustered, sample, percentage, signature, panel1, panel3, panel5, fig):
	sample = 'All'
	if 'roman' in matplotlib.font_manager.weight_dict:
		del matplotlib.font_manager.weight_dict['roman']
		matplotlib.font_manager._rebuild()

	mutations = dict()
	total_count = []
	with open (matrix_path) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
							 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0]
			mut_type = line[0][2:5]

			if percentage:
				mutCount = float(sum([int(x) for x in line[1:]]))
			else:
				mutCount = sum([int(x) for x in line[1:]])

			mutations[sample][mut_type][nuc] = mutCount


	total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
	xlabels = []

	x = 0.4
	ymax = 0
	colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
	i = 0
	for key in mutations[sample]:
		for seq in mutations[sample][key]:
			xlabels.append(seq[0]+seq[2]+seq[6])
			if signature:
				if percentage:
					panel1.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]*100 > ymax:
						ymax = mutations[sample][key][seq]*100
				else:	
					panel1.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]/total_count*100 > ymax:
						ymax = mutations[sample][key][seq]/total_count*100
			else:
				panel1.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
				if mutations[sample][key][seq] > ymax:
						ymax = mutations[sample][key][seq]
			x += 1
		i += 1

	x = .077
	y3 = .895
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 6, 1):
		panel1.add_patch(plt.Rectangle((x,y3), .088, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .09

	yText = y3 + .015
	plt.text(.11, yText, 'C>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.2, yText, 'C>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.288, yText, 'C>T', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.376, yText, 'T>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.466, yText, 'T>C', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	plt.text(.554, yText, 'T>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

	while y%4 != 0:
		y += 1
	ytick_offest = int(y/4)
	if signature:
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
	else:
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]		

	labs = np.arange(0.375,96.375,1)

	panel1.set_xlim([0, 96])
	panel1.set_ylim([0, y])
	panel1.set_xticks(labs)
	panel1.set_yticks(ylabs)
	count = 0
	m = 0
	for i in range (0, 96, 1):
		plt.text(i/177 + .075, .052, xlabels[i][0], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .06, xlabels[i][1], fontsize=6, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .068, xlabels[i][2], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		count += 1
		if count == 16:
			count = 0
			m += 1	


	panel1.set_yticklabels(ylabels, fontsize=8, color='b')
	panel1.set_xlabel('')
	panel1.set_ylabel('')
	panel1.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
	panel1.set_ylabel("All Mutations", fontname='Times New Roman', fontsize=15, fontweight='bold')
	panel1.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=True, labelleft=True,\
					   right=True, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=10, colors='lightgray', width=1)
	[i.set_color("black") for i in panel1.get_yticklabels()]



	#############################################################################################
	# Plot 96 bar plot (clustered)
	#############################################################################################
	mutations = dict()
	total_count = []
	with open (matrix_path_clustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
							 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0]
			mut_type = line[0][2:5]

			if percentage:
				mutCount = float(sum([int(x) for x in line[1:]]))
			else:
				mutCount = sum([int(x) for x in line[1:]])
			mutations[sample][mut_type][nuc] = mutCount

	total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
	xlabels = []

	x = 0.4
	ymax = 0
	colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
	i = 0
	for key in mutations[sample]:
		for seq in mutations[sample][key]:
			xlabels.append(seq[0]+seq[2]+seq[6])
			if signature:
				if percentage:
					panel3.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]*100 > ymax:
						ymax = mutations[sample][key][seq]*100
				else:	
					panel3.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]/total_count*100 > ymax:
						ymax = mutations[sample][key][seq]/total_count*100
			else:
				panel3.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
				if mutations[sample][key][seq] > ymax:
						ymax = mutations[sample][key][seq]
			x += 1
		i += 1

	x = .077
	y3 = .598
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 6, 1):
		panel3.add_patch(plt.Rectangle((x,y3), .088, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .09

	yText = y3 + .015
	panel3.text(.11, yText, 'C>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.2, yText, 'C>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.288, yText, 'C>T', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.376, yText, 'T>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.466, yText, 'T>C', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel3.text(.554, yText, 'T>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

	panel3.set_ylabel("Clustered", fontname='Times New Roman', fontsize=15, fontweight='bold')


	while y%4 != 0:
		y += 1
	ytick_offest = int(y/4)

	if signature:
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
	else:
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]		

	labs = np.arange(0.375,96.375,1)

	panel3.set_xlim([0, 96])
	panel3.set_yticklabels(ylabels, fontsize=8, color='b')
	panel3.set_ylim([0, y])
	panel3.set_xticks(labs)
	panel3.set_yticks(ylabs)
	count = 0
	m = 0
	for i in range (0, 96, 1):
		plt.text(i/177 + .075, .3495, xlabels[i][0], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .3575, xlabels[i][1], fontsize=6, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .3655, xlabels[i][2], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		count += 1
		if count == 16:
			count = 0
			m += 1

	panel3.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=True, labelleft=True,\
					   right=True, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=10, colors='lightgray', width=1)

	panel3.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)
	[i.set_color("black") for i in panel3.get_yticklabels()]


	#############################################################################################
	# Plot 96 bar plot (non-clustered)
	#############################################################################################
	mutations = dict()
	total_count = []
	with open (matrix_path_nonClustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]
		mutations[sample] = {'C>A':OrderedDict(), 'C>G':OrderedDict(), 'C>T':OrderedDict(),
							 'T>A':OrderedDict(), 'T>C':OrderedDict(), 'T>G':OrderedDict()}

		for lines in f:
			line = lines.strip().split()
			nuc = line[0]
			mut_type = line[0][2:5]

			if percentage:
				mutCount = float(sum([int(x) for x in line[1:]]))
			else:
				mutCount = sum([int(x) for x in line[1:]])
			mutations[sample][mut_type][nuc] = mutCount

	total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
	xlabels = []

	x = 0.4
	ymax = 0
	colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],[228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
	i = 0
	for key in mutations[sample]:
		for seq in mutations[sample][key]:
			xlabels.append(seq[0]+seq[2]+seq[6])
			if signature:
				if percentage:
					panel5.bar(x, mutations[sample][key][seq]*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]*100 > ymax:
						ymax = mutations[sample][key][seq]*100
				else:	
					panel5.bar(x, mutations[sample][key][seq]/total_count*100,width=0.4,color=colors[i],align='center', zorder=1000)
					if mutations[sample][key][seq]/total_count*100 > ymax:
						ymax = mutations[sample][key][seq]/total_count*100
			else:
				panel5.bar(x, mutations[sample][key][seq],width=0.4,color=colors[i],align='center', zorder=1000)
				if mutations[sample][key][seq] > ymax:
						ymax = mutations[sample][key][seq]
			x += 1
		i += 1

	x = .077
	y3 = .301
	y = int(ymax*1.25)
	y2 = y+2
	for i in range(0, 6, 1):
		panel5.add_patch(plt.Rectangle((x,y3), .088, .01, facecolor=colors[i], clip_on=False, transform=plt.gcf().transFigure)) 
		x += .09

	yText = y3 + .015
	panel5.text(.11, yText, 'C>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.2, yText, 'C>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.288, yText, 'C>T', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.376, yText, 'T>A', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.466, yText, 'T>C', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)
	panel5.text(.554, yText, 'T>G', fontsize=12, fontweight='bold', fontname='Arial', transform=plt.gcf().transFigure)

	panel5.set_ylabel("Non-Clustered", fontname='Times New Roman', fontsize=15, fontweight='bold')

	while y%4 != 0:
		y += 1
	ytick_offest = int(y/4)

	if signature:
		ylabs = [0, round(ytick_offest, 1), round(ytick_offest*2, 1), round(ytick_offest*3, 1), round(ytick_offest*4, 1)]
		ylabels= [str(0), str(round(ytick_offest, 1)) + "%", str(round(ytick_offest*2, 1)) + "%", 
				  str(round(ytick_offest*3, 1)) + "%", str(round(ytick_offest*4, 1)) + "%"]
	else:
		ylabs = [0, ytick_offest, ytick_offest*2, ytick_offest*3, ytick_offest*4]
		ylabels= [0, ytick_offest, ytick_offest*2, 
			  	  ytick_offest*3, ytick_offest*4]		

	labs = np.arange(0.375,96.375,1)

	panel5.set_xlim([0, 96])
	panel5.set_yticklabels(ylabels, fontsize=8, color='b')
	panel5.set_ylim([0, y])
	panel5.set_xticks(labs)
	panel5.set_yticks(ylabs)
	panel5.grid(which='major', axis='y', color=[0.93,0.93,0.93], zorder=1)

	count = 0
	m = 0
	for i in range (0, 96, 1):
		plt.text(i/177 + .075, .647, xlabels[i][0], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .655, xlabels[i][1], fontsize=6, color=colors[m], rotation='vertical', verticalalignment='center', fontname='Courier New', fontweight='bold',transform=plt.gcf().transFigure)
		plt.text(i/177 + .075, .663, xlabels[i][2], fontsize=6, color='black', rotation='vertical', verticalalignment='center', fontname='Courier New', transform=plt.gcf().transFigure)
		count += 1
		if count == 16:
			count = 0
			m += 1

	panel5.tick_params(axis='both',which='both',\
					   bottom=False, labelbottom=False,\
					   left=True, labelleft=True,\
					   right=True, labelright=False,\
					   top=False, labeltop=False,\
					   direction='in', length=10, colors='lightgray', width=1)

	[i.set_color("black") for i in panel5.get_yticklabels()]

	fig.suptitle(sample, fontsize=30, fontname="Times New Roman")
	plt.text(.008, .55, 'Mutation Counts', rotation='vertical', fontsize=20, fontweight='bold', fontname="Times New Roman", transform=plt.gcf().transFigure)


def first_run (distances, distances_orig, original_vcf, vcf_path_clust, vcf_path_nonClust, sample, original, sim_count, matrix_path_clustered, matrix_path_nonClustered, cutoff, vcf_path_all_clust, vcf_path_all_nonClust, project, distance_path, genome, clustering_vaf):

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
	upper_conf = []
	lower_conf = []
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

		upper_conf.append(bin_std_dev[CI-1]-int(round(median(bin_std_dev),0)))
		lower_conf.append(int(round(median(bin_std_dev),0)) - bin_std_dev[-CI])
		upper_CI.append(bin_std_dev[CI-1])
		lower_CI.append(bin_std_dev[-CI])


	y, binEdges = np.histogram(avg_bin_counts, bins = bin1)
	y2, binEdges2 = np.histogram(distances_orig, bins = bin1)

	bincenters2 = binEdges2[1:-1]
	y2 = y2[1:]
	avg_bin_counts = avg_bin_counts[1:]
	mean_counts = mean_counts[1:]
	std_dev = std_dev[1:]

	upper_conf=upper_conf[1:]
	lower_conf=lower_conf[1:]
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

	cluster = 0
	if True:
		interval_line = 10
		distance_cut = bincenters2[interval_line]
		with open (original_vcf) as f, open(vcf_path_clust + project + "_clustered.txt", 'a') as clust, open(vcf_path_nonClust + project + "_nonClustered.txt", 'a') as non, open(distance_path + "distances_" + sample + ".txt", 'a') as dist_out:
			initial_line = True
			initial_mutation = True
			first_sample_occurence = True
			previous_line = None

			prev_dist = 0

			for lines in f:
				line = lines.strip().split()
				sample2 = line[0]
				chrom = line[1]
				start = int(line[2])
				ref = line[3]
				mut = line[4]
				if clustering_vaf:
					vaf = line[5]

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
				if sample2 != sample:
					if first_sample_occurence:
						continue
					else:
						break
				else:
					first_sample_occurence = False
					if initial_line:
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
						if chrom != init_chrom:
							if prev_dist <= distance_cut:
								cluster += 1
								print('\t'.join([x for x in previous_line]), file=clust)
							else:
								print('\t'.join([x for x in previous_line]), file=non)
							print(str(prev_dist), file=dist_out)

							init_chrom = chrom
							init_start = start
							init_ref = ref
							init_mut = mut
							previous_line = line
							record_prev = True
							centro_start = centromeres[genome][init_chrom][0]
							centro_end = centromeres[genome][init_chrom][1]
							initial_mutation = True
							continue

						else:
							if init_start < centro_start and start > centro_start:
								if prev_dist <= distance_cut:
									cluster += 1
									print('\t'.join([x for x in previous_line]), file=clust)
								else:
									print('\t'.join([x for x in previous_line]), file=non)
								print(str(prev_dist), file=dist_out)
								initial_line = True
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
								if dist <= distance_cut:
									cluster += 1
									print('\t'.join([x for x in previous_line]), file=clust)
								else:
									print('\t'.join([x for x in previous_line]), file=non)

								print(str(dist), file=dist_out)
								initial_mutation = False
							else:
								dist_min = min(dist, prev_dist)
								if dist_min <= distance_cut:
									cluster += 1
									print('\t'.join([x for x in previous_line]), file=clust)

								else:
									print('\t'.join([x for x in previous_line]), file=non)

								print(str(dist_min), file=dist_out)

							previous_line = line
							prev_dist = dist
				init_start = start
			if prev_dist <= distance_cut:
				cluster += 1
				print('\t'.join([x for x in previous_line]), file=clust)

			else:
				print('\t'.join([x for x in previous_line]), file=non)

			print(str(prev_dist), file=dist_out)



def hotSpotAnalysis (project, genome, contexts, simContext, ref_dir, original=False, signature=False, percentage=False, firstRun=False, clustering_vaf=False):
	height = 8
	width = 13
	scaled_width = (width/1.75 *.95)/width
	scaled_height = (height/4.5)/height

	if ref_dir[-1] != "/":
		ref_dir += "/"

	simContext = sorted(simContext, reverse=True)
	simContext = "_".join(simContext)
	file_context = contexts
	if file_context == 'INDEL':
		cutoff = 25
	else:
		cutoff = 50

	if contexts == '96':
		matrix_file_suffix = '.SBS96.'
	elif contexts ==  'INDEL':
		matrix_file_suffix = '.DBS94.'


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

	else:
		vcf_path_clust = ref_dir + "/references/vcf_files/" + project + "_clustered/INDEL/"
		vcf_path_nonClust = ref_dir + "/references/vcf_files/" + project + "_nonClustered/INDEL/"
		vcf_path_all_clust = ref_dir + "/references/vcf_files/" + project + "_all_clustered/INDEL/"
		vcf_path_all_nonClust = ref_dir + "/references/vcf_files/" + project + "_all_nonClustered/INDEL/"	
	

	matrix_path = ref_dir + "output/SBS/" + project + matrix_file_suffix + "all"
	matrix_path_clustered = ref_dir + "output/vcf_files/" + project + "_clustered/SNV/output/SBS/" + project + "_clustered" + matrix_file_suffix + "all" 
	matrix_path_nonClustered = ref_dir + "output/vcf_files/" + project + "_nonClustered/SNV/output/SBS/" + project + "_nonClustered" + matrix_file_suffix + "all" 
	matrix_path_all_clustered = ref_dir + "output/vcf_files/" + project + "_all_clustered/SNV/output/SBS" + project + "_all_clustered" + matrix_file_suffix + "all"
	matrix_path_all_nonClustered = ref_dir + "/references/matrix/" + project + "_all_nonClustered/" + project + "_all_nonClustered" + matrix_file_suffix + "all"
	output_path = directory_out + project + '_intradistance_plots_' + contexts + '.pdf'


	if os.path.exists(directory_out) == False:
		os.mkdir(directory_out)

	if firstRun:
		if os.path.exists(vcf_path_clust):
			shutil.rmtree(vcf_path_clust)

	if not os.path.exists(vcf_path_clust):
		os.makedirs(vcf_path_clust)

	if not os.path.exists(vcf_path_nonClust):
		os.makedirs(vcf_path_nonClust)


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
			overall_distances = []
			distances_orig = []
			distances_avg = []
			distances_orig_adj = []

			if original:
				if first_sort:
					original_vcfs = os.listdir(ref_dir + "output/vcf_files/single/SNV/")
					with open(ref_dir + "output/vcf_files/single/" + project + "_all.txt", "w") as out:
						for f in original_vcfs:
							if f == '.DS_Store':
								continue
							with open(ref_dir + "output/vcf_files/single/SNV/" + f) as fd:
								shutil.copyfileobj(fd, out)

					with open(original_vcf) as f:
					    lines = [line.strip().split() for line in f]
					output = open(original_vcf, 'w')
					for line in sorted(lines, key = lambda x: (x[0], ['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'].index(x[1]), int(x[2]))):
					    print('\t'.join(line), file=output)
					output.close()
					first_sort = False

				with open (directory_orig + sample + "_" + contexts + "_intradistance.txt") as f2:
					for lines in f2:
						line = lines.strip().split()
						if int(line[0]) >= 1:
							distances_orig.append(int(line[0]))
							distances_orig_adj.append(int(line[1]))

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
				overall_distances.append(distances)
			first_run(overall_distances, distances_orig, original_vcf, vcf_path_clust, vcf_path_nonClust, sample, original, sim_count, matrix_path_clustered, matrix_path_nonClustered, cutoff, vcf_path_all_clust, vcf_path_all_nonClust, project, distance_path, genome, clustering_vaf)
		
		with open(vcf_path_clust + project + "_clustered.txt") as clust:	
			clust_lines = clust.readlines()

		with open(vcf_path_nonClust + project + "_nonClustered.txt") as nonclust:
			nonclust_lines = nonclust.readlines()

		f = open(vcf_path_clust + project + "_clustered.txt", "w")
		for lines in clust_lines:
			line = lines.strip().split("\t")
			sample = line[0]
			chrom = line[1]
			pos = line[2]
			ref  = line[3]
			mut = line[4]
			if clustering_vaf:
				vaf = line[5]
				print("\t".join([project,sample,".",genome,"SNP",chrom,pos,pos,ref,mut,"SOMATIC",vaf]),file=f)
			else:
				print("\t".join([project,sample,".",genome,"SNP",chrom,pos,pos,ref,mut,"SOMATIC"]),file=f)

		f.close()

		f = open(vcf_path_nonClust + project + "_nonClustered.txt", "w")
		for lines in nonclust_lines:
			line = lines.strip().split("\t")
			sample = line[0]
			chrom = line[1]
			pos = line[2]
			ref  = line[3]
			mut = line[4]
			print("\t".join([project,sample,".",genome,"SNP",chrom,pos,pos,ref,mut,"SOMATIC"]),file=f)
		f.close()

		print("Completed!", flush=True)

		if file_context == 'INDEL':
			os.system("python3 sigProfilerMatrixGenerator_bi_2_dinuc_context.py -p " + project + "_clustered -g " + genome + " -i")
			os.system("python3 sigProfilerMatrixGenerator_bi_2_dinuc_context.py -p " + project + "_nonClustered -g " + genome + " -i")
		else:
			print("\nAnalyzing clustered mutations...", flush=True)
			matGen.SigProfilerMatrixGeneratorFunc(project + "_clustered", genome, vcf_path_clust,plot=False)
			print("\nAnalyzing non-clustered mutations...", flush=True)
			matGen.SigProfilerMatrixGeneratorFunc(project + "_nonClustered", genome, vcf_path_nonClust,plot=False)
			print(flush=True)

	with open(matrix_path_clustered) as f:
		first_line = f.readline()
		samples = first_line.strip().split()
		samples = samples[1:]

	pp = PdfPages(directory_out + project + '_intradistance_plots_' + simContext + '.pdf')

	l = 0
	all_distances_orig = []
	all_distances_sim = []
	histo = True

	print("Plotting SigProfilerHotSpot Results...", end='', flush=True)

	for folder in folders:
		if folder == '.DS_Store_intradistance.txt' or folder == '.DS_Store':
			continue
		if folder not in samples:
			histo = False
		sample = folder
		files = os.listdir(directory + sample + "/")
		overall_distances = []
		distances_orig = []
		distances_avg = []
		# distances_orig_adj = []

		fig = plt.figure(figsize = (width, height))
		panel1=plt.axes([0.075, 0.225 + scaled_height*2, scaled_width, scaled_height])
		panel2=plt.axes([0.125 + scaled_width, 0.225 + scaled_height*2, 0.3, scaled_height])
		panel3=plt.axes([0.075, 0.15 + scaled_height, scaled_width, scaled_height])
		panel4=plt.axes([0.125 + scaled_width, 0.15 + scaled_height, 0.3, scaled_height])
		panel5=plt.axes([0.075, 0.075, scaled_width, scaled_height])
		panel6=plt.axes([0.125 + scaled_width, 0.075, 0.3, scaled_height])


		if original:
			if first_sort:
				with open(original_vcf) as f:
				    lines = [line.strip().split() for line in f]
				output = open(original_vcf, 'w')
				for line in sorted(lines, key = lambda x: (x[1], ['X','Y','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'].index(x[5]), int(x[6]))):
				    print('\t'.join(line), file=output)
				output.close()
				first_sort = False
			with open (distance_path + "distances_" + sample + ".txt") as f2:
				for lines in f2:
					line = lines.strip().split()
					if int(line[0]) >= 1:
						distances_orig.append(int(line[0]))
						# distances_orig_adj.append(int(line[1]))
						all_distances_orig.append(int(line[0]))


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
			overall_distances.append(distances)
			all_distances_sim.append(distances)


		if histo:
			clustered = plot_hist(overall_distances, distances_orig, distances_orig_adj, original_vcf, vcf_path_clust, vcf_path_nonClust, sample, pp, original, sim_count, matrix_path_clustered, matrix_path_nonClustered, panel2, panel3, panel4, panel5, panel6, cutoff)
			if clustered:
				if file_context == '96':
					plot96_same (matrix_path, matrix_path_clustered, matrix_path_nonClustered, sample, percentage, signature, panel1, panel3, panel5, fig)
				else:
					plotINDEL_same (matrix_path, matrix_path_clustered, matrix_path_nonClustered, sample, percentage, signature, panel1, panel3, panel5, fig)
				l +=1
				pp.savefig()
		
			plt.close()
		
		histo = True

	plt.close()
	pp.close()
	print("Completed!", flush=True)

if __name__ == '__main__':
	main()
