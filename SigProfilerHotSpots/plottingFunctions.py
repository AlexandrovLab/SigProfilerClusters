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