import os
import numpy as np
import os
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import shutil
from numpy import median
import pickle


def pullVaf (project, project_path, sanger=True, TCGA=False):
	
	vcf_path = project_path + "input/"
	clusteredMutsPath = project_path + "output/vcf_files/" + project + "_clustered/"
	clusteredMutsFile = project_path + "output/vcf_files/" + project + "_clustered" + "/SNV/" + project + "_clustered.txt"
	vcf_files = [x for x in os.listdir(vcf_path) if x != ".DS_Store"]


	if sanger:
		vafs = {}
		for vcfFile in vcf_files:
			sample = vcfFile.split(".")[0]
			vafs[sample] = {}
			with open(vcf_path + vcfFile) as f:
				for lines in f:
					if lines[0] == "#":
						continue
					lines = lines.strip().split()
					chrom = lines[0]
					pos = lines[1]
					ref = lines[3]
					alt = lines[4]
					try:
						vaf = float(lines[-1].split(":")[-1])
					except:
						print("There does not seem to be VAF scores in this Sanger-produced file.")
						break
					keyLine = ":".join([chrom, pos, ref, alt])
					vafs[sample][keyLine] = vaf

		with open(clusteredMutsFile) as f, open(clusteredMutsPath + project + "_clustered_vaf.txt", "w") as out:
			next(f)
			print("HEADER", file=out)
			for lines in f:
				lines = lines.strip().split()
				sample = lines[1]
				newKey = ":".join([lines[5], lines[6], lines[8], lines[9]])

				vaf = vafs[sample][newKey]

				lines.append(str(vaf))
				print("\t".join([x for x in lines]), file=out)

	elif TCGA:
		vafs = {}
		for vcfFile in vcf_files:
			sample = vcfFile.split(".")[0]
			vafs[sample] = {}
			with open(vcf_path + vcfFile) as f:
				for lines in f:
					if lines[0] == "#":
						continue
					lines = lines.strip().split()
					chrom = lines[0]
					pos = lines[1]
					ref = lines[3]
					alt = lines[4]
					try:
						vaf = float(lines[7].split("VAF=")[1].split(";")[0])
					except:
						vaf = -1.5
					keyLine = ":".join([chrom, pos, ref, alt])
					vafs[sample][keyLine] = vaf

		with open(clusteredMutsFile) as f, open(clusteredMutsPath + project + "_clustered_vaf.txt", "w") as out:
			next(f)
			print("HEADER", file=out)
			for lines in f:
				lines = lines.strip().split()
				sample = lines[1]
				newKey = ":".join([lines[5], lines[6], lines[8], lines[9]])

				vaf = vafs[sample][newKey]

				lines.append(str(vaf))
				print("\t".join([x for x in lines]), file=out)







def findClustersOfClusters (project, chrom_based, project_parent_path):
	project_path = project_parent_path + "output/vcf_files/ESCC_clustered/"
	file  = project_path + project + '_clustered_vaf.txt'
	out_file = project_path + project + '_clusters_of_clusters.txt'
	out_file2 = project_path + project + '_clusters_of_clusters_imd.txt'

	path_suffix = ''
	if chrom_based:
		path_suffix = "_chrom"

	out_file3 = project_path + 'subclasses' + path_suffix + '/class1/' + project + '_clustered_class1.txt'
	out_file4 = project_path + 'subclasses' + path_suffix + '/class2/' + project + '_clustered_class2.txt'
	out_file5 = project_path + 'subclasses' + path_suffix + '/class3/' + project + '_clustered_class3.txt'
	out_file6 = project_path + 'subclasses' + path_suffix + '/class1a/' + project + '_clustered_class1a.txt'
	out_file7 = project_path + 'subclasses' + path_suffix + '/class1b/' + project + '_clustered_class1b.txt'
	out_file8 = project_path + 'subclasses' + path_suffix + '/class1c/' + project + '_clustered_class1c.txt'


	if os.path.exists(project_path + 'subclasses' + path_suffix + '/class1/'):
		shutil.rmtree(project_path + 'subclasses' + path_suffix + '/class1/')
	os.makedirs(project_path + 'subclasses' + path_suffix + '/class1/')
	if os.path.exists(project_path + 'subclasses' + path_suffix + '/class2/'):
		shutil.rmtree(project_path + 'subclasses' + path_suffix + '/class2/')
	os.makedirs(project_path + 'subclasses' + path_suffix + '/class2/')
	if os.path.exists(project_path + 'subclasses' + path_suffix + '/class3/'):
		shutil.rmtree(project_path + 'subclasses' + path_suffix + '/class3/')
	os.makedirs(project_path + 'subclasses' + path_suffix + '/class3/')
	if os.path.exists(project_path + 'subclasses' + path_suffix + '/class1a/'):
		shutil.rmtree(project_path + 'subclasses' + path_suffix + '/class1a/')
	os.makedirs(project_path + 'subclasses' + path_suffix + '/class1a/')
	if os.path.exists(project_path + 'subclasses' + path_suffix + '/class1b/'):
		shutil.rmtree(project_path + 'subclasses' + path_suffix + '/class1b/')
	os.makedirs(project_path + 'subclasses' + path_suffix + '/class1b/')
	if os.path.exists(project_path + 'subclasses' + path_suffix + '/class1c/'):
		shutil.rmtree(project_path + 'subclasses' + path_suffix + '/class1c/')
	os.makedirs(project_path + 'subclasses' + path_suffix + '/class1c/')

	cutoff = 10001
	imd = 21
	vaf_cut = 0.1

	if chrom_based:
		with open(project_parent_path + "output/simulations/data/imds_chrom.pickle", "rb") as handle:
			imdsData = pickle.load(handle) 
	else:
		with open(project_parent_path + "output/simulations/data/imds.pickle", "rb") as handle:
			imdsData = pickle.load(handle) 

	first_pass = True
	total_muts = 0
	if first_pass:
		mnv_length = 0
		len_mnvs = {}
		total_mnvs = {}
		distances = []
		count = 1
		out = open(out_file, 'w')
		with open(file) as f:
			next(f)
			lines = [line.strip().split() for line in f]

		for i in range(1, len(lines), 1):
			prev_chrom = lines[i-1][5]
			prev_pos = int(lines[i-1][6])
			prev_samp = lines[i-1][1]
			chrom = lines[i][5]
			pos = int(lines[i][6])
			samp = lines[i][1]

			if prev_samp == samp:
				if prev_chrom == chrom:
					if pos - prev_pos < cutoff:
						distances.append(pos-prev_pos)
						mnv_length += 1
						lines[i-1] = [str(count)] + lines[i-1]
						print("\t".join([x for x in lines[i-1]]), file=out)
						total_muts += 1
					else:
						mnv_length += 1
						if samp not in len_mnvs:
							len_mnvs[samp] = {}
						if str(mnv_length) not in len_mnvs[samp]:
							len_mnvs[samp][str(mnv_length)] = 1
						else:
							len_mnvs[samp][str(mnv_length)] += 1
						if str(mnv_length) not in total_mnvs:
							total_mnvs[str(mnv_length)] = 1
						else:
							total_mnvs[str(mnv_length)] += 1
						mnv_length = 0
						lines[i-1] = [str(count)] + lines[i-1]
						print("\t".join([x for x in lines[i-1]]), file=out)
						total_muts += 1
						count += 1
						print("\n\n", file=out)
				else:
					mnv_length += 1
					if samp not in len_mnvs:
						len_mnvs[samp] = {}
					if str(mnv_length) not in len_mnvs[samp]:
						len_mnvs[samp][str(mnv_length)] = 1
					else:
						len_mnvs[samp][str(mnv_length)] += 1
					if str(mnv_length) not in total_mnvs:
						total_mnvs[str(mnv_length)] = 1
					else:
						total_mnvs[str(mnv_length)] += 1
					mnv_length = 0

					lines[i-1] = [str(count)] + lines[i-1]
					print("\t".join([x for x in lines[i-1]]), file=out)
					total_muts += 1
					count += 1
					print("\n\n", file=out)
			else:
				mnv_length += 1
				if prev_samp not in len_mnvs:
					len_mnvs[prev_samp] = {}
				if str(mnv_length) not in len_mnvs[prev_samp]:
					len_mnvs[prev_samp][str(mnv_length)] = 1
				else:
					len_mnvs[prev_samp][str(mnv_length)] += 1
				if str(mnv_length) not in total_mnvs:
					total_mnvs[str(mnv_length)] = 1
				else:
					total_mnvs[str(mnv_length)] += 1
				mnv_length = 0

				lines[i-1] = [str(count)] + lines[i-1]
				print("\t".join([x for x in lines[i-1]]), file=out)
				total_muts += 1
				count += 1
				count = 1
				print("\n\n################ New Sample #################", file=out)
				print("\n\n", file=out)

		lines[i] = [str(count)] + lines[i]
		print("\t".join([x for x in lines[i]]), file=out)
		total_muts += 1
		out.close()


	if True:
		len_mnvs = {'I':{}, 'II':{},'III':{},'Ia':{},'Ib':{},'Ic':{}}
		distances = []
		distances_mnv = {}
		lines = []
		count = 1
		with open(out_file) as f, open(out_file2, 'w') as out2, open(out_file3, 'w') as out3, open(out_file4, 'w') as out4, open(out_file5, 'w') as out5, open(out_file6, 'w') as out6, open(out_file7, 'w') as out7, open(out_file8, 'w') as out8:
			for line in f:
				line = line.strip().split()[1:]
				if line != []:
					lines.append(line)
				else:
					write_out = False
					category = None
					writeClassI = False
					writeClassII = False
					writeClassIII = False
					writeClassII = False
					writeClassIb = False
					writeClassIc = False
					writeClassIa = False

					if len(lines) > 0:
						if lines[-1][0] == "New":
							lines = lines[1:]
							count = 1
							write_out = False
						if len(lines) == 1 or len(lines) == 0:
							lines = []
							continue
						else:
							# imd = imdsData[lines[0][2]]
							distancesLine = len([int(y[7])-int(x[7]) for x,y in zip(lines, lines[1:]) if int(y[7])-int(x[7]) > 1])
							vafs = len([abs(float(y[11]) - float(x[11])) for x,y in zip(lines, lines[1:]) if round(abs(float(y[11]) - float(x[11])),4) > vaf_cut])
							unknownVafs = len([abs(float(y[11]) - float(x[11])) for x,y in zip(lines, lines[1:]) if abs(float(y[11]) - float(x[11])) < 0])
							distancesActual = [int(y[7])-int(x[7]) for x,y in zip(lines, lines[1:])]
							vafsActual = [abs(float(y[11]) - float(x[11])) for x,y in zip(lines, lines[1:])]
							# imdsDistances = len([int(y[7])-int(x[7]) for x,y in zip(lines, lines[1:]) if int(y[7])-int(x[7]) > imd])

							if unknownVafs == 0 and vafs == 0:
								# Class I or II
								if distancesLine > 1:
									writeClassII = True
									if str(len(lines)) not in len_mnvs['II']:
										len_mnvs['II'][str(len(lines))] = [1,distancesActual, vafsActual]
									else:
										len_mnvs['II'][str(len(lines))][0] += 1
										len_mnvs['II'][str(len(lines))][1] += distancesActual
										len_mnvs['II'][str(len(lines))][2] += vafsActual
								else:
									writeClassI = True
									if str(len(lines)) not in len_mnvs['I']:
										len_mnvs['I'][str(len(lines))] = [1,distancesActual, vafsActual]
									else:
										len_mnvs['I'][str(len(lines))][0] += 1
										len_mnvs['I'][str(len(lines))][1] += distancesActual
										len_mnvs['I'][str(len(lines))][2] += vafsActual
									if distancesLine == 0:
										if len(lines) == 2:
											writeClassIa = True
											if str(len(lines)) not in len_mnvs['Ia']:
												len_mnvs['Ia'][str(len(lines))] = [1,distancesActual, vafsActual]
											else:
												len_mnvs['Ia'][str(len(lines))][0] += 1
												len_mnvs['Ia'][str(len(lines))][1] += distancesActual
												len_mnvs['Ia'][str(len(lines))][2] += vafsActual
										else:
											writeClassIb = True
											if str(len(lines)) not in len_mnvs['Ib']:
												len_mnvs['Ib'][str(len(lines))] = [1,distancesActual, vafsActual]
											else:
												len_mnvs['Ib'][str(len(lines))][0] += 1
												len_mnvs['Ib'][str(len(lines))][1] += distancesActual
												len_mnvs['Ib'][str(len(lines))][2] += vafsActual
									else:
										writeClassIc = True
										if str(len(lines)) not in len_mnvs['Ic']:
											len_mnvs['Ic'][str(len(lines))] = [1,distancesActual, vafsActual]
										else:
											len_mnvs['Ic'][str(len(lines))][0] += 1
											len_mnvs['Ic'][str(len(lines))][1] += distancesActual
											len_mnvs['Ic'][str(len(lines))][2] += vafsActual																			
							else:
								writeClassIII = True

						
						if writeClassII:
							for i in range(0, len(lines), 1):
								lines[i].append("ClassII")
								print("\t".join([x for x in lines[i]]), file=out4)
								lines[i] = [str(count)] + lines[i]
								print("\t".join([x for x in lines[i]]), file=out2)
							count += 1
							print("\n\n", file=out2)
						else:
							if writeClassI:
								# Writes Class I (Single Events)
								try:
									for i in range(0, len(lines), 1):
										lines[i].append("ClassI")
										print("\t".join([x for x in lines[i]]), file=out3)
								except:
									print(lines)

								if writeClassIc:
									# Writes Class Ic (extended MBSs) 
									try:
										for i in range(0, len(lines), 1):
											lines[i][-1] = "ClassIC"
											print("\t".join([x for x in lines[i]]), file=out8)
									except:
										print(lines)	
								elif writeClassIa:
									# Writes Class Ia (DBSs) 
									try:
										for i in range(0, len(lines), 1):
											lines[i][-1] = "ClassIA"
											print("\t".join([x for x in lines[i]]), file=out6)
									except:
										print(lines)	
								elif writeClassIb:
									# print("yes")
									# Writes Class 1b (MBSs)
									try:
										for i in range(0, len(lines), 1):
											lines[i][-1] = "ClassIB"
											# lines[i].append(category)
											print("\t".join([x for x in lines[i]]), file=out7)
										# print("yes", lines)
									except:
										print(lines)
							elif writeClassIII:
								# Writes Class III (all other mutations - leftovers) 
								linesSubClass = lines[:]
								while len(linesSubClass) > 1:
									writeClassI = False
									writeClassII = False
									writeClassIII = False
									writeClassII = False
									writeClassIb = False
									writeClassIc = False
									writeClassIa = False
									saveNewEvent = [linesSubClass[0]]
									lineRef = linesSubClass[0]
									for i in range(1, len(linesSubClass), 1):
										# try:
											if chrom_based:
												if abs(float(linesSubClass[i][11]) - float(lineRef[11])) < vaf_cut and int(linesSubClass[i][7])-int(lineRef[7]) <= imdsData[lineRef[1]][lineRef[5]]:
													saveNewEvent.append(linesSubClass[i])
													lineRef = linesSubClass[i]										
											else:
										# except:
												if abs(float(linesSubClass[i][11]) - float(lineRef[11])) < vaf_cut and int(linesSubClass[i][7])-int(lineRef[7]) <= imdsData[lineRef[1]]:
													saveNewEvent.append(linesSubClass[i])
													lineRef = linesSubClass[i]	
											

									if len(saveNewEvent) > 1:
										distancesLine = len([int(y[7])-int(x[7]) for x,y in zip(saveNewEvent, saveNewEvent[1:]) if int(y[7])-int(x[7]) > 1])
										vafs = len([abs(float(y[11]) - float(x[11])) for x,y in zip(saveNewEvent, saveNewEvent[1:]) if round(abs(float(y[11]) - float(x[11])),4) > vaf_cut])
										unknownVafs = len([abs(float(y[11]) - float(x[11])) for x,y in zip(saveNewEvent, saveNewEvent[1:]) if abs(float(y[11]) - float(x[11])) < 0])
										distancesActual = [int(y[7])-int(x[7]) for x,y in zip(saveNewEvent, saveNewEvent[1:])]
										vafsActual = [abs(float(y[11]) - float(x[11])) for x,y in zip(saveNewEvent, saveNewEvent[1:])]
										# Class I or II
										if distancesLine > 1:
											writeClassII = True
										else:
											writeClassI = True
											if str(len(saveNewEvent)) not in len_mnvs['I']:
												len_mnvs['I'][str(len(saveNewEvent))] = [1,distancesActual, vafsActual]
											else:
												len_mnvs['I'][str(len(saveNewEvent))][0] += 1
												len_mnvs['I'][str(len(saveNewEvent))][1] += distancesActual
												len_mnvs['I'][str(len(saveNewEvent))][2] += vafsActual
											if distancesLine == 0:
												if len(saveNewEvent) == 2:
													writeClassIa = True
													if str(len(saveNewEvent)) not in len_mnvs['Ia']:
														len_mnvs['Ia'][str(len(saveNewEvent))] = [1,distancesActual, vafsActual]
													else:
														len_mnvs['Ia'][str(len(saveNewEvent))][0] += 1
														len_mnvs['Ia'][str(len(saveNewEvent))][1] += distancesActual
														len_mnvs['Ia'][str(len(saveNewEvent))][2] += vafsActual
												else:
													writeClassIb = True
													if str(len(saveNewEvent)) not in len_mnvs['Ib']:
														len_mnvs['Ib'][str(len(saveNewEvent))] = [1,distancesActual, vafsActual]
													else:
														len_mnvs['Ib'][str(len(saveNewEvent))][0] += 1
														len_mnvs['Ib'][str(len(saveNewEvent))][1] += distancesActual
														len_mnvs['Ib'][str(len(saveNewEvent))][2] += vafsActual
											else:
												writeClassIc = True
												if str(len(saveNewEvent)) not in len_mnvs['Ic']:
													len_mnvs['Ic'][str(len(saveNewEvent))] = [1,distancesActual, vafsActual]
												else:
													len_mnvs['Ic'][str(len(saveNewEvent))][0] += 1
													len_mnvs['Ic'][str(len(saveNewEvent))][1] += distancesActual
													len_mnvs['Ic'][str(len(saveNewEvent))][2] += vafsActual
									else:
										writeClassIII = True
										if unknownVafs > 0:
											category = 'unknown'
										else:
											category = "vaf"
										if str(len(saveNewEvent)) not in len_mnvs['III']:
											len_mnvs['III'][str(len(saveNewEvent))] = [1,distancesActual, vafsActual]
										else:
											len_mnvs['III'][str(len(saveNewEvent))][0] += 1
											len_mnvs['III'][str(len(saveNewEvent))][1] += distancesActual
											len_mnvs['III'][str(len(saveNewEvent))][2] += vafsActual
									for line in saveNewEvent:
										linesSubClass.remove(line)

									if writeClassII:
										for i in range(0, len(saveNewEvent), 1):
											saveNewEvent[i].append("ClassII")
											print("\t".join([x for x in saveNewEvent[i]]), file=out4)
											saveNewEvent[i] = [str(count)] + saveNewEvent[i]
											print("\t".join([x for x in saveNewEvent[i]]), file=out2)
										count += 1
										print("\n\n", file=out2)
									else:
										if writeClassI:
											# Writes Class I (Single Events)
											try:
												for i in range(0, len(saveNewEvent), 1):
													saveNewEvent[i].append("ClassI")
													print("\t".join([x for x in saveNewEvent[i]]), file=out3)
											except:
												print(saveNewEvent)

											if writeClassIc:
												# Writes Class Ic (extended MBSs) 
												try:
													for i in range(0, len(saveNewEvent), 1):
														saveNewEvent[i][-1] = "ClassIC"
														print("\t".join([x for x in saveNewEvent[i]]), file=out8)
												except:
													print(saveNewEvent)	
											elif writeClassIa:
												# Writes Class Ia (DBSs) 
												try:
													for i in range(0, len(saveNewEvent), 1):
														saveNewEvent[i][-1] = "ClassIA"
														print("\t".join([x for x in saveNewEvent[i]]), file=out6)
												except:
													print(saveNewEvent)	
											elif writeClassIb:
												# print("yes")
												# Writes Class 1b (MBSs)
												try:
													for i in range(0, len(saveNewEvent), 1):
														saveNewEvent[i][-1] = "ClassIB"
														# lines[i].append(category)
														print("\t".join([x for x in saveNewEvent[i]]), file=out7)
													# print("yes", saveNewEvent)
												except:
													print(saveNewEvent)

										elif writeClassIII:
											# print("yes")
											try:
												for i in range(0, len(saveNewEvent), 1):
													saveNewEvent[i].append("ClassIII")
													saveNewEvent[i].append(category)
													print("\t".join([x for x in saveNewEvent[i]]), file=out5)
											except:
												print(saveNewEvent)	
								if len(linesSubClass) != 0:
									try:
										category = "vaf"
										for i in range(0, len(linesSubClass), 1):
											linesSubClass[i].append("ClassIII")
											linesSubClass[i].append(category)
											print("\t".join([x for x in linesSubClass[i]]), file=out5)
									except:
										print(linesSubClass)	


					lines = []


		with open(project_path + "clusteredStats" + path_suffix + ".pickle", "wb") as f:
			pickle.dump(len_mnvs, f)
		# for classes in len_mnvs:
		# 	# print(classes)
		# 	for leng in len_mnvs[classes]:
		# 		print("\t",leng, len_mnvs[classes][leng][0], median(len_mnvs[classes][leng][1]))


	try:
		print("Generating matrices for Class 1 mutations:")
		matrices = matGen.SigProfilerMatrixGeneratorFunc("class1", "GRCh37", project_path + 'subclasses'+ path_suffix + '/class1/', seqInfo=True)#, plot=True)
		print()
	except:
		pass
	try:
		print("Generating matrices for Class 2 mutations:")
		matrices = matGen.SigProfilerMatrixGeneratorFunc("class2", "GRCh37", project_path + 'subclasses'+ path_suffix +'/class2/', seqInfo=True)#, plot=True)
		print()
	except:
		pass
	try:
		print("Generating matrices for Class 3 mutations:")
		matrices = matGen.SigProfilerMatrixGeneratorFunc("class3", "GRCh37", project_path + 'subclasses'+ path_suffix +'/class3/', seqInfo=True)#, plot=True)
		print()
	except:
		pass
	try:
		print("Generating matrices for Class 1a mutations:")
		matrices = matGen.SigProfilerMatrixGeneratorFunc("class1a", "GRCh37", project_path + 'subclasses'+ path_suffix +'/class1a/', seqInfo=True)#, plot=True)
		print()
	except:
		pass
	try:	
		print("Generating matrices for Class 1b mutations:")
		matrices = matGen.SigProfilerMatrixGeneratorFunc("class1b", "GRCh37", project_path + 'subclasses'+ path_suffix +'/class1b/', seqInfo=True)#, plot=True)
		print()
	except:
		pass
	try:
		print("Generating matrices for Class 1c mutations:")
		matrices = matGen.SigProfilerMatrixGeneratorFunc("class1c", "GRCh37", project_path + 'subclasses'+ path_suffix +'/class1c/', seqInfo=True)#, plot=True)
		print()
	except:
		pass











