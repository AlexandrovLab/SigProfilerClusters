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
import pandas as pd
import matplotlib.patches as mplpatches
import pickle
import bisect


def plot96_same(
    matrix_path,
    matrix_path_clustered,
    matrix_path_nonClustered,
    sample,
    percentage,
    signature,
    panel1,
    panel3,
    panel5,
    fig,
    chrom=None,
):
    """
    Plots the SBS96 catalogs.

    Parameters:
                                     matrix_path	-> 	path to SBS96 matrix (string)
               matrix_path_clustered	-> 	path to SBS96 matrix of clustered mutations (string)
            matrix_path_nonClustered	->	path to SBS96 matrix of non-clustered mutations (string)
                                              sample	->	name of current sample (string)
                                      percentage	->	optional parameter to normalize plots (boolean; default=False)
                                       signature	-> 	optional parameter to generate plots using floats (boolean; default=False)
                                              panel1	->	axis information for the SBS96 plot of all mutations
                                              panel3	->	axis information for the SBS96 plot of clustered mutations
                                              panel5	->	axis information for the SBS96 plot of non-clustered mutations
                                                     fig	->	the matplotlib object for the current figure
                                               chrom	-> 	DO NOT USE

    Returns:
            None

    Outputs:
            Plots for all of the SBS96 axes
    """
    # if 'roman' in matplotlib.font_manager.weight_dict:
    # 	del matplotlib.font_manager.weight_dict['roman']
    # 	matplotlib.font_manager._rebuild()

    mutations = dict()
    total_count = []
    with open(matrix_path) as f:
        first_line = f.readline()
        samples = first_line.strip().split()
        samples = samples[1:]
        sample_index = samples.index(sample) + 1
        mutations[sample] = {
            "C>A": OrderedDict(),
            "C>G": OrderedDict(),
            "C>T": OrderedDict(),
            "T>A": OrderedDict(),
            "T>C": OrderedDict(),
            "T>G": OrderedDict(),
        }

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
    colors = [
        [3 / 256, 189 / 256, 239 / 256],
        [1 / 256, 1 / 256, 1 / 256],
        [228 / 256, 41 / 256, 38 / 256],
        [203 / 256, 202 / 256, 202 / 256],
        [162 / 256, 207 / 256, 99 / 256],
        [236 / 256, 199 / 256, 197 / 256],
    ]
    i = 0
    for key in mutations[sample]:
        for seq in mutations[sample][key]:
            xlabels.append(seq[0] + seq[2] + seq[6])
            if signature:
                if percentage:
                    panel1.bar(
                        x,
                        mutations[sample][key][seq] * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if mutations[sample][key][seq] * 100 > ymax:
                        ymax = mutations[sample][key][seq] * 100
                else:
                    panel1.bar(
                        x,
                        mutations[sample][key][seq] / total_count * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if mutations[sample][key][seq] / total_count * 100 > ymax:
                        ymax = mutations[sample][key][seq] / total_count * 100
            else:
                panel1.bar(
                    x,
                    mutations[sample][key][seq],
                    width=0.4,
                    color=colors[i],
                    align="center",
                    zorder=1000,
                )
                if mutations[sample][key][seq] > ymax:
                    ymax = mutations[sample][key][seq]
            x += 1
        i += 1

    x = 0.077
    y3 = 0.895
    y = int(ymax * 1.25)
    y2 = y + 2
    for i in range(0, 6, 1):
        panel1.add_patch(
            plt.Rectangle(
                (x, y3),
                0.088,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.09

    yText = y3 + 0.015
    plt.text(
        0.11,
        yText,
        "C>A",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2,
        yText,
        "C>G",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.288,
        yText,
        "C>T",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.376,
        yText,
        "T>A",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.466,
        yText,
        "T>C",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.554,
        yText,
        "T>G",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )

    while y % 4 != 0:
        y += 1
    ytick_offest = int(y / 4)
    if signature:
        ylabs = [
            0,
            round(ytick_offest, 1),
            round(ytick_offest * 2, 1),
            round(ytick_offest * 3, 1),
            round(ytick_offest * 4, 1),
        ]
        ylabels = [
            str(0),
            str(round(ytick_offest, 1)) + "%",
            str(round(ytick_offest * 2, 1)) + "%",
            str(round(ytick_offest * 3, 1)) + "%",
            str(round(ytick_offest * 4, 1)) + "%",
        ]
    else:
        ylabs = [0, ytick_offest, ytick_offest * 2, ytick_offest * 3, ytick_offest * 4]
        ylabels = [
            0,
            ytick_offest,
            ytick_offest * 2,
            ytick_offest * 3,
            ytick_offest * 4,
        ]

    labs = np.arange(0.375, 96.375, 1)
    panel1.set_xlim([0, 96])
    panel1.set_ylim([0, y])
    panel1.set_xticks(labs)
    panel1.set_yticks(ylabs)
    count = 0
    m = 0
    for i in range(0, 96, 1):
        plt.text(
            i / 177 + 0.075,
            0.052,
            xlabels[i][0],
            fontsize=6,
            color="black",
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            i / 177 + 0.075,
            0.06,
            xlabels[i][1],
            fontsize=6,
            color=colors[m],
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            fontweight="bold",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            i / 177 + 0.075,
            0.068,
            xlabels[i][2],
            fontsize=6,
            color="black",
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            transform=plt.gcf().transFigure,
        )
        count += 1
        if count == 16:
            count = 0
            m += 1

    panel1.set_yticklabels(ylabels, fontsize=8, color="b")
    panel1.set_xlabel("")
    panel1.set_ylabel("")
    panel1.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
    panel1.set_ylabel(
        "All Mutations", fontname="Times New Roman", fontsize=15, fontweight="bold"
    )
    panel1.tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=True,
        labelleft=True,
        right=True,
        labelright=False,
        top=False,
        labeltop=False,
        direction="in",
        length=10,
        colors="lightgray",
        width=1,
    )
    [i.set_color("black") for i in panel1.get_yticklabels()]

    #############################################################################################
    # Plot 96 bar plot (clustered)
    #############################################################################################
    mutations = dict()
    total_count = []
    with open(matrix_path_clustered) as f:
        first_line = f.readline()
        samples = first_line.strip().split()
        samples = samples[1:]
        sample_index = samples.index(sample) + 1
        mutations[sample] = {
            "C>A": OrderedDict(),
            "C>G": OrderedDict(),
            "C>T": OrderedDict(),
            "T>A": OrderedDict(),
            "T>C": OrderedDict(),
            "T>G": OrderedDict(),
        }

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
    colors = [
        [3 / 256, 189 / 256, 239 / 256],
        [1 / 256, 1 / 256, 1 / 256],
        [228 / 256, 41 / 256, 38 / 256],
        [203 / 256, 202 / 256, 202 / 256],
        [162 / 256, 207 / 256, 99 / 256],
        [236 / 256, 199 / 256, 197 / 256],
    ]
    i = 0
    for key in mutations[sample]:
        for seq in mutations[sample][key]:
            xlabels.append(seq[0] + seq[2] + seq[6])
            if signature:
                if percentage:
                    panel3.bar(
                        x,
                        mutations[sample][key][seq] * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if mutations[sample][key][seq] * 100 > ymax:
                        ymax = mutations[sample][key][seq] * 100
                else:
                    panel3.bar(
                        x,
                        mutations[sample][key][seq] / total_count * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if mutations[sample][key][seq] / total_count * 100 > ymax:
                        ymax = mutations[sample][key][seq] / total_count * 100
            else:
                panel3.bar(
                    x,
                    mutations[sample][key][seq],
                    width=0.4,
                    color=colors[i],
                    align="center",
                    zorder=1000,
                )
                if mutations[sample][key][seq] > ymax:
                    ymax = mutations[sample][key][seq]
            x += 1
        i += 1

    x = 0.077
    y3 = 0.598
    y = int(ymax * 1.25)
    y2 = y + 2
    for i in range(0, 6, 1):
        panel3.add_patch(
            plt.Rectangle(
                (x, y3),
                0.088,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.09

    yText = y3 + 0.015
    panel3.text(
        0.11,
        yText,
        "C>A",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel3.text(
        0.2,
        yText,
        "C>G",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel3.text(
        0.288,
        yText,
        "C>T",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel3.text(
        0.376,
        yText,
        "T>A",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel3.text(
        0.466,
        yText,
        "T>C",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel3.text(
        0.554,
        yText,
        "T>G",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )

    panel3.set_ylabel(
        "Clustered", fontname="Times New Roman", fontsize=15, fontweight="bold"
    )

    while y % 4 != 0:
        y += 1
    ytick_offest = int(y / 4)

    if signature:
        ylabs = [
            0,
            round(ytick_offest, 1),
            round(ytick_offest * 2, 1),
            round(ytick_offest * 3, 1),
            round(ytick_offest * 4, 1),
        ]
        ylabels = [
            str(0),
            str(round(ytick_offest, 1)) + "%",
            str(round(ytick_offest * 2, 1)) + "%",
            str(round(ytick_offest * 3, 1)) + "%",
            str(round(ytick_offest * 4, 1)) + "%",
        ]
    else:
        ylabs = [0, ytick_offest, ytick_offest * 2, ytick_offest * 3, ytick_offest * 4]
        ylabels = [
            0,
            ytick_offest,
            ytick_offest * 2,
            ytick_offest * 3,
            ytick_offest * 4,
        ]

    labs = np.arange(0.375, 96.375, 1)

    panel3.set_xlim([0, 96])
    panel3.set_yticklabels(ylabels, fontsize=8, color="b")
    panel3.set_ylim([0, y])
    panel3.set_xticks(labs)
    panel3.set_yticks(ylabs)
    count = 0
    m = 0
    for i in range(0, 96, 1):
        plt.text(
            i / 177 + 0.075,
            0.3495,
            xlabels[i][0],
            fontsize=6,
            color="black",
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            i / 177 + 0.075,
            0.3575,
            xlabels[i][1],
            fontsize=6,
            color=colors[m],
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            fontweight="bold",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            i / 177 + 0.075,
            0.3655,
            xlabels[i][2],
            fontsize=6,
            color="black",
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            transform=plt.gcf().transFigure,
        )
        count += 1
        if count == 16:
            count = 0
            m += 1

    panel3.tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=True,
        labelleft=True,
        right=True,
        labelright=False,
        top=False,
        labeltop=False,
        direction="in",
        length=10,
        colors="lightgray",
        width=1,
    )

    panel3.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
    [i.set_color("black") for i in panel3.get_yticklabels()]

    #############################################################################################
    # Plot 96 bar plot (non-clustered)
    #############################################################################################
    mutations = dict()
    total_count = []
    with open(matrix_path_nonClustered) as f:
        first_line = f.readline()
        samples = first_line.strip().split()
        samples = samples[1:]
        sample_index = samples.index(sample) + 1
        mutations[sample] = {
            "C>A": OrderedDict(),
            "C>G": OrderedDict(),
            "C>T": OrderedDict(),
            "T>A": OrderedDict(),
            "T>C": OrderedDict(),
            "T>G": OrderedDict(),
        }

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
    colors = [
        [3 / 256, 189 / 256, 239 / 256],
        [1 / 256, 1 / 256, 1 / 256],
        [228 / 256, 41 / 256, 38 / 256],
        [203 / 256, 202 / 256, 202 / 256],
        [162 / 256, 207 / 256, 99 / 256],
        [236 / 256, 199 / 256, 197 / 256],
    ]
    i = 0
    for key in mutations[sample]:
        for seq in mutations[sample][key]:
            xlabels.append(seq[0] + seq[2] + seq[6])
            if signature:
                if percentage:
                    panel5.bar(
                        x,
                        mutations[sample][key][seq] * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if mutations[sample][key][seq] * 100 > ymax:
                        ymax = mutations[sample][key][seq] * 100
                else:
                    panel5.bar(
                        x,
                        mutations[sample][key][seq] / total_count * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if mutations[sample][key][seq] / total_count * 100 > ymax:
                        ymax = mutations[sample][key][seq] / total_count * 100
            else:
                panel5.bar(
                    x,
                    mutations[sample][key][seq],
                    width=0.4,
                    color=colors[i],
                    align="center",
                    zorder=1000,
                )
                if mutations[sample][key][seq] > ymax:
                    ymax = mutations[sample][key][seq]
            x += 1
        i += 1

    x = 0.077
    y3 = 0.301
    y = int(ymax * 1.25)
    y2 = y + 2
    for i in range(0, 6, 1):
        panel5.add_patch(
            plt.Rectangle(
                (x, y3),
                0.088,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.09

    yText = y3 + 0.015
    panel5.text(
        0.11,
        yText,
        "C>A",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel5.text(
        0.2,
        yText,
        "C>G",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel5.text(
        0.288,
        yText,
        "C>T",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel5.text(
        0.376,
        yText,
        "T>A",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel5.text(
        0.466,
        yText,
        "T>C",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )
    panel5.text(
        0.554,
        yText,
        "T>G",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
        transform=plt.gcf().transFigure,
    )

    panel5.set_ylabel(
        "Non-Clustered", fontname="Times New Roman", fontsize=15, fontweight="bold"
    )

    while y % 4 != 0:
        y += 1
    ytick_offest = int(y / 4)

    if signature:
        ylabs = [
            0,
            round(ytick_offest, 1),
            round(ytick_offest * 2, 1),
            round(ytick_offest * 3, 1),
            round(ytick_offest * 4, 1),
        ]
        ylabels = [
            str(0),
            str(round(ytick_offest, 1)) + "%",
            str(round(ytick_offest * 2, 1)) + "%",
            str(round(ytick_offest * 3, 1)) + "%",
            str(round(ytick_offest * 4, 1)) + "%",
        ]
    else:
        ylabs = [0, ytick_offest, ytick_offest * 2, ytick_offest * 3, ytick_offest * 4]
        ylabels = [
            0,
            ytick_offest,
            ytick_offest * 2,
            ytick_offest * 3,
            ytick_offest * 4,
        ]

    labs = np.arange(0.375, 96.375, 1)

    panel5.set_xlim([0, 96])
    panel5.set_yticklabels(ylabels, fontsize=8, color="b")
    panel5.set_ylim([0, y])
    panel5.set_xticks(labs)
    panel5.set_yticks(ylabs)
    panel5.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
    count = 0
    m = 0
    for i in range(0, 96, 1):
        plt.text(
            i / 177 + 0.075,
            0.647,
            xlabels[i][0],
            fontsize=6,
            color="black",
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            i / 177 + 0.075,
            0.655,
            xlabels[i][1],
            fontsize=6,
            color=colors[m],
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            fontweight="bold",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            i / 177 + 0.075,
            0.663,
            xlabels[i][2],
            fontsize=6,
            color="black",
            rotation="vertical",
            verticalalignment="center",
            fontname="Courier New",
            transform=plt.gcf().transFigure,
        )
        count += 1
        if count == 16:
            count = 0
            m += 1

    panel5.tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=True,
        labelleft=True,
        right=True,
        labelright=False,
        top=False,
        labeltop=False,
        direction="in",
        length=10,
        colors="lightgray",
        width=1,
    )

    [i.set_color("black") for i in panel5.get_yticklabels()]
    if chrom:
        fig.suptitle(
            sample + "(" + chrom + ")", fontsize=30, fontname="Times New Roman"
        )
    else:
        fig.suptitle(sample, fontsize=30, fontname="Times New Roman")
    plt.text(
        0.008,
        0.55,
        "Mutation Counts",
        rotation="vertical",
        fontsize=20,
        fontweight="bold",
        fontname="Times New Roman",
        transform=plt.gcf().transFigure,
    )


def plotINDEL_same(
    matrix_path,
    matrix_path_clustered,
    matrix_path_nonClustered,
    sample,
    percentage,
    signature,
    panel1,
    panel3,
    panel5,
    fig,
):
    """
    Plots the ID83 catalogs.

    Parameters:
                                     matrix_path	-> 	path to SBS96 matrix (string)
               matrix_path_clustered	-> 	path to SBS96 matrix of clustered mutations (string)
            matrix_path_nonClustered	->	path to SBS96 matrix of non-clustered mutations (string)
                                              sample	->	name of current sample (string)
                                      percentage	->	optional parameter to normalize plots (boolean; default=False)
                                       signature	-> 	optional parameter to generate plots using floats (boolean; default=False)
                                              panel1	->	axis information for the SBS96 plot of all mutations
                                              panel3	->	axis information for the SBS96 plot of clustered mutations
                                              panel5	->	axis information for the SBS96 plot of non-clustered mutations
                                                     fig	->	the matplotlib object for the current figure

    Returns:
            None

    Outputs:
            Plots for all of the ID83 axes
    """
    # if 'roman' in matplotlib.font_manager.weight_dict:
    # 	del matplotlib.font_manager.weight_dict['roman']
    # 	matplotlib.font_manager._rebuild()

    indel_types = [
        "1:Del:C:1",
        "1:Del:C:2",
        "1:Del:C:3",
        "1:Del:C:4",
        "1:Del:C:5",
        "1:Del:C:6" "1:Del:T:1",
        "1:Del:T:2",
        "1:Del:T:3",
        "1:Del:T:4",
        "1:Del:T:5",
        "1:Del:T:6" "1:Ins:C:0",
        "1:Ins:C:1",
        "1:Ins:C:2",
        "1:Ins:C:3",
        "1:Ins:C:4",
        "1:Ins:C:5",
        "1:Ins:T:0",
        "1:Ins:T:1",
        "1:Ins:T:2",
        "1:Ins:T:3",
        "1:Ins:T:4",
        "1:Ins:T:5",
        # >1bp INDELS
        "2:Del:R:0",
        "2:Del:R:1",
        "2:Del:R:2",
        "2:Del:R:3",
        "2:Del:R:4",
        "2:Del:R:5",
        "3:Del:R:0",
        "3:Del:R:1",
        "3:Del:R:2",
        "3:Del:R:3",
        "3:Del:R:4",
        "3:Del:R:5",
        "4:Del:R:0",
        "4:Del:R:1",
        "4:Del:R:2",
        "4:Del:R:3",
        "4:Del:R:4",
        "4:Del:R:5",
        "5:Del:R:0",
        "5:Del:R:1",
        "5:Del:R:2",
        "5:Del:R:3",
        "5:Del:R:4",
        "5:Del:R:5",
        "2:Ins:R:0",
        "2:Ins:R:1",
        "2:Ins:R:2",
        "2:Ins:R:3",
        "2:Ins:R:4",
        "2:Ins:R:5",
        "3:Ins:R:0",
        "3:Ins:R:1",
        "3:Ins:R:2",
        "3:Ins:R:3",
        "3:Ins:R:4",
        "3:Ins:R:5",
        "4:Ins:R:0",
        "4:Ins:R:1",
        "4:Ins:R:2",
        "4:Ins:R:3",
        "4:Ins:R:4",
        "4:Ins:R:5",
        "5:Ins:R:0",
        "5:Ins:R:1",
        "5:Ins:R:2",
        "5:Ins:R:3",
        "5:Ins:R:4",
        "5:Ins:R:5",
        # MicroHomology INDELS
        "2:Del:M:1",
        "3:Del:M:1",
        "3:Del:M:2",
        "4:Del:M:1",
        "4:Del:M:2",
        "4:Del:M:3",
        "5:Del:M:1",
        "5:Del:M:2",
        "5:Del:M:3",
        "5:Del:M:4",
        "5:Del:M:5",
        "2:Ins:M:1",
        "3:Ins:M:1",
        "3:Ins:M:2",
        "4:Ins:M:1",
        "4:Ins:M:2",
        "4:Ins:M:3",
        "5:Ins:M:1",
        "5:Ins:M:2",
        "5:Ins:M:3",
        "5:Ins:M:4",
        "5:Ins:M:5",
    ]

    mutations = dict()
    with open(matrix_path) as f:
        first_line = f.readline()
        samples = first_line.strip().split()
        samples = samples[1:]
        sample_index = samples.index(sample) + 1
        mutations[sample] = {
            "1DelC": [0, 0, 0, 0, 0, 0],
            "1DelT": [0, 0, 0, 0, 0, 0],
            "1InsC": [0, 0, 0, 0, 0, 0],
            "1InsT": [0, 0, 0, 0, 0, 0],
            "2DelR": [0, 0, 0, 0, 0, 0],
            "3DelR": [0, 0, 0, 0, 0, 0],
            "4DelR": [0, 0, 0, 0, 0, 0],
            "5DelR": [0, 0, 0, 0, 0, 0],
            "2InsR": [0, 0, 0, 0, 0, 0],
            "3InsR": [0, 0, 0, 0, 0, 0],
            "4InsR": [0, 0, 0, 0, 0, 0],
            "5InsR": [0, 0, 0, 0, 0, 0],
            "2DelM": [0],
            "3DelM": [0, 0],
            "4DelM": [0, 0, 0],
            "5DelM": [0, 0, 0, 0, 0],
        }

        for lines in f:
            line = lines.strip().split()
            categories = line[0].split(":")
            mut_type = categories[0] + categories[1] + categories[2]
            repeat_size = int(categories[3])
            if categories[2] == "M":
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
    colors = [
        [253 / 256, 190 / 256, 111 / 256],
        [255 / 256, 128 / 256, 2 / 256],
        [176 / 256, 221 / 256, 139 / 256],
        [54 / 256, 161 / 256, 46 / 256],
        [253 / 256, 202 / 256, 181 / 256],
        [252 / 256, 138 / 256, 106 / 256],
        [241 / 256, 68 / 256, 50 / 256],
        [188 / 256, 25 / 256, 26 / 256],
        [208 / 256, 225 / 256, 242 / 256],
        [148 / 256, 196 / 256, 223 / 256],
        [74 / 256, 152 / 256, 201 / 256],
        [23 / 256, 100 / 256, 171 / 256],
        [226 / 256, 226 / 256, 239 / 256],
        [182 / 256, 182 / 256, 216 / 256],
        [134 / 256, 131 / 256, 189 / 256],
        [98 / 256, 64 / 256, 155 / 256],
    ]

    i = 0
    for key in mutations[sample]:
        l = 1
        for seq in mutations[sample][key]:
            xlabels.append(l)
            if signature:
                if percentage:
                    panel1.bar(
                        x,
                        seq * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if seq * 100 > ymax:
                        ymax = seq * 100

                else:
                    panel1.bar(
                        x,
                        seq / total_count * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if seq / total_count * 100 > ymax:
                        ymax = seq / total_count * 100
            else:
                panel1.bar(
                    x, seq, width=0.4, color=colors[i], align="center", zorder=1000
                )
                if seq > ymax:
                    ymax = seq
            x += 1
            l += 1
        i += 1

    x = 0.0757
    y_top = 0.895
    y_bottom = 0.656
    y = int(ymax * 1.25)
    y2 = y + 2
    for i in range(0, 12, 1):
        panel1.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.037,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.037,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.0393

    panel1.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.0035,
            0.01,
            facecolor=colors[12],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel1.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.0035,
            0.01,
            facecolor=colors[12],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0058
    panel1.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.011,
            0.01,
            facecolor=colors[13],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel1.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.011,
            0.01,
            facecolor=colors[13],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0133
    panel1.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.018,
            0.01,
            facecolor=colors[14],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel1.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.018,
            0.01,
            facecolor=colors[14],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0203
    panel1.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.03,
            0.01,
            facecolor=colors[15],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel1.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.03,
            0.01,
            facecolor=colors[15],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )

    yText = y_top + 0.0015
    plt.text(
        0.092,
        yText,
        "C",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.1313,
        yText,
        "T",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.1706,
        yText,
        "C",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2099,
        yText,
        "T",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2492,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2885,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.3278,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.3671,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.4064,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.4457,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.485,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5243,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5467,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5565,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.573,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5977,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )

    yText_labels_top = yText + 0.015
    yText_labels_bottom = y_bottom + 0.002
    yText_labels_bottom_sec = yText_labels_bottom - 0.015

    plt.text(
        0.09,
        yText_labels_top,
        "1bp Deletion",
        fontsize=7,
        fontname="Times New Roman",
        weight="bold",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.167,
        yText_labels_top,
        "1bp Insertion",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.266,
        yText_labels_top,
        ">1bp Deletion at Repeats\n      (Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.42,
        yText_labels_top,
        ">1bp Insertions at Repeats\n       (Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.55,
        yText_labels_top,
        " Mircohomology\n(Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )

    plt.text(
        0.079,
        yText_labels_bottom_sec,
        "Homopolymer Length",
        fontsize=6.5,
        fontname="Times New Roman",
        weight="bold",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.156,
        yText_labels_bottom_sec,
        "Homopolymer Length",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.275,
        yText_labels_bottom_sec,
        "Number of Repeat Units",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.43,
        yText_labels_bottom_sec,
        "Number of Repeat Units",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5475,
        yText_labels_bottom_sec,
        "Mircohomology Length",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )

    x = 0.0767
    for i in range(0, 8, 1):
        if i != 2 and i != 3:
            if i == 1 or i == 7:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5  6+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5  6+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
        else:
            if i == 3:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "0  1  2  3  4  5+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "0  1  2  3  4  5+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

        x += 0.03925

    for i in range(0, 4, 1):
        if i == 3:
            plt.text(
                x,
                yText_labels_bottom,
                "0  1  2  3  4  5+",
                fontsize=5.7,
                fontweight="bold",
                fontname="Times New Roman",
                color="white",
                transform=plt.gcf().transFigure,
            )
        else:
            plt.text(
                x,
                yText_labels_bottom,
                "0  1  2  3  4  5+",
                fontsize=5.7,
                fontweight="bold",
                fontname="Times New Roman",
                color="black",
                transform=plt.gcf().transFigure,
            )

        x += 0.03925

    plt.text(
        x,
        yText_labels_bottom,
        "1",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.006
    plt.text(
        x,
        yText_labels_bottom,
        "1  2",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.0135
    plt.text(
        x,
        yText_labels_bottom,
        "1  2  3",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.02
    plt.text(
        x,
        yText_labels_bottom,
        "1  2  3  4  5+",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )

    while y % 4 != 0:
        y += 1
    if signature:
        ytick_offest = int(y / 4)
        ylabs = [
            0,
            round(ytick_offest, 1),
            round(ytick_offest * 2, 1),
            round(ytick_offest * 3, 1),
            round(ytick_offest * 4, 1),
        ]
        ylabels = [
            str(0),
            str(round(ytick_offest, 1)) + "%",
            str(round(ytick_offest * 2, 1)) + "%",
            str(round(ytick_offest * 3, 1)) + "%",
            str(round(ytick_offest * 4, 1)) + "%",
        ]

    else:
        ytick_offest = int(y / 4)
        ylabs = [0, ytick_offest, ytick_offest * 2, ytick_offest * 3, ytick_offest * 4]
        ylabels = [
            0,
            ytick_offest,
            ytick_offest * 2,
            ytick_offest * 3,
            ytick_offest * 4,
        ]

    labs = np.arange(0.375, 83.375, 1)
    panel1.set_xlim([0, 83])
    panel1.set_yticklabels(ylabels, fontsize=8, color="b")
    panel1.set_ylim([0, y])
    panel1.set_xticks(labs)
    panel1.set_yticks(ylabs)
    panel1.set_yticklabels(ylabels, fontsize=8)
    panel1.grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
    panel1.set_xlabel("")
    panel1.set_ylabel("")
    panel1.tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=False,
        labelleft=True,
        right=False,
        labelright=False,
        top=False,
        labeltop=False,
        direction="in",
        length=25,
        colors="gray",
        width=2,
    )

    [i.set_color("black") for i in panel1.get_yticklabels()]
    panel1.set_ylabel(
        "All Mutations", fontname="Times New Roman", fontsize=15, fontweight="bold"
    )

    ########################################################
    # Non-clustered
    ########################################################
    mutations = dict()
    with open(matrix_path_clustered) as f:
        first_line = f.readline()
        samples = first_line.strip().split()
        samples = samples[1:]
        sample_index = samples.index(sample) + 1
        mutations[sample] = {
            "1DelC": [0, 0, 0, 0, 0, 0],
            "1DelT": [0, 0, 0, 0, 0, 0],
            "1InsC": [0, 0, 0, 0, 0, 0],
            "1InsT": [0, 0, 0, 0, 0, 0],
            "2DelR": [0, 0, 0, 0, 0, 0],
            "3DelR": [0, 0, 0, 0, 0, 0],
            "4DelR": [0, 0, 0, 0, 0, 0],
            "5DelR": [0, 0, 0, 0, 0, 0],
            "2InsR": [0, 0, 0, 0, 0, 0],
            "3InsR": [0, 0, 0, 0, 0, 0],
            "4InsR": [0, 0, 0, 0, 0, 0],
            "5InsR": [0, 0, 0, 0, 0, 0],
            "2DelM": [0],
            "3DelM": [0, 0],
            "4DelM": [0, 0, 0],
            "5DelM": [0, 0, 0, 0, 0],
        }

        for lines in f:
            line = lines.strip().split()
            categories = line[0].split(":")
            mut_type = categories[0] + categories[1] + categories[2]
            repeat_size = int(categories[3])
            if categories[2] == "M":
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
                    panel3.bar(
                        x,
                        seq * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if seq * 100 > ymax:
                        ymax = seq * 100

                else:
                    panel3.bar(
                        x,
                        seq / total_count * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if seq / total_count * 100 > ymax:
                        ymax = seq / total_count * 100
            else:
                panel3.bar(
                    x, seq, width=0.4, color=colors[i], align="center", zorder=1000
                )
                if seq > ymax:
                    ymax = seq
            x += 1
            l += 1
        i += 1

    x = 0.0757
    y_top = 0.597
    y_bottom = 0.36
    y = int(ymax * 1.25)
    y2 = y + 2
    for i in range(0, 12, 1):
        panel1.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.037,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.037,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.0393

    panel3.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.0035,
            0.01,
            facecolor=colors[12],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel3.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.0035,
            0.01,
            facecolor=colors[12],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0058
    panel3.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.011,
            0.01,
            facecolor=colors[13],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel3.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.011,
            0.01,
            facecolor=colors[13],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0133
    panel3.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.018,
            0.01,
            facecolor=colors[14],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel3.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.018,
            0.01,
            facecolor=colors[14],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0203
    panel3.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.03,
            0.01,
            facecolor=colors[15],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel3.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.03,
            0.01,
            facecolor=colors[15],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )

    yText = y_top + 0.00138
    plt.text(
        0.092,
        yText,
        "C",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.1313,
        yText,
        "T",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.1706,
        yText,
        "C",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2099,
        yText,
        "T",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2492,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2885,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.3278,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.3671,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.4064,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.4457,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.485,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5243,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5467,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5565,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.573,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5977,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )

    yText_labels_top = yText + 0.015
    yText_labels_bottom = y_bottom + 0.002
    yText_labels_bottom_sec = yText_labels_bottom - 0.015

    plt.text(
        0.09,
        yText_labels_top,
        "1bp Deletion",
        fontsize=7,
        fontname="Times New Roman",
        weight="bold",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.167,
        yText_labels_top,
        "1bp Insertion",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.266,
        yText_labels_top,
        ">1bp Deletion at Repeats\n      (Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.42,
        yText_labels_top,
        ">1bp Insertions at Repeats\n       (Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.55,
        yText_labels_top,
        " Mircohomology\n(Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )

    plt.text(
        0.079,
        yText_labels_bottom_sec,
        "Homopolymer Length",
        fontsize=6.5,
        fontname="Times New Roman",
        weight="bold",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.156,
        yText_labels_bottom_sec,
        "Homopolymer Length",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.275,
        yText_labels_bottom_sec,
        "Number of Repeat Units",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.43,
        yText_labels_bottom_sec,
        "Number of Repeat Units",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5475,
        yText_labels_bottom_sec,
        "Mircohomology Length",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )

    x = 0.0767
    for i in range(0, 8, 1):
        if i != 2 and i != 3:
            if i == 1 or i == 7:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5  6+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5  6+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
        else:
            if i == 3:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "0  1  2  3  4  5+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "0  1  2  3  4  5+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

        x += 0.03925

    for i in range(0, 4, 1):
        if i == 3:
            plt.text(
                x,
                yText_labels_bottom,
                "0  1  2  3  4  5+",
                fontsize=5.7,
                fontweight="bold",
                fontname="Times New Roman",
                color="white",
                transform=plt.gcf().transFigure,
            )
        else:
            plt.text(
                x,
                yText_labels_bottom,
                "0  1  2  3  4  5+",
                fontsize=5.7,
                fontweight="bold",
                fontname="Times New Roman",
                color="black",
                transform=plt.gcf().transFigure,
            )
        x += 0.03925
    plt.text(
        x,
        yText_labels_bottom,
        "1",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.006
    plt.text(
        x,
        yText_labels_bottom,
        "1  2",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.0135
    plt.text(
        x,
        yText_labels_bottom,
        "1  2  3",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.02
    plt.text(
        x,
        yText_labels_bottom,
        "1  2  3  4  5+",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )

    while y % 4 != 0:
        y += 1
    if signature:
        ytick_offest = int(y / 4)
        ylabs = [
            0,
            round(ytick_offest, 1),
            round(ytick_offest * 2, 1),
            round(ytick_offest * 3, 1),
            round(ytick_offest * 4, 1),
        ]
        ylabels = [
            str(0),
            str(round(ytick_offest, 1)) + "%",
            str(round(ytick_offest * 2, 1)) + "%",
            str(round(ytick_offest * 3, 1)) + "%",
            str(round(ytick_offest * 4, 1)) + "%",
        ]

    else:
        ytick_offest = int(y / 4)
        ylabs = [0, ytick_offest, ytick_offest * 2, ytick_offest * 3, ytick_offest * 4]
        ylabels = [
            0,
            ytick_offest,
            ytick_offest * 2,
            ytick_offest * 3,
            ytick_offest * 4,
        ]

    panel3.set_xlim([0, 83])
    panel3.set_yticklabels(ylabels, fontsize=8, color="b")
    panel3.set_ylim([0, y])
    panel3.set_xticks(labs)
    panel3.set_yticks(ylabs)
    panel3.set_yticklabels(ylabels, fontsize=8)
    panel3.grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
    panel3.set_xlabel("")
    panel3.set_ylabel("")
    panel3.tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=False,
        labelleft=True,
        right=False,
        labelright=False,
        top=False,
        labeltop=False,
        direction="in",
        length=25,
        colors="gray",
        width=2,
    )

    [i.set_color("black") for i in panel3.get_yticklabels()]
    panel3.set_ylabel(
        "Clustered", fontname="Times New Roman", fontsize=15, fontweight="bold"
    )

    ########################################################
    # Non-clustered
    ########################################################
    mutations = dict()
    with open(matrix_path_nonClustered) as f:
        first_line = f.readline()
        samples = first_line.strip().split()
        samples = samples[1:]
        try:
            sample_index = samples.index(sample) + 1
        except:
            pass
        mutations[sample] = {
            "1DelC": [0, 0, 0, 0, 0, 0],
            "1DelT": [0, 0, 0, 0, 0, 0],
            "1InsC": [0, 0, 0, 0, 0, 0],
            "1InsT": [0, 0, 0, 0, 0, 0],
            "2DelR": [0, 0, 0, 0, 0, 0],
            "3DelR": [0, 0, 0, 0, 0, 0],
            "4DelR": [0, 0, 0, 0, 0, 0],
            "5DelR": [0, 0, 0, 0, 0, 0],
            "2InsR": [0, 0, 0, 0, 0, 0],
            "3InsR": [0, 0, 0, 0, 0, 0],
            "4InsR": [0, 0, 0, 0, 0, 0],
            "5InsR": [0, 0, 0, 0, 0, 0],
            "2DelM": [0],
            "3DelM": [0, 0],
            "4DelM": [0, 0, 0],
            "5DelM": [0, 0, 0, 0, 0],
        }

        for lines in f:
            line = lines.strip().split()
            categories = line[0].split(":")
            mut_type = categories[0] + categories[1] + categories[2]
            repeat_size = int(categories[3])
            if categories[2] == "M":
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
                    panel5.bar(
                        x,
                        seq * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if seq * 100 > ymax:
                        ymax = seq * 100

                else:
                    panel5.bar(
                        x,
                        seq / total_count * 100,
                        width=0.4,
                        color=colors[i],
                        align="center",
                        zorder=1000,
                    )
                    if seq / total_count * 100 > ymax:
                        ymax = seq / total_count * 100
            else:
                panel5.bar(
                    x, seq, width=0.4, color=colors[i], align="center", zorder=1000
                )
                if seq > ymax:
                    ymax = seq
            x += 1
            l += 1
        i += 1

    x = 0.0757
    y_top = 0.3
    y_bottom = 0.0625
    y = int(ymax * 1.25)
    y2 = y + 2
    for i in range(0, 12, 1):
        panel5.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.037,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel5.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.037,
                0.01,
                facecolor=colors[i],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.0393

    panel5.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.0035,
            0.01,
            facecolor=colors[12],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel5.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.0035,
            0.01,
            facecolor=colors[12],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0058
    panel5.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.011,
            0.01,
            facecolor=colors[13],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel5.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.011,
            0.01,
            facecolor=colors[13],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0133
    panel5.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.018,
            0.01,
            facecolor=colors[14],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel5.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.018,
            0.01,
            facecolor=colors[14],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    x += 0.0203
    panel5.add_patch(
        plt.Rectangle(
            (x, y_top),
            0.03,
            0.01,
            facecolor=colors[15],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )
    panel5.add_patch(
        plt.Rectangle(
            (x, y_bottom),
            0.03,
            0.01,
            facecolor=colors[15],
            clip_on=False,
            transform=plt.gcf().transFigure,
        )
    )

    yText = y_top + 0.00138
    plt.text(
        0.092,
        yText,
        "C",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.1313,
        yText,
        "T",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.1706,
        yText,
        "C",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2099,
        yText,
        "T",
        fontsize=7,
        fontname="Times New Roman",
        fontweight="bold",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2492,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.2885,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.3278,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.3671,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.4064,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.4457,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.485,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5243,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5467,
        yText,
        "2",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5565,
        yText,
        "3",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.573,
        yText,
        "4",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5977,
        yText,
        "5+",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )

    yText_labels_top = yText + 0.015
    yText_labels_bottom = y_bottom + 0.002
    yText_labels_bottom_sec = yText_labels_bottom - 0.015

    plt.text(
        0.09,
        yText_labels_top,
        "1bp Deletion",
        fontsize=7,
        fontname="Times New Roman",
        weight="bold",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.167,
        yText_labels_top,
        "1bp Insertion",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.266,
        yText_labels_top,
        ">1bp Deletion at Repeats\n      (Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.42,
        yText_labels_top,
        ">1bp Insertions at Repeats\n       (Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.55,
        yText_labels_top,
        " Mircohomology\n(Deletion Length)",
        fontsize=7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )

    plt.text(
        0.079,
        yText_labels_bottom_sec,
        "Homopolymer Length",
        fontsize=6.5,
        fontname="Times New Roman",
        weight="bold",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.156,
        yText_labels_bottom_sec,
        "Homopolymer Length",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.275,
        yText_labels_bottom_sec,
        "Number of Repeat Units",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.43,
        yText_labels_bottom_sec,
        "Number of Repeat Units",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.5475,
        yText_labels_bottom_sec,
        "Mircohomology Length",
        fontsize=6.5,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )

    x = 0.0767
    for i in range(0, 8, 1):
        if i != 2 and i != 3:
            if i == 1 or i == 7:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5  6+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5  6+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
        else:
            if i == 3:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "0  1  2  3  4  5+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="white",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "0  1  2  3  4  5+",
                    fontsize=5.7,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

        x += 0.03925

    for i in range(0, 4, 1):
        if i == 3:
            plt.text(
                x,
                yText_labels_bottom,
                "0  1  2  3  4  5+",
                fontsize=5.7,
                fontweight="bold",
                fontname="Times New Roman",
                color="white",
                transform=plt.gcf().transFigure,
            )
        else:
            plt.text(
                x,
                yText_labels_bottom,
                "0  1  2  3  4  5+",
                fontsize=5.7,
                fontweight="bold",
                fontname="Times New Roman",
                color="black",
                transform=plt.gcf().transFigure,
            )
        x += 0.03925
    plt.text(
        x,
        yText_labels_bottom,
        "1",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.006
    plt.text(
        x,
        yText_labels_bottom,
        "1  2",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.0135
    plt.text(
        x,
        yText_labels_bottom,
        "1  2  3",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="black",
        transform=plt.gcf().transFigure,
    )
    x += 0.02
    plt.text(
        x,
        yText_labels_bottom,
        "1  2  3  4  5+",
        fontsize=5.7,
        fontweight="bold",
        fontname="Times New Roman",
        color="white",
        transform=plt.gcf().transFigure,
    )

    while y % 4 != 0:
        y += 1
    if signature:
        ytick_offest = int(y / 4)
        ylabs = [
            0,
            round(ytick_offest, 1),
            round(ytick_offest * 2, 1),
            round(ytick_offest * 3, 1),
            round(ytick_offest * 4, 1),
        ]
        ylabels = [
            str(0),
            str(round(ytick_offest, 1)) + "%",
            str(round(ytick_offest * 2, 1)) + "%",
            str(round(ytick_offest * 3, 1)) + "%",
            str(round(ytick_offest * 4, 1)) + "%",
        ]

    else:
        ytick_offest = int(y / 4)
        ylabs = [0, ytick_offest, ytick_offest * 2, ytick_offest * 3, ytick_offest * 4]
        ylabels = [
            0,
            ytick_offest,
            ytick_offest * 2,
            ytick_offest * 3,
            ytick_offest * 4,
        ]

    panel5.set_xlim([0, 83])
    panel5.set_yticklabels(ylabels, fontsize=8, color="b")
    panel5.set_ylim([0, y])
    panel5.set_xticks(labs)
    panel5.set_yticks(ylabs)
    panel5.set_yticklabels(ylabels, fontsize=8)
    panel5.grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)

    fig.suptitle(sample, fontsize=30, fontname="Times New Roman")
    plt.text(
        0.008,
        0.55,
        "Mutation Counts",
        rotation="vertical",
        fontsize=20,
        fontweight="bold",
        fontname="Times New Roman",
        transform=plt.gcf().transFigure,
    )
    panel5.tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=False,
        labelleft=True,
        right=False,
        labelright=False,
        top=False,
        labeltop=False,
        direction="in",
        length=25,
        colors="gray",
        width=2,
    )
    panel5.set_ylabel(
        "Non-Clustered", fontname="Times New Roman", fontsize=15, fontweight="bold"
    )
    [i.set_color("black") for i in panel5.get_yticklabels()]


def rainfall(
    chrom_based_IMD,
    project,
    project_path,
    chrom_path,
    chromLengths,
    centromeres,
    contexts,
    includedVAFs,
    includedCCFs,
    correction=True,
    windowSize=10000000,
    bedRanges=None,
):
    """
    Generates rainfall plots when subClassify is True.

    Parameters:
                    chrom_based_IMD	->	optional parameter for generating plots w.r.t chromosome-based IMDs (boolea; default=False)
                                    project	->	user provided project name (string)
                       project_path	-> 	the directory for the given project (string)
                             chrom_path	->	path to the reference genome chromosome files used for SigProfilerMatrixGenerator (string)
                             correction	-> 	optional parameter that corrects for mutational density (boolean; default=True)
                             windowSize	->	the window size used to calculate the mutation densities across the genome (integer; default=None)

    Returns:
            None

    Outputs:
            Rainfall plots for each sample.
    """
    aggregated = False
    chrom_based = False
    windowSize = 10000000
    path_suffix = ""
    if chrom_based_IMD:
        path_suffix = "_chrom"
    path_suffix2 = ""
    if correction:
        path_suffix2 += "_corrected"

    projectPath_parent2 = project_path
    projectPath_parent = project_path
    if contexts != "ID":
        projectPath = (
            project_path
            + "output/vcf_files"
            + path_suffix2
            + "/"
            + project
            + "_clustered/subclasses"
            + path_suffix
            + "/"
        )
    else:
        projectPath = project_path + "output/vcf_files" + path_suffix2 + "/"
    outputPath = projectPath_parent + "output/plots/"

    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    plot_suffix = project
    if chrom_based:
        plot_suffix += "_chrom_based"

    if chrom_based_IMD:
        plot_suffix += "_chrom"

    if correction:
        plot_suffix += "_corrected"

    if contexts != "ID":
        classes = [
            "Class IA",
            "Class II",
            "Class IB",
            "Class IC",
            "Class III",
            "Non-clust",
        ]
        mutationsPath = [
            projectPath + "class1a/" + project + "_clustered_class1a.txt",
            projectPath + "class2/" + project + "_clustered_class2.txt",
            projectPath + "class1b/" + project + "_clustered_class1b.txt",
            projectPath + "class1c/" + project + "_clustered_class1c.txt",
            projectPath + "class3/" + project + "_clustered_class3.txt",
            projectPath_parent
            + "output/vcf_files"
            + path_suffix2
            + "/"
            + project
            + "_nonClustered/SNV/"
            + project
            + "_nonClustered.txt",
        ]
    else:
        classes = ["Clust", "Non-clust"]
        mutationsPath = [
            projectPath + project + "_clustered/INDEL/" + project + "_clustered.txt",
            projectPath
            + project
            + "_nonClustered/INDEL/"
            + project
            + "_nonClustered.txt",
        ]

    with open(
        projectPath_parent2 + "output/simulations/data/imds" + path_suffix + ".pickle",
        "rb",
    ) as handle:
        imdsData = pickle.load(handle)

    with open(
        projectPath_parent2 + "output/simulations/data/imds.pickle", "rb"
    ) as handle:
        imdsDataSample = pickle.load(handle)
    if correction:
        with open(
            projectPath_parent2 + "output/simulations/data/imds_corrected.pickle", "rb"
        ) as handle:
            imdsDataSample_corrected = pickle.load(handle)

    if bedRanges:
        bedFile = pd.read_csv(bedRanges, sep="\t", header=0, index_col=3)

    headerFields = [
        "project",
        "samples",
        "ID",
        "genome",
        "mutType",
        "chr",
        "start",
        "end",
        "ref",
        "alt",
        "mutClass",
        "IMDplot",
        "group",
        "IMD",
        "vaf",
        "class",
        "failedReason",
    ]
    if not includedVAFs and not includedCCFs:
        headerFields = [
            "project",
            "samples",
            "ID",
            "genome",
            "mutType",
            "chr",
            "start",
            "end",
            "ref",
            "alt",
            "mutClass",
            "group",
            "IMDplot",
            "IMD",
            "class",
            "failedReason",
        ]
    # centromeres = {'GRCh38':{'1': [122026460,125184587],'10': [39686683,41593521],'11':[51078349,54425074],'12':[34769408,37185252],'13': [16000001,18051248],
    # 			'14': [16000001,18173523],'15': [17000001,19725254],'16': [36311159,38280682],'17': [22813680,26885980],'18': [15460900,20861206],
    # 			'19': [24498981,27190874],'2': [92188146,94090557],'20': [26436233,30038348],'21': [10864561,12915808],'22': [12954789,15054318],
    # 			'3': [90772459,93655574],'4': [49708101,51743951],'5': [46485901,50059807],'6': [58553889,59829934],'7': [58169654,60828234],
    # 			'8': [44033745,45877265],'9': [43236168,45518558],'X': [58605580,62412542],'Y': [10316945,10544039]},

    # 			'GRCh37':{'1': [121535434,124535434],'10': [39254935,42254935],'11':[51644205,54644205],'12':[34856694,37856694],'13': [16000000,19000000],
    # 			'14': [16000000,19000000],'15': [17000000,20000000],'16': [35335801,38335801],'17': [22263006,25263006],'18': [15460898,18460898],
    # 			'19': [24681782,27681782],'2': [92326171,95326171],'20': [26369569,29369569],'21': [11288129,14288129],'22': [13000000,16000000],
    # 			'3': [90504854,93504854],'4': [49660117,52660117],'5': [46405641,49405641],'6': [58830166,61830166],'7': [58054331,61054331],
    # 			'8': [43838887,46838887],'9': [47367679,50367679],'X': [58632012,61632012],'Y': [10316945,10544039]},

    # 			'mm10':{'1': [110000,3000000],'10': [110000,3000000],'11':[110000,3000000],'12':[110000,3000000],'13': [110000,3000000],
    # 			'14': [110000,3000000],'15': [110000,3000000],'16': [110000,3000000],'17': [110000,3000000],'18': [110000,3000000],
    # 			'19': [110000,3000000],'2': [110000,3000000],'3': [110000,3000000],'4': [110000,3000000],'5': [110000,3000000],
    # 			'6': [110000,3000000],'7': [110000,3000000],'8': [110000,3000000],'9': [110000,3000000],'X': [110000,3000000],'Y': [110000,3000000]},

    # 			'mm9':{'1': [0,3000000],'10': [0,3000000],'11':[0,3000000],'12':[0,3000000],'13': [0,3000000],
    # 			'14': [0,3000000],'15': [0,3000000],'16': [0,3000000],'17': [0,3000000],'18': [0,3000000],
    # 			'19': [0,3000000],'2': [0,3000000],'3': [0,3000000],'4': [0,3000000],'5': [0,3000000],
    # 			'6': [0,3000000],'7': [0,3000000],'8': [0,3000000],'9': [0,3000000],'X': [0,3000000],'Y': [0,3000000]}}

    chromLengths = {
        "GRCh37": {
            "1": 249250619,
            "2": 243199368,
            "3": 198022427,
            "4": 191154274,
            "5": 180915259,
            "6": 171115063,
            "7": 159138662,
            "8": 146364019,
            "9": 141213427,
            "10": 135534745,
            "11": 135006514,
            "12": 133851894,
            "13": 115169878,
            "14": 107349537,
            "15": 102531391,
            "16": 90354750,
            "17": 81195205,
            "18": 78077248,
            "19": 59128982,
            "20": 63025519,
            "21": 48129895,
            "22": 51304566,
            "X": 155270559,
            "Y": 59373566,
        }
    }

    # Gather all mutations based upon user input
    samplesSet = set()
    samples = {}
    allMutations = {}
    if classes[0] == "Class III":
        mutations = pd.read_csv(
            mutationsPath[0], sep="\t", names=headerFields, header=0
        )
    elif classes[0] == "Non-clust":
        mutations = pd.read_csv(
            mutationsPath[0], sep="\t", names=headerFields[:-3], header=0
        )  # , skiprows=[0])
    elif classes[0] == "Simulation":
        mutations = pd.read_csv(
            mutationsPath[0],
            sep="\t",
            names=[
                "project",
                "ID",
                "sim",
                "genome",
                "chr",
                "start",
                "end",
                "strand",
                "placeHolder",
                "mutType",
                "ref",
                "alt1",
                "alt",
                "placeHolder2",
                "placeHolder3",
                "samples",
                "seq",
            ],
            header=None,
            skiprows=[0],
        )
    elif classes[0] == "Clust":
        mutations = pd.read_csv(
            mutationsPath[0],
            sep="\t",
            names=[
                "project",
                "samples",
                "ID",
                "genome",
                "mutType",
                "chr",
                "start",
                "end",
                "ref",
                "alt",
                "mutClass",
                "IMDplot",
                "IMD",
            ],
            header=0,
            skiprows=[0],
        )
    else:
        mutations = pd.read_csv(
            mutationsPath[0], sep="\t", names=headerFields[:-1], header=0
        )
    # genome = mutations.head(1).loc[0,'genome']
    genome = None
    mutations["chr"] = mutations["chr"].astype(str)
    mutations = mutations.set_index(["samples", "chr"])
    mutations = mutations.sort_index()
    allMutations[classes[0]] = mutations
    samples[classes[0]] = list(set(mutations.index.droplevel(1)))
    samplesSet.update(list(set(mutations.index.droplevel(1))))
    if len(mutationsPath) > 1:
        for i in range(1, len(mutationsPath), 1):
            if classes[i] == "Class III":
                newMutations = pd.read_csv(
                    mutationsPath[i], sep="\t", names=headerFields, header=0
                )  # , skiprows=[0],)
            elif classes[i] == "Non-clust":
                newMutations = pd.read_csv(
                    mutationsPath[i],
                    sep="\t",
                    names=[
                        "project",
                        "samples",
                        "ID",
                        "genome",
                        "mutType",
                        "chr",
                        "start",
                        "end",
                        "ref",
                        "alt",
                        "mutClass",
                        "IMDplot",
                        "IMD",
                    ],
                    header=0,
                    skiprows=[0],
                    engine="python",
                )
            else:
                newMutations = pd.read_csv(
                    mutationsPath[i], sep="\t", names=headerFields[:-1], header=0
                )
            newMutations["chr"] = newMutations["chr"].astype(str)
            newMutations = newMutations.set_index(["samples", "chr"])
            newMutations = newMutations.sort_index()
            allMutations[classes[i]] = newMutations
            samples[classes[i]] = list(set(newMutations.index.droplevel(1)))
            samplesSet.update(list(set(newMutations.index.droplevel(1))))
            try:
                if not genome:
                    genome = newMutations.head(1).iloc[0, 2]
            except:
                continue

    samplesSet = set([str(x) for x in list(samplesSet) if x != "."])

    # Set up chromosome parameters:
    # chrom_path = "/anaconda3/lib/python3.6/site-packages/SigProfilerMatrixGenerator/references/chromosomes/tsb/" + genome + "/"
    chroms = [
        x.split(".")[0]
        for x in os.listdir(chrom_path)
        if x != ".DS_Store"
        and x[0] != "G"
        and x[0] != "M"
        and x != "BED_" + genome + "_proportions.txt"
        and "proportions" not in x
    ]
    chroms = sorted(
        chroms,
        key=lambda x: [
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "X",
            "Y",
        ].index(x[0:]),
    )
    pp = PdfPages(outputPath + "rainfallPlots_clustered_" + plot_suffix + ".pdf")
    colors = ["red", "orange", "black", "green", "blue", "grey"]
    colors = {
        "Class IA": "red",
        "Class IB": "black",
        "Class IC": "green",
        "Class II": "orange",
        "Class III": "blue",
        "Non-clust": "grey",
        "Simulation": "grey",
        "Clust": "orange",
    }

    subClassConversion = {
        "Class IA": "DBS",
        "Class IB": "MBS",
        "Class IC": "OMIKLI",
        "Class II": "KATAEGIS",
        "Class III": "OTHER",
        "Non-clust": "Non-clust",
        "Simulation": "Simulation",
        "Clust": "Clust",
    }

    if genome not in chromLengths:
        chromLengths[genome] = {}
        for chrom in chroms:
            with open(chrom_path + chrom + ".txt", "rb") as f:
                chromLengths[genome][chrom] = len(f.read())

    if correction:
        chromLengths2 = {}
        chroms = [
            x.split(".")[0]
            for x in os.listdir(chrom_path)
            if x != ".DS_Store"
            and x[0] != "G"
            and x[0] != "M"
            and x != "BED_" + genome + "_proportions.txt"
            and "proportions" not in x
        ]
        chroms = sorted(
            chroms,
            key=lambda x: [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "7",
                "8",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22",
                "X",
                "Y",
            ].index(x[0:]),
        )
        chromLengths2[genome] = {}
        totalLength = 0
        for i in range(0, len(chroms), 1):
            with open(chrom_path + chroms[i] + ".txt", "rb") as f:
                chromLengths2[genome][chroms[i]] = totalLength
                chrom_length = len(f.read())
                totalLength += chrom_length
        binsDensity = []
        for i in range(0, totalLength, windowSize):
            binsDensity.append(i)
        binsDensity.append(totalLength)

    falsePositives = {}
    falsePositives_1000 = {}

    tokenize = re.compile(r"(\d+)|(\D+)").findall

    def natural_sortkey(string):
        return tuple(int(num) if num else alpha for num, alpha in tokenize(string))

    try:
        sortedSampleSet = sorted(samplesSet, key=natural_sortkey)
    except:
        sortedSampleSet = samplesSet

    totalMutsIMD = {}
    totalMuts1000 = {}

    if not chrom_based:
        genomeLength = sum([x for x in chromLengths[genome].values()])
        xtick_pos = []
        totalPos = 0
        for chrom in chroms:
            newPos = chromLengths[genome][chrom]
            xtick_pos.append(newPos / 2 + totalPos)
            totalPos += newPos
        bins = []
        for i in range(0, genomeLength, windowSize):
            bins.append(i)
        bins.append(genomeLength)
        count = 1
        for sample in sortedSampleSet:
            chrom_startsAll = {}
            # if correction:
            # 	# regions = list(imdsDataSample_corrected[sample.split("_")[0]].keys())
            # 	regions = list(imdsDataSample_corrected[sample.split(".")[0]].keys())
            zorderPlot = 0
            # falsePositives[sample.split("_")[0]] = []
            # falsePositives_1000[sample.split("_")[0]] = []
            falsePositives[sample.split(".")[0]] = []
            falsePositives_1000[sample.split(".")[0]] = []
            falsePositivesChrom = []
            falsePositivesChrom_1000 = []
            totalMutsIMDChrom = []
            totalMuts1000Chrom = []
            densityMuts = []
            count += 1
            plot1 = plt.figure(figsize=(12, 8))
            plt.rc("axes", edgecolor="lightgray")
            panel1 = plt.axes([0.07, 0.09, 0.9, 0.6])
            panel2 = plt.axes([0.07, 0.75, 0.9, 0.15])
            if not aggregated and not chrom_based_IMD:
                try:
                    imd_line = imdsData[sample.split("_")[0]]
                except:
                    imd_line = imdsData[sample]
                panel1.axhline(imd_line, color="red")

            ymax = 0
            classColor = 0
            firstClass = True
            firstChrom = True
            chrom_prop_start = 0.075
            for subclass in classes:
                chrom_start = 0
                if sample not in samples[subclass]:
                    classColor += 1
                    continue
                for chrom in chroms:
                    chrom_startsAll[chrom] = chrom_start
                    if firstClass:
                        centromere_start = centromeres[genome][chrom][0]
                        centromere_end = centromeres[genome][chrom][1]
                        centromere_startProp = (
                            (chrom_start + centromere_start) / genomeLength
                        ) * 0.9 + 0.07
                        centromereLength = (
                            (centromere_end - centromere_start) / genomeLength * 0.9
                        )
                        chrom_prop_start = (chrom_start / genomeLength) * 0.9 + 0.07
                        chrom_prop = (
                            (chromLengths[genome][chrom] + 10000) / genomeLength
                        ) * 0.9 - 0.003
                        panel1.add_patch(
                            plt.Rectangle(
                                (chrom_prop_start, 0.06),
                                chrom_prop,
                                0.01,
                                facecolor="lightgrey",
                                clip_on=False,
                                transform=plt.gcf().transFigure,
                            )
                        )
                        panel1.add_patch(
                            plt.Rectangle(
                                (centromere_startProp, 0.06),
                                centromereLength,
                                0.01,
                                facecolor="red",
                                clip_on=False,
                                transform=plt.gcf().transFigure,
                            )
                        )
                        if chrom_based_IMD:
                            try:
                                imd_line = imdsData[sample.split(".")[0]][chrom]
                                # imd_line = imdsData[sample.split("_")[0]][chrom]
                            except:
                                imd_line = imdsDataSample[sample.split(".")[0]]
                                # imd_line = imdsDataSample[sample.split("_")[0]]
                            panel1.plot(
                                [
                                    chrom_start,
                                    chrom_start + chromLengths[genome][chrom],
                                ],
                                [imd_line, imd_line],
                                color="red",
                            )

                    if str(chrom) in [
                        str(x) for x in allMutations[subclass].loc[sample].index
                    ]:
                        try:
                            pos = allMutations[subclass].loc[sample, str(chrom)]
                        except:
                            pos = allMutations[subclass].loc[sample, int(chrom)]

                        # try:
                        if True:
                            starts = pos.loc[(sample, chrom), "start"]
                            IMDsPlotRecorded = pos.loc[(sample, chrom), "IMDplot"]
                            imds_recorded = pos.loc[(sample, chrom), "IMD"]
                            if type(starts) == np.int64:
                                starts = [starts]
                                imds_recorded = [imds_recorded]
                                IMDsPlotRecorded = [IMDsPlotRecorded]
                            else:
                                starts = list(starts)
                                imds_recorded = list(imds_recorded)
                                IMDsPlotRecorded = list(IMDsPlotRecorded)
                            if contexts != "ID":
                                newGroup = [0]
                                if (
                                    subclass != "Non-clust"
                                ):  # and subclass != "Class III":
                                    newGroup = pos.loc[(sample, chrom), "group"]
                                    if type(newGroup) == np.int64:
                                        newGroup = np.array([newGroup])
                                    else:
                                        newGroup = np.array(newGroup)
                                    if newGroup.size > 0:
                                        newGroup = [
                                            x + 1
                                            for x in list(
                                                np.where(newGroup[:-1] != newGroup[1:])[
                                                    0
                                                ]
                                            )
                                        ]
                                    newGroup = [0] + newGroup

                            imds = [y - x for x, y in zip(starts, starts[1:])]
                            minIMDs = [
                                min(x, y)
                                for x, y in zip(
                                    [float("inf")] + imds, imds + [float("inf")]
                                )
                            ]
                            if not aggregated:
                                # if subclass == "Class III":
                                # 	if chrom_based_IMD:
                                # 		try:
                                # 			plotIMDs = [[x,y] for x,y in zip(starts, minIMDs) if y <= imdsData[sample.split(".")[0]][chrom] or y < 10000]
                                # 			# plotIMDs = [[x,y] for x,y in zip(starts, minIMDs) if y <= imdsData[sample.split("_")[0]][chrom] or y < 10000]
                                # 		except:
                                # 			plotIMDs = [[x,y] for x,y in zip(starts, minIMDs) if y <= imdsDataSample[sample.split(".")[0]] or y < 10000]
                                # 			# plotIMDs = [[x,y] for x,y in zip(starts, minIMDs) if y <= imdsDataSample[sample.split("_")[0]] or y < 10000]
                                # 	else:
                                # 		plotIMDs = [[x,y] for x,y in zip(starts, imds_recorded)]# if y <= imdsData[sample.split("_")[0]] or y < 10000]
                                # 		# plotIMDs = [[x,y] for x,y in zip(starts, minIMDs) if y <= imdsData[sample.split("_")[0]] or y < 10000]
                                # else:
                                # 	# plotIMDs = [[x,y] for x,y in zip(starts, imds) if y <= imdsDataSample[sample.split(".")[0]] or y < 10000]
                                # 	plotIMDs = [[x,y] for x,y in zip(starts, minIMDs)]
                                if contexts != "ID":
                                    if subclass == "Class III":
                                        plotIMDs = [
                                            [x, int(y)]
                                            for z, (x, y) in enumerate(
                                                zip(starts, imds_recorded)
                                            )
                                            if (y != "c" and y != "d")
                                            and (
                                                z not in newGroup
                                                or len(starts) == 1
                                                or len(starts) == len(newGroup)
                                            )
                                        ]
                                    else:
                                        # plotIMDs = [[x,int(y)] for x,y in zip(starts, IMDsPlotRecorded) if y != 'c' ]
                                        plotIMDs = [
                                            [x, int(y)]
                                            for z, (x, y) in enumerate(
                                                zip(starts, IMDsPlotRecorded)
                                            )
                                            if (y != "c" and y != "d")
                                            and (
                                                z not in newGroup
                                                or len(starts) == 1
                                                or len(starts) == len(newGroup)
                                            )
                                        ]
                                else:
                                    plotIMDs = [
                                        [x, y] for x, y in zip(starts, imds_recorded)
                                    ]  # if y <= imdsData[sample.split("_")[0]] or y < 10000]
                                if len(plotIMDs) == 0:
                                    continue
                                plotX, plotY = zip(*plotIMDs)
                                newMax = max(plotY)
                                # if chrom_based_IMD:
                                # 	# falsePositivesChrom.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <=imdsData[sample.split("_")[0]][chrom]]))
                                # 	falsePositivesChrom.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <=imdsData[sample.split(".")[0]][chrom]]))
                                # 	falsePositivesChrom_1000.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <=1000]))
                                # else:
                                # 	falsePositivesChrom.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <=imdsData[sample.split(".")[0]]]))
                                # 	# falsePositivesChrom.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <=imdsData[sample.split("_")[0]]]))
                                # 	falsePositivesChrom_1000.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <=1000]))

                                # try:
                                # 	totalMutsIMDChrom.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <= imdsData[sample.split("_")[0]][chrom]]))
                                # except:
                                # 	totalMutsIMDChrom.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <= imdsDataSample[sample.split("_")[0]]]))
                                # totalMuts1000Chrom.append(len([[x,y] for x,y in zip(starts, minIMDs) if y <= 1000]))
                                # if correction and subclass != "Non-clust":
                                # 	correctionMuts = []
                                # 	for imd, start in zip(minIMDs, starts):
                                # 	# 	position = int(start) + chromLengths2[genome][chrom]
                                # 	# 	try:
                                # 	# 		bisectRegion = regions[bisect.bisect_left(regions, position)]
                                # 	# 	except:
                                # 	# 		continue
                                # 	# 	if bisectRegion - position < windowSize:
                                # 	# 		imdCorrected = imdsDataSample_corrected[sample.split("_")[0]][bisectRegion]
                                # 		# if imd < imdCorrected and imd > imdsDataSample[sample.split("_")[0]]:
                                # 		correctionMuts.append(1)
                                # 	totalMutsIMDChrom.append(len(correctionMuts))
                            else:
                                plotIMDs = [
                                    [x, y]
                                    for x, y in zip(starts, minIMDs)
                                    if y <= 10000
                                ]
                                plotX, plotY = zip(*plotIMDs)
                                newMax = max(plotY)
                        # except:
                        # 	plotY = 0
                        # 	try:
                        # 		plotX = [int(pos.loc[(sample, chrom), 'start'])]
                        # 	except:
                        # 		plotY = [0 for x in pos.loc[(sample, chrom), 'start']]
                        # 		plotX = [int(x) for x in pos.loc[(sample, chrom), 'start']]
                        # 	newMax = 0

                        if newMax > ymax:
                            ymax = newMax
                        plotX = [x + chrom_start for x in plotX]
                        densityMuts += plotX
                        if subclass == "Non-clust":
                            zorderPlot = -10
                        panel1.scatter(
                            plotX,
                            plotY,
                            color=colors[subclass],
                            clip_on=False,
                            s=12,
                            label=subClassConversion[subclass],
                            zorder=zorderPlot,
                        )

                    chrom_start += chromLengths[genome][chrom] + 10000
                classColor += 1
                firstClass = False
                zorderPlot -= 1
            totalMutsIMD[sample] = sum(totalMutsIMDChrom)
            totalMuts1000[sample] = sum(totalMuts1000Chrom)

            # falsePositives[sample.split("_")[0]].append(sum(falsePositivesChrom))
            # falsePositives_1000[sample.split("_")[0]].append(sum(falsePositivesChrom_1000))
            falsePositives[sample.split(".")[0]].append(sum(falsePositivesChrom))
            falsePositives_1000[sample.split(".")[0]].append(
                sum(falsePositivesChrom_1000)
            )
            hist, bins = np.histogram(densityMuts, bins=bins, density=True)
            panel2.fill_between(bins[:-1], hist, zorder=10)
            if ymax > 10**5:
                panel1.set_ylim([1, ymax])
            else:
                panel1.set_ylim([1, 10**5])

            ############################################
            ############################################
            if bedRanges:
                if sample in bedFile.index:
                    currentRanges = bedFile.loc[sample]
                    prev_start = 0
                    prev_end = 0
                    prevIndex = 0
                    allRangesCurrent = currentRanges.to_numpy().tolist()
                    if not any(isinstance(i, list) for i in allRangesCurrent):
                        allRangesCurrent = [allRangesCurrent]
                    for ranges in allRangesCurrent:
                        start = ranges[1]
                        end = ranges[2]
                        currentChrom = str(ranges[0])
                        ampIndex = ranges[5]
                        if currentChrom in chrom_startsAll:
                            # print(currentChrom, start, ranges[-1])
                            panel1.axvline(
                                x=start + chrom_startsAll[currentChrom],
                                color=ranges[-1],
                                lw=1,
                                zorder=-20,
                                alpha=0.6,
                            )
                            panel1.axvline(
                                x=end + chrom_startsAll[currentChrom],
                                color=ranges[-1],
                                lw=1,
                                zorder=-20,
                                alpha=0.6,
                            )
                            # kw = dict(linestyle=None, lw=1,color=ranges[-1],connectionstyle="arc3,rad=0.2",clip_on=False)
                            # arc = mplpatches.FancyArrowPatch((start+chrom_startsAll[currentChrom], ymax), (end+chrom_startsAll[currentChrom],ymax) , **kw)
                            # panel1.add_patch(arc)
                            # if prevIndex != 0:
                            # 	if ampIndex == prevIndex:
                            # 		if start+chrom_startsAll[currentChrom] > prev_end:
                            # 			kw = dict(linestyle=None, lw=1,color=ranges[-1],connectionstyle="arc3,rad=0.2", clip_on=False)
                            # 			arc = mplpatches.FancyArrowPatch((start+chrom_startsAll[currentChrom], ymax), (prev_end,ymax) , **kw)
                            # 			panel1.add_patch(arc)
                            # 		elif start+chrom_startsAll[currentChrom] < prev_start:
                            # 			kw = dict(linestyle=None, lw=1,color=ranges[-1],connectionstyle="arc3,rad=0.2", clip_on=False)
                            # 			arc = mplpatches.FancyArrowPatch((start+chrom_startsAll[currentChrom], ymax), (prev_start,ymax) , **kw)
                            # 			panel1.add_patch(arc)
                            prev_start = start + chrom_startsAll[currentChrom]
                            prev_end = end + chrom_startsAll[currentChrom]
                            prevIndex = ampIndex

            ############################################
            ############################################

            panel1.set_yscale("log")
            panel1.set_title("Clustered mutations - " + sample, pad=150, fontsize=20)
            panel1.set_ylabel(
                "Distance between mutations in a single event (log10)", fontsize=12
            )
            panel1.spines["right"].set_visible(False)
            panel1.spines["top"].set_visible(False)
            panel1.spines["bottom"].set_visible(False)
            panel1.set_xticks(xtick_pos)
            panel1.set_xticklabels(["chr" + x for x in chroms], rotation=45)
            panel1.set_xlim([10000, genomeLength])
            panel2.set_xlim([10000, genomeLength])
            panel2.set_ylabel("Density", labelpad=10)
            panel1.tick_params(axis="y", which="major", labelsize=10)
            panel1.tick_params(axis="x", which="major", labelsize=10, pad=20, length=0)
            panel2.set_xticks([])
            panel2.set_xticklabels([])
            panel2.set_ylim([0, max(hist) * 1.2])
            panel2.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
            handles, labels = panel1.get_legend_handles_labels()
            goodLabels = [labels.index(x) for x in sorted(list(set(labels)))]
            panel1.legend(
                [handles[x] for x in goodLabels],
                [labels[x] for x in goodLabels],
                loc="best",
                prop={"size": 10},
            )
            pp.savefig(plot1)
            plt.close()

        # with open(projectPath_parent + "output/simulations/data/originalCounts" + path_suffix2 + ".pickle", "wb") as f:
        # 	pickle.dump(totalMutsIMD, f)
        with open(
            projectPath_parent
            + "output/simulations/data/originalCounts1000"
            + path_suffix2
            + ".pickle",
            "wb",
        ) as f:
            pickle.dump(totalMuts1000, f)
    else:
        for chrom in chroms:
            genomeLength = chromLengths[genome][chrom]
            xtick_pos = []
            bins = []
            for i in range(0, genomeLength, windowSize):
                bins.append(i)
            bins.append(genomeLength)
            count = 1
            for sample in samplesSet:
                densityMuts = []
                count += 1
                plot1 = plt.figure(figsize=(12, 8))
                plt.rc("axes", edgecolor="lightgray")
                panel1 = plt.axes([0.07, 0.09, 0.9, 0.6])
                panel2 = plt.axes([0.07, 0.75, 0.9, 0.15])
                ymax = 0
                classColor = 0
                firstClass = True
                firstChrom = True
                chrom_prop_start = 0.075
                for subclass in classes:
                    chrom_start = 0
                    if sample not in samples[subclass]:
                        classColor += 1
                        continue
                    if firstClass:
                        centromere_start = centromeres[genome][chrom][0]
                        centromere_end = centromeres[genome][chrom][1]
                        centromere_startProp = (
                            centromere_start / genomeLength
                        ) * 0.9 + 0.07
                        centromereLength = (
                            (centromere_end - centromere_start) / genomeLength * 0.9
                        )
                        chrom_prop_start = (chrom_start / genomeLength) * 0.9 + 0.07
                        chrom_prop = (
                            (chromLengths[genome][chrom] + 10000) / genomeLength
                        ) * 0.9 - 0.003
                        panel1.add_patch(
                            plt.Rectangle(
                                (chrom_prop_start, 0.06),
                                chrom_prop,
                                0.01,
                                facecolor="lightgrey",
                                clip_on=False,
                                transform=plt.gcf().transFigure,
                            )
                        )
                        panel1.add_patch(
                            plt.Rectangle(
                                (centromere_startProp, 0.06),
                                centromereLength,
                                0.01,
                                facecolor="red",
                                clip_on=False,
                                transform=plt.gcf().transFigure,
                            )
                        )

                    if chrom in allMutations[subclass].loc[sample].index:
                        pos = allMutations[subclass].loc[sample, chrom]
                        try:
                            starts = list(pos.loc[(sample, chrom), "start"])
                            imds = [y - x for x, y in zip(starts, starts[1:])]
                            minIMDs = [
                                min(x, y)
                                for x, y in zip(
                                    [float("inf")] + imds, imds + [float("inf")]
                                )
                            ]
                            if not aggregated:
                                if subclass == "Class III":
                                    plotIMDs = [
                                        [x, y]
                                        for x, y in zip(starts, minIMDs)
                                        if y <= imdsData[sample]
                                    ]
                                else:
                                    plotIMDs = [[x, y] for x, y in zip(starts, minIMDs)]
                                plotX, plotY = zip(*plotIMDs)
                                newMax = max(plotY)
                            else:
                                plotIMDs = [
                                    [x, y]
                                    for x, y in zip(starts, minIMDs)
                                    if y <= 10000
                                ]
                                plotX, plotY = zip(*plotIMDs)
                                newMax = max(plotY)

                        except:
                            plotY = 0
                            try:
                                plotX = [int(pos.loc[(sample, chrom), "start"])]
                            except:
                                plotY = [0 for x in pos.loc[(sample, chrom), "start"]]
                                plotX = [
                                    int(x) for x in pos.loc[(sample, chrom), "start"]
                                ]
                            newMax = 0

                        if newMax > ymax and newMax < float("inf"):
                            ymax = newMax
                        plotX = [x + chrom_start for x in plotX]
                        densityMuts += plotX
                        panel1.scatter(
                            plotX,
                            plotY,
                            color=colors[subclass],
                            clip_on=False,
                            s=12,
                            label=subclass,
                        )
                    chrom_start += chromLengths[genome][chrom] + 10000
                classColor += 1
                firstClass = False
                hist, bins = np.histogram(densityMuts, bins=bins, density=True)
                panel2.fill_between(bins[:-1], hist, zorder=10)
                if ymax > 10**5:
                    panel1.set_ylim([1, ymax])
                else:
                    panel1.set_ylim([1, 10**5])
                panel1.set_yscale("log")
                panel1.set_title(
                    "Clustered mutations chromosome " + chrom + " - " + sample,
                    pad=150,
                    fontsize=20,
                )
                panel1.set_ylabel(
                    "Distance between mutations in a single event (log10)", fontsize=12
                )
                panel1.spines["right"].set_visible(False)
                panel1.spines["top"].set_visible(False)
                panel1.spines["bottom"].set_visible(False)
                panel1.set_xticks(xtick_pos)
                panel1.set_xticklabels(["chr" + x for x in chroms], rotation=45)
                panel1.set_xlim([10000, genomeLength])
                panel2.set_xlim([10000, genomeLength])
                panel2.set_ylabel("Density", labelpad=10)
                panel1.tick_params(axis="y", which="major", labelsize=10)
                panel1.tick_params(
                    axis="x", which="major", labelsize=10, pad=20, length=0
                )
                panel2.set_xticks([])
                panel2.set_xticklabels([])
                try:
                    panel2.set_ylim([0, max(hist) * 1.2])
                except:
                    pass
                panel2.grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
                handles, labels = panel1.get_legend_handles_labels()
                goodLabels = [labels.index(x) for x in sorted(list(set(labels)))]
                panel1.legend(
                    [handles[x] for x in goodLabels],
                    [labels[x] for x in goodLabels],
                    loc="best",
                    prop={"size": 10},
                )
                pp.savefig(plot1)
                plt.close()
    pp.close()
