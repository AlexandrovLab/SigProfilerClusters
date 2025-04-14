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
import pickle
import bisect
from scipy.signal import find_peaks
import random
import multiprocessing as mp

# import seaborn as sns
# from collections import Iterable
from itertools import chain
import sys


warnings.filterwarnings("ignore", message="invalid value encountered in long_scalars")
warnings.filterwarnings(
    "ignore", message="Data has no positive values, and therefore cannot be log-scaled."
)


def z_test(x, mu, sigma):
    """
    Performs a z-test for statistical comparisons of simulated and original data.

    Parameters:
                    x	->	observed number in original sample (mutation count; int)
               mu	->	average number observed in simulations (average mutation count; float)
            sigma	->	standard deviation in simulations (float)

    Returns:
            z	->	z-score (float)
            p	->	associated p_value (float)
    """
    z = (x - mu) / sigma
    p = 2 * (1 - norm.cdf(z))
    return (z, p)


def plot_clustered(orig_bins, sim_bins, bincenters2, panel4, lower_CI, upper_CI):
    """
    Plots the histogram of IMDs in the middle right-hand column of the intradistance plots for
    clustered mutations.

    Parameters:
              orig_bins	->	histogram heights for clustered mutations in original sample (numpy array)
               sim_bins	->	histogram heights for clustered mutations in simulated sample (numpy array)
            bincenters2	->	center bounds of each histogram bin  (numpy array)
                     panel4	->	axis location for the given clustered plot (matplotlib panel/axes object)
               lower_CI	->	lower 95% confidence interval values for each bin (list)
               upper_CI	->	upper 95% confidence interval values for each bin (list)

    Returns:
            None

    Outputs:
            panel4	->	clustered histogram of IMDs onto a given figure
    """
    (sim,) = panel4.plot(
        bincenters2[: len(sim_bins)],
        sim_bins,
        "-",
        marker="o",
        markersize=2,
        color="red",
    )
    panel4.fill_between(
        bincenters2[: len(sim_bins)],
        upper_CI[: len(sim_bins)],
        lower_CI[: len(sim_bins)],
        alpha=0.5,
        color="red",
        zorder=1000,
    )
    (orig,) = panel4.plot(
        bincenters2[: len(orig_bins)],
        orig_bins,
        "-",
        marker="o",
        markersize=2,
        color="green",
    )

    panel4.set_yscale("log")
    panel4.set_xscale("log")


def plot_non_clustered(orig_bins, sim_bins, bincenters2, panel6, lower_CI, upper_CI):
    """
    Plots the histogram of IMDs in the lower right-hand column of the intradistance plots for
    nonclustered mutations.

    Parameters:
              orig_bins	->	histogram heights for nonclustered mutations in original sample (numpy array)
               sim_bins	->	histogram heights for nonclustered mutations in simulated sample (numpy array)
            bincenters2	->	center bounds of each histogram bin  (numpy array)
                     panel6	->	axis location for the given nonclustered plot (matplotlib panel/axes object)
               lower_CI	->	lower 95% confidence interval values for each bin (list)
               upper_CI	->	upper 95% confidence interval values for each bin (list)

    Returns:
            None

    Outputs:
            panel6	->	nonclustered histogram of IMDs onto a given figure
    """
    (sim,) = panel6.plot(
        bincenters2[-len(sim_bins) :],
        sim_bins,
        "-",
        marker="o",
        markersize=2,
        color="red",
    )
    panel6.fill_between(
        bincenters2[-len(sim_bins) :],
        upper_CI[-len(sim_bins) :],
        lower_CI[-len(sim_bins) :],
        alpha=0.5,
        color="red",
        zorder=1000,
    )
    (orig,) = panel6.plot(
        bincenters2[-len(orig_bins) :],
        orig_bins,
        "-",
        marker="o",
        markersize=2,
        color="green",
    )

    panel6.set_yscale("log")
    panel6.set_xscale("log")
    panel6.set_xlabel("Inter-Mutational Distance (IMD)", fontweight="bold", fontsize=12)


def plot_hist(
    y2,
    bincenters2,
    q_vals,
    interval_line,
    orig_mutations,
    avgSimCounts,
    stdSimCounts,
    imd,
    lower_CI,
    upper_CI,
    lower_CI_refined,
    upper_CI_refined,
    avg_bin_counts,
    sample,
    original,
    panel2,
    panel3,
    panel4,
    panel5,
    panel6,
):
    """
    Plots the histogram of IMDs in the upper right-hand column of the intradistance plots for
    all mutations and calls the functions to generate the plots for solely clustered and nonclustered mutations.

    Parameters:
                                            y2	->	histogram heights for all mutations in original sample (numpy array)
                       bincenters2	->	center bounds of each histogram bin; uses a 2^n binning  (numpy array)
                                    q_vals	->	corrected p_value of signifance using FDR correction for a given sample at the given IMD cutoff (float)
                     interval_line	->	the index for the IMD value used for the clustering cutoff; must be a valid index 0<interval_line<len(bincenters2) (integer)
                    orig_mutations	->	the number of clustered mutations observed in the original sample using the refinded IMD (int)
                      avgSimCounts	->	the average number of clustered mutations observed in the simulated sample using the refinded IMD (float)
                      stdSimCounts	->	the standard deviation of clustered mutations observed in the simulated sample using the refinded IMD (int)
                                       imd	->	the refined IMD (int)
                              lower_CI	->	lower 95% confidence interval values for each bin (list)
                              upper_CI	->	upper 95% confidence interval values for each bin (list)
              lower_CI_refined	->	lower 95% confidence interval value (int)
              upper_CI_refined	->	upper 95% confidence interval value (int)
                    avg_bin_counts	-> 	histogram heights for all mutations in simulated sample (numpy array)
                                    sample	->	the current sample that is being plotted
                              original	->  option to plot original sample or not (boolean; default=True)
                                    panel2	->	axis location for the given histogram plot of all IMDs (matplotlib panel/axes object)
                                    panel3	->	axis location for the given clustered SBS96 plot (matplotlib panel/axes object)
                                    panel4	->	axis location for the given clustered histogram plot (matplotlib panel/axes object)
                                    panel5	->	axis location for the given nonclustered SBS96 plot (matplotlib panel/axes object)
                                    panel6	->	axis location for the given nonclustered histogram plot (matplotlib panel/axes object)

    Returns:
            True	->	boolean to mark the plotting completion

    Outputs:
            Histogram plots on a newly generated figure. The figure is completed in subsequent steps in the main function which will plot the SBS96 plots.
    """
    axes = plt.gca()
    distance_cutoff = bincenters2[interval_line]
    if original:
        panel2.axvline(x=imd, linestyle="--", linewidth=0.5, color="darkred")
        extra = Rectangle(
            (0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0
        )
        panel3.text(
            0.08, 0.87, "Total=" + str(sum(y2)), transform=plt.gcf().transFigure
        )
        panel3.text(
            0.08,
            0.575,
            "Clustered=" + str(orig_mutations) + "(IMD<" + str(imd) + ")",
            transform=plt.gcf().transFigure,
        )
        panel3.text(
            0.08,
            0.56,
            "Expected="
            + str(round(avgSimCounts, 0))
            + "("
            + str(lower_CI_refined)
            + "-"
            + str(upper_CI_refined)
            + ")",
            transform=plt.gcf().transFigure,
        )
        # panel3.text(0.08, 0.56, "Expected=" + str(statistics.mean([sum(lower_CI[:interval_line]), sum(upper_CI[:interval_line])])) + "(" +str(sum(lower_CI[:interval_line])) + "-" + str(sum(upper_CI[:interval_line])) + ")",transform=plt.gcf().transFigure)
        panel5.text(
            0.08,
            0.28,
            "Non-Clustered=" + str(sum(y2) - orig_mutations),
            transform=plt.gcf().transFigure,
        )
    panel2.set_yscale("log")
    panel2.set_xscale("log")
    l = 0

    y2 = y2[l:]
    avg_bin_counts = avg_bin_counts[l:]
    bincenters2 = bincenters2[l:]
    upper_CI = upper_CI[l:]
    lower_CI = lower_CI[l:]

    (sim,) = panel2.plot(
        bincenters2, avg_bin_counts, "-", marker="o", markersize=2, color="red"
    )
    panel2.fill_between(
        bincenters2, upper_CI, lower_CI, alpha=0.5, color="red", zorder=1000
    )
    if original:
        (orig,) = panel2.plot(
            bincenters2, y2, "-", marker="o", markersize=2, color="green"
        )

    if original:
        panel2.legend(
            [sim, orig, extra],
            [
                "simulated",
                "real samples",
                "q_value = " + "{:.2E}".format(Decimal(q_vals)),
            ],
        )
    else:
        panel2.legend([sim], ["simulated"])

    plot_clustered(
        y2[: interval_line - l],
        avg_bin_counts[: interval_line - l],
        bincenters2,
        panel4,
        lower_CI,
        upper_CI,
    )
    plot_non_clustered(
        y2[interval_line - l :],
        avg_bin_counts[interval_line - l :],
        bincenters2,
        panel6,
        lower_CI,
        upper_CI,
    )

    return True


def refineIMD(
    distances,
    distances_orig,
    origCounts,
    avgSimCounts,
    interval_line,
    left_bound,
    right_bound,
    CI,
    lower_CI,
    upper_CI,
    sigValue=0.01,
):
    """
    Refines the discretized IMD cutoff using a binary search algorithm

    Parameters:
             distances	->	the IMD distances across all simulations (list; nested list of lists: a list of 100 lists with one nested list per iteration of simulations)
    distances_orig	->	the IMD distances across the original sample (list)
            origCounts	->	the counts in each bin from the original sample (numpy array)
      avgSimCounts	->	the counts in each bin from the simulated sample (numpy array)
     interval_line	->	the index for the IMD value used for the clustering cutoff; must be a valid index 0<interval_line<len(bincenters2) (integer)
            left_bound	->	the current IMD value; the last significant bin (int)
       right_bound	->	the IMD value immediately following the last subsequent bin (int; i.e. left:256, right:512)
                            CI	->	95% confidence interval range based upon the number of simulations provided (integer)
              lower_CI	->	lower 95% confidence interval values for each bin (list)
              upper_CI	->	upper 95% confidence interval values for each bin (list)
              sigValue	->	significant value threshold (float; default=0.01)


    Returns:
                     middle	->	the refined IMD (int)
              q_vals[0]	->	the refined q_value using the new IMD (float)
              avgSim[0]	->	the average number of clustered mutations in the simulated data (float)
            stdevSim[0]	->	the standard deviation of clustered mutations in the simulated data (float)
                    upperCI	->	the upper 95% value of clustered mutations in the simulated data (int)
                    lowerCI	->	the lower 95% value of clustered mutations in the simulated data (int)
    """
    left = left_bound
    right = right_bound
    middle = 0
    cumulativeSumOrig = sum(origCounts[:interval_line])
    cumulativeSumSim = sum(avgSimCounts[:interval_line])
    cumulativeSumUpperCI = sum(upper_CI[:interval_line])
    cumulativeSumLowerCI = sum(lower_CI[:interval_line])
    while left <= right:
        middle = (left + right) // 2
        bins = [left_bound, middle, right]
        ys = np.histogram(distances_orig, bins=bins)[0]
        ys[1] += ys[0] + cumulativeSumOrig
        ys[0] += cumulativeSumOrig
        total_bin_counts = []
        for dist in distances:
            ysSims = np.histogram(dist, bins=bins)[0]
            total_bin_counts.append(
                [ysSims[0] + cumulativeSumSim, ysSims[0] + cumulativeSumSim + ysSims[1]]
            )
        avgSim = np.mean(total_bin_counts, axis=0)
        stdevSim = np.std(total_bin_counts, axis=0)
        total_bin_counts.sort()
        upperCI = total_bin_counts[CI - 1][0]
        lowerCI = total_bin_counts[-CI][0]
        z, p = z_test(ys, avgSim, stdevSim)
        q_vals = sm.fdrcorrection(p)[1]

        if middle == right or middle == left:
            return (middle, q_vals[0], avgSim[0], stdevSim[0], upperCI, lowerCI)
        if q_vals[0] < sigValue and ys[0] / (avgSim[0] + ys[0]) > 0.9:
            left = middle
        else:
            right = middle


def first_run(
    distances,
    distances_orig_all,
    distances_orig,
    vcf_path_clust,
    vcf_path_nonClust,
    sample,
    original,
    sim_count,
    project,
    genome,
    clustering_vaf,
    correctionData=False,
    correction=False,
    regions=None,
    imds_corrected=None,
    windowSize=None,
    chromLengths=None,
):
    """
    Calculate the IMD cutoffs on a per sample basis. Newest addition refines this IMD cutoff using an adpated binary search algorithm rather than using a discretized binning method.

    Parameters:
                             distances	->	the IMD distances across all simulations (list; nested list of lists: a list of 100 lists with one nested list per iteration of simulations)
            distances_orig_all	-> 	all of the IMD distances along with the paired mutation information for each distance (list; suggest adding a print statement: print(distances_orig_all[0]) to see the format)
                    distances_orig	-> 	the IMD distances across the original sample (list)
                    vcf_path_clust	->	path to save the clustered mutations (string)
             vcf_path_nonClust	->	path to save the nonclustered mutations (string)
                                    sample	->	the current sample that is being plotted
                              original	->	option to plot original sample or not (boolean; default=True)
                             sim_count	->	the total number of simulations provided
                               project	->	user provided project name (string)
                                    genome	->	the reference genome used for the given analysis (string)
                    clustering_vaf	->	optional paramater to save the VAF from the original mutation files. Required if performing subclassification (default=False; switched to True if subclassify=True)
                    correctionData	->	optional parameter that informs the function whether to save the clustered mutations or to simply calculate a region-specific IMD (boolean; default=False)
                            correction	->	optional parameter to perform a genome-wide mutational density correction (boolean; default=False)
                               regions	->	a list of regions that contain higher than 1.25 mutation density compared with simulations using a 10MB window. The regions are individual integers recorded as a cumulative sum across the whole genome (list; default=None; ie chr2:10000 -> 10000+len(chrom1))
                    imds_corrected	->	a dictionary of corrected IMDs per sample per region (dictionary; default=None)
                            windowSize	->	the window size used to calculate the mutation densities across the genome (integer; default=None)
                      chromLengths	->	a dictionary of the cumulative chromosome lengths for a given reference genome (dictionary)

    Returns:
                                            y2	->	the counts per bin of the original sample's mutations (numpy array)
                       bincenters2	->  the center bounds for each bin (numpy array)
                                     q_val	->	the correct p_value using an FDR correction at a given IMD cutoff (float)
                     interval_line	->	index of the IMD cutoff used to index the bincenters and y2 data structures (int)
       len(clustered_muts)	->	the number of clustered mutations filtered out for a given sample (int)
                                    avgSim	->	the number of clustered mutations filtered out of the simulations for a given sample (float)
                              stdevSim	->	the standard deviation of mutations filtered out of the simulations for a given sample (float)
                      distance_cut	->	the refined IMD cutoff (int)
                              lower_CI	->	lower 95% confidence interval values for each bin (list)
                              upper_CI	->	upper 95% confidence interval values for each bin (list)
              lower_CI_refined	->	lower 95% confidence interval value (int)
              upper_CI_refined	->	upper 95% confidence interval value (int)
                    avg_bin_counts	->  the bin counts for the average number of mutations in each bin across the simulations

    Outputs:
            A clustered and non-clustered mutation file
    """
    # combinedSimDensity = [ele for ele in densityMutsSim.values()]
    # combinedSimDensity = list(chain.from_iterable(distances))
    # distances_orig_all_only = [int(x[0]) for x in distances_orig_all]
    # # print(distances_orig_all)
    # x,y = sns.kdeplot(np.log10(combinedSimDensity), bw_adjust=20, gridsize=200).get_lines()[0].get_data()
    # plt.close()
    # distances_orig_all_only = np.log10(distances_orig_all_only)
    # imdsProbabilities = []
    # for dist in distances_orig_all_only:
    # 	# print(dist, y[bisect.bisect_left(x, dist)])
    # 	imdsProbabilities.append(y[bisect.bisect_left(x, dist)])

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
        total_bin_counts.append(np.histogram(dist, bins=bin1)[0])

    avg_bin_counts = []
    upper_CI = []
    lower_CI = []
    mean_counts = []
    std_dev = []

    cumulativeBackground = [np.cumsum(x) for x in total_bin_counts]
    mean_counts_cum = []
    std_dev_cum = []

    CI = int(round(0.95 * sim_count, 0))

    for i in range(0, len(total_bin_counts[0]), 1):
        count = 0
        bin_std_dev = []
        bin_std_dev_cum = []
        for l in range(0, len(total_bin_counts), 1):
            count += total_bin_counts[l][i]
            bin_std_dev.append(total_bin_counts[l][i])
            bin_std_dev_cum.append(cumulativeBackground[l][i])

        bin_std_dev.sort()
        std_dev.append(stdev([int(x) for x in bin_std_dev]))
        mean_counts.append(mean(bin_std_dev))
        avg_bin_counts.append(int(round(median(bin_std_dev), 0)))
        upper_CI.append(bin_std_dev[CI - 1])
        lower_CI.append(bin_std_dev[-CI])

        bin_std_dev_cum.sort()
        std_dev_cum.append(stdev([int(x) for x in bin_std_dev_cum]))
        mean_counts_cum.append(median(bin_std_dev_cum))

    y2, binEdges2 = np.histogram(distances_orig, bins=bin1)

    bincenters2 = binEdges2[1:-1]
    y2 = y2[1:]
    avg_bin_counts = avg_bin_counts[1:]
    mean_counts = mean_counts[1:]
    std_dev = std_dev[1:]
    upper_CI = upper_CI[1:]
    lower_CI = lower_CI[1:]

    total_mutations = 0
    sim_mutations = 0
    orig_mutations = 0

    interval_line = 0
    previous_interval = 0
    interval_percent = 0

    p_vals = []
    percent_clustered = 0.90
    if correctionData:
        percent_clustered = 0.95

    current_orig_mutations_cum = 0
    current_mean_cum = 0
    current_stdev_cum = 0

    for i in range(0, len(avg_bin_counts), 1):
        current_orig_mutations = y2[i]
        current_mean = mean_counts[i]
        current_stdev = std_dev[i]

        current_orig_mutations_cum += y2[i]
        current_mean_cum += mean_counts[i]
        current_stdev_cum += std_dev[i]

        if current_stdev == 0:
            p = 0
        else:
            # z, p = z_test (current_orig_mutations, current_mean, current_stdev)
            z, p = z_test(
                current_orig_mutations_cum, current_mean_cum, current_stdev_cum
            )
        p_vals.append(p)
    q_vals = sm.fdrcorrection(p_vals)[1]

    for i in range(0, len(avg_bin_counts), 1):
        sim_mutations += avg_bin_counts[i]
        orig_mutations += y2[i]
        total_mutations = total_mutations + avg_bin_counts[i] + y2[i]

        current_orig_mutations = y2[i]
        current_mean = mean_counts[i]
        current_stdev = std_dev[i]
        current_q = q_vals[i]
        if orig_mutations / total_mutations < percent_clustered or current_q > 0.01:
            # if correctionData:
            # 	print(orig_mutations/total_mutations, previous_interval)
            if abs((orig_mutations / total_mutations) - percent_clustered) < abs(
                previous_interval - percent_clustered
            ):
                if current_q > 0.01:
                    interval_line = i - 1
                    orig_mutations -= y2[i]
                    interval_percent = previous_interval * 100
                else:
                    interval_line = i
                    interval_percent = orig_mutations / total_mutations * 100
            else:
                interval_line = i - 1
                orig_mutations -= y2[i]
                interval_percent = previous_interval * 100
            break
        previous_interval = orig_mutations / total_mutations
    interval_line = abs(interval_line)
    distance_cut = bincenters2[interval_line]
    sigValue = 0.01
    distance_cut, q_val, avgSim, stdevSim, upper_CI_refined, lower_CI_refined = (
        refineIMD(
            distances,
            distances_orig,
            y2,
            avg_bin_counts,
            interval_line,
            bincenters2[interval_line],
            bincenters2[interval_line + 1],
            CI,
            lower_CI,
            upper_CI,
            sigValue,
        )
    )
    if distance_cut > 10000:
        distance_cut = 10000

    if not correctionData:
        if correction:
            if len(regions) > 0:
                clustered_muts = [
                    x[1:] + [x[0]]
                    for x in distances_orig_all
                    if int(x[0]) <= distance_cut
                    or (
                        regions[catch(x, regions, chromLengths, genome)]
                        - (int(x[3]) + chromLengths[genome][x[2]])
                        < windowSize
                        and cutoffCatch(
                            x,
                            imds_corrected,
                            regions,
                            catch(x, regions, chromLengths, genome),
                            distance_cut,
                            chromLengths[genome],
                            windowSize,
                        )
                    )
                ]  # int(x[0]) < imds_corrected[x[1]][regions[catch(x, regions, chromLengths, genome, imds_corrected[x[1]])]]))]
                nonClustered_muts = [
                    x[1:] + [x[0]]
                    for x in distances_orig_all
                    if x[1:] + [x[0]] not in clustered_muts
                ]
            else:
                clustered_muts = [
                    x[1:] + [x[0]]
                    for x in distances_orig_all
                    if int(x[0]) <= distance_cut
                ]
                nonClustered_muts = [
                    x[1:] + [x[0]]
                    for x in distances_orig_all
                    if x[1:] + [x[0]] not in clustered_muts
                ]
        else:
            clustered_muts = [
                x[1:] + [x[0]] for x in distances_orig_all if int(x[0]) <= distance_cut
            ]
            nonClustered_muts = [
                x[1:] + [x[0]]
                for x in distances_orig_all
                if x[1:] + [x[0]] not in clustered_muts
            ]
    else:
        clustered_muts = [
            x[1:] for x in distances_orig_all if int(x[0]) <= distance_cut
        ]
        nonClustered_muts = [
            x[1:] for x in distances_orig_all if int(x[0]) > distance_cut
        ]

    if not correctionData:
        with open(vcf_path_clust + project + "_clustered.txt", "a") as clust:
            for muts in clustered_muts:
                sample = muts[0]
                chrom = muts[1]
                pos = muts[2]
                ref = muts[3]
                alt = muts[4]
                plotIMD = muts[5]
                imd_recorded = muts[-1]
                if clustering_vaf:
                    vaf = muts[5]
                    print(
                        "\t".join(
                            [
                                project,
                                sample,
                                ".",
                                genome,
                                "SNP",
                                chrom,
                                pos,
                                pos,
                                ref,
                                alt,
                                "SOMATIC",
                                vaf,
                                plotIMD,
                                imd_recorded,
                            ]
                        ),
                        file=clust,
                    )
                else:
                    print(
                        "\t".join(
                            [
                                project,
                                sample,
                                ".",
                                genome,
                                "SNP",
                                chrom,
                                pos,
                                pos,
                                ref,
                                alt,
                                "SOMATIC",
                                plotIMD,
                                imd_recorded,
                            ]
                        ),
                        file=clust,
                    )

        with open(vcf_path_nonClust + project + "_nonClustered.txt", "a") as nonclust:
            for muts in nonClustered_muts:
                sample = muts[0]
                chrom = muts[1]
                pos = muts[2]
                ref = muts[3]
                alt = muts[4]
                plotIMD = muts[5]
                imd_recorded = muts[-1]
                print(
                    "\t".join(
                        [
                            project,
                            sample,
                            ".",
                            genome,
                            "SNP",
                            chrom,
                            pos,
                            pos,
                            ref,
                            alt,
                            "SOMATIC",
                            plotIMD,
                            imd_recorded,
                        ]
                    ),
                    file=nonclust,
                )

    return (
        y2,
        bincenters2,
        q_val,
        interval_line,
        len(clustered_muts),
        avgSim,
        stdevSim,
        distance_cut,
        lower_CI,
        upper_CI,
        lower_CI_refined,
        upper_CI_refined,
        avg_bin_counts,
    )


def catch(x, regions, chromLengths, genome, imds_corrected=None):
    bisectPos = bisect.bisect_left(regions, int(x[3]) + chromLengths[genome][x[2]])
    if bisectPos >= len(regions):
        bisectPos -= 1
    return bisectPos


def cutoffCatch(
    x, imds_corrected, regions, bisectPos, distance_cut, chroms, windowSize
):
    if x[1] not in imds_corrected:
        return int(x[0]) < distance_cut
    else:
        if regions[bisectPos] not in imds_corrected[x[1]]:
            return int(x[0]) < distance_cut
        else:
            if int(x[3]) + chroms[x[2]] < windowSize:
                return int(x[0]) < imds_corrected[x[1]][regions[bisectPos]]
            else:
                return False


def densityCorrectionOriginal(densityMuts, densityMutsSim, binsDensity):
    """
    Searched the genome for regions with a higher mutation density than expected by chance (1.25 times greater)

    Parameters:
               densityMuts	->	distances with cumulative genomic position added (list; print(densityMuts[0]) for an example)
            densityMutsSim	->	distances with cumulative genomic position added across each simulation (list; nested list of lists)
               binsDensity	->	bins that are generated using a given window size (list)

    Returns:
       correctionFolds	->	the cumulative genomic position of regions that have a higher mutation density than what is expected by chance (list)
    """
    hist, bins = np.histogram(densityMuts, bins=binsDensity, density=True)
    hist_muts, bins_muts = np.histogram(densityMuts, bins=binsDensity, density=False)
    histSims = []
    binsSims = []
    for sim in densityMutsSim:
        histSim, binsSim = np.histogram(
            densityMutsSim[sim], bins=binsDensity, density=True
        )
        histSims.append(histSim)
    p_Vals = []
    meanSims = np.mean(histSims, axis=0)
    stdevSims = np.mean(histSims, axis=0)
    for i in range(0, len(hist), 1):
        z, p = z_test(hist[i], meanSims[i], stdevSims[i])
        p_Vals.append(p)
    p_Vals = np.nan_to_num(p_Vals)
    q_vals = sm.fdrcorrection(p_Vals)[1]
    histSimFinal = np.median(histSims, axis=0)
    fold = hist / histSimFinal
    correctionFolds = [
        y for x, y, z in zip(fold, bins[:], hist_muts) if x > 1.25 and z > 50
    ]  # or x < 1/1.25]
    correctionFolds.append(6000000000)
    return correctionFolds


def moving_density(muts, wS, wE):

    return list(y for y in muts if wS <= y < wE)


def flatten(lis):
    for item in lis:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item


def densityCorrection(densityMuts, densityMutsSim, windowSize):
    """
    Searched the genome for regions with a higher mutation density than expected by chance (1.25 times greater)

    Parameters:
               densityMuts	->	distances with cumulative genomic position added (list; print(densityMuts[0]) for an example)
            densityMutsSim	->	distances with cumulative genomic position added across each simulation (list; nested list of lists)
               binsDensity	->	bins that are generated using a given window size (list)

    Returns:
       correctionFolds	->	the cumulative genomic position of regions that have a higher mutation density than what is expected by chance (list)
    """
    sims = random.sample(list(densityMutsSim.keys()), 10)
    densityMutsSimSample = {x: densityMutsSim[x] for x in sims}

    threshold = 9.0
    try:
        maxDistance = (
            max(
                max(densityMuts),
                max([int(densityMutsSimSample[x][-1]) for x in densityMutsSimSample]),
            )
            + windowSize
        )
    except:
        if len(densityMuts) == 0:
            maxDistance = (
                max([int(densityMutsSimSample[x][-1]) for x in densityMutsSimSample])
                + windowSize
            )
        else:
            maxDistance = max(densityMuts) + windowSize
    # maxDistance = max(max(densityMuts), max([int(densityMutsSim[x][-1]) for x in densityMutsSim])) + windowSize
    simDensities = []
    simWindowParitions = []
    slideSize = windowSize
    imdsProbabilities = []
    for x in densityMutsSimSample:
        currentWindow = 0
        densities = []
        while currentWindow < maxDistance:
            currentMuts = moving_density(
                densityMutsSimSample[x], currentWindow, currentWindow + slideSize
            )
            densities.append(len(currentMuts))
            simWindowParitions.append(currentMuts)

            currentWindow += slideSize
        simDensities.append(densities)

    finalSimDensities = np.mean(simDensities, axis=0)
    finalSimDensities[finalSimDensities == 0] = 1

    densities = []
    currentWindow = 0

    while currentWindow < maxDistance:
        currentMuts = moving_density(
            densityMutsSimSample[x], currentWindow, currentWindow + slideSize
        )
        densities.append(len(currentMuts))
        currentWindow += slideSize

    foldChanges = densities / finalSimDensities
    peaks, heights = find_peaks([0] + foldChanges + [0], threshold)
    regions = [slideSize * (x - 1) for x in peaks]
    return regions


def calculateSampleIMDs(
    project,
    folders,
    directory,
    directory_orig,
    vcf_path_clust,
    vcf_path_nonClust,
    original,
    genome,
    windowSize,
    clustering_vaf,
    contexts,
    exomeSuffix,
    chrom_based,
    correction,
    chromLengths,
):
    imds = {}
    y2s = {}
    bincenters2s = {}
    q_values = {}
    interval_lines = {}
    orig_mutations_samps = {}
    lower_CIs = {}
    upper_CIs = {}
    avg_bin_counts_samp = {}
    avg_simCounts = {}
    std_simCounts = {}
    upper_CIs_refined = {}
    lower_CIs_refined = {}

    imds_corrected = {}
    y2s_corrected = {}
    bincenters2s_corrected = {}
    q_values_corrected = {}
    interval_lines_corrected = {}
    orig_mutations_samps_corrected = {}
    lower_CIs_corrected = {}
    upper_CIs_corrected = {}
    avg_bin_counts_samp_corrected = {}
    regions = []
    regionsSamps = {}
    avg_simCounts_corrected = {}
    std_simCounts_corrected = {}
    upper_CIs_refined_corrected = {}
    lower_CIs_refined_corrected = {}

    for folder in folders:
        regionsSamps[folder] = []
        chromosomes = []
        if folder == ".DS_Store_intradistance.txt" or folder == ".DS_Store":
            continue
        sample = folder
        files = os.listdir(directory + sample + "/")

        if chrom_based:
            overall_distances_all = {}
            distances_orig_all = {}
            distances_orig_all_samps = {}
        else:
            overall_distances_all = []
            distances_orig_all = []
            distances_orig_all_samps = []
        if correction:
            densityMuts = []
            densityMutsSim = {}

        if original:  # Calculate data/plots for original sample
            # Gather the distances for the original sample
            try:
                if not os.path.exists(
                    directory_orig
                    + sample
                    + "_"
                    + contexts
                    + exomeSuffix
                    + "_intradistance.txt"
                ):
                    continue
                with open(
                    directory_orig
                    + sample
                    + "_"
                    + contexts
                    + exomeSuffix
                    + "_intradistance.txt"
                ) as f2:
                    for lines in f2:
                        line = lines.strip().split()
                        if int(line[0]) >= 1:
                            if chrom_based:
                                if line[2] not in distances_orig_all_samps:
                                    distances_orig_all_samps[line[2]] = [line]
                                    distances_orig_all[line[2]] = [int(line[0])]
                                else:
                                    distances_orig_all_samps[line[2]].append(line)
                                    distances_orig_all[line[2]].append(int(line[0]))
                            else:
                                distances_orig_all_samps.append(line)
                                distances_orig_all.append(int(line[0]))
                            if correction:
                                densityMuts.append(
                                    int(line[3]) + chromLengths[genome][line[2]]
                                )
            except:
                print(
                    sample
                    + " does not have nearby IDs to one another. Skipping this sample."
                )
                continue

        # Collect simulation distances
        sim_count = len(files)
        for file in files:
            if correction:
                if "_exome_" in file:
                    sim = file.split("_")[-4]
                else:
                    sim = file.split("_")[-3]
                densityMutsSim[sim] = []
            chroms = []
            if file == ".DS_Store":
                continue
            if chrom_based:
                distances = {}
            else:
                distances = []
            with open(directory + sample + "/" + file) as f:
                for lines in f:
                    line = lines.strip().split()
                    if int(line[0]) >= 1:
                        if chrom_based:
                            if line[2] not in chroms:
                                chroms.append(line[2])
                                distances[line[2]] = [int(line[0])]
                            else:
                                distances[line[2]].append(int(line[0]))
                        else:
                            distances.append(int(line[0]))
                        if correction:
                            densityMutsSim[sim].append(
                                int(line[3]) + chromLengths[genome][line[2]]
                            )
            if chrom_based:
                for chrom in chroms:
                    if chrom not in chromosomes:
                        chromosomes.append(chrom)
                    if chrom not in overall_distances_all:
                        overall_distances_all[chrom] = [distances[chrom]]
                    else:
                        overall_distances_all[chrom].append(distances[chrom])
            else:
                overall_distances_all.append(distances)

        # Perform the regional mutation density corrections
        if correction:
            regions = densityCorrection(densityMuts, densityMutsSim, windowSize)
            regionsSamps[folder] = regions
            try:
                densityCorrectDistances = {}
                densityCorrectDistances_samps = {}
                densityCorrectDistancesSims = {}
                for region in regions:
                    densityCorrectDistances[region] = []
                    densityCorrectDistancesSims[region] = []
                    densityCorrectDistances_samps[region] = []
                with open(
                    directory_orig
                    + sample
                    + "_"
                    + contexts
                    + exomeSuffix
                    + "_intradistance.txt"
                ) as f2:
                    for lines in f2:
                        line = lines.strip().split()
                        if int(line[0]) >= 1:
                            position = int(line[3]) + chromLengths[genome][line[2]]
                            try:
                                bisectRegion = regions[
                                    bisect.bisect_left(regions, position)
                                ]
                            except:
                                continue
                            if bisectRegion - position < windowSize:
                                densityCorrectDistances_samps[bisectRegion].append(line)
                                densityCorrectDistances[bisectRegion].append(
                                    int(line[0])
                                )
            except:
                print(
                    sample
                    + " does not have nearby IDs to one another. Skipping this sample."
                )
                continue

            for file in files:
                if file == ".DS_Store":
                    continue
                distances = {}
                for region in regions:
                    distances[region] = []
                with open(directory + sample + "/" + file) as f:
                    for lines in f:
                        line = lines.strip().split()
                        if int(line[0]) >= 1:
                            position = int(line[3]) + chromLengths[genome][line[2]]
                            try:
                                bisectRegion = regions[
                                    bisect.bisect_left(regions, position)
                                ]
                            except:
                                continue
                            if bisectRegion - position < windowSize:
                                distances[bisectRegion].append(int(line[0]))

                for region in regions:
                    densityCorrectDistancesSims[region].append(distances[region])

            for region in regions:
                correctionData = True
                if sample not in y2s_corrected:
                    y2s_corrected[sample] = {}
                    bincenters2s_corrected[sample] = {}
                    q_values_corrected[sample] = {}
                    interval_lines_corrected[sample] = {}
                    orig_mutations_samps_corrected[sample] = {}
                    lower_CIs_corrected[sample] = {}
                    upper_CIs_corrected[sample] = {}
                    avg_bin_counts_samp_corrected[sample] = {}
                    imds_corrected[sample] = {}
                    avg_simCounts_corrected[sample] = {}
                    std_simCounts_corrected[sample] = {}
                    upper_CIs_refined_corrected[sample] = {}
                    lower_CIs_refined_corrected[sample] = {}
                try:
                    if (
                        len(densityCorrectDistancesSims[region]) == 0
                        or len(densityCorrectDistances_samps[region]) == 0
                        or len(densityCorrectDistances[region]) == 0
                    ):
                        continue
                    (
                        y2s_corrected[sample][region],
                        bincenters2s_corrected[sample][region],
                        q_values_corrected[sample][region],
                        interval_lines_corrected[sample][region],
                        orig_mutations_samps_corrected[sample][region],
                        avg_simCounts_corrected[sample][region],
                        std_simCounts_corrected[sample][region],
                        imds_corrected[sample][region],
                        lower_CIs_corrected[sample][region],
                        upper_CIs_corrected[sample][region],
                        lower_CIs_refined_corrected[sample][region],
                        upper_CIs_refined_corrected[sample][region],
                        avg_bin_counts_samp_corrected[sample][region],
                    ) = first_run(
                        densityCorrectDistancesSims[region],
                        densityCorrectDistances_samps[region],
                        densityCorrectDistances[region],
                        vcf_path_clust,
                        vcf_path_nonClust,
                        sample,
                        original,
                        sim_count,
                        project,
                        genome,
                        clustering_vaf,
                        correctionData,
                        chromLengths,
                    )
                except:
                    continue

        correctionData = False
        if chrom_based:
            y2s[sample] = {}
            bincenters2s[sample] = {}
            q_values[sample] = {}
            interval_lines[sample] = {}
            orig_mutations_samps[sample] = {}
            lower_CIs[sample] = {}
            upper_CIs[sample] = {}
            avg_bin_counts_samp[sample] = {}
            imds[sample] = {}
            avg_simCounts[sample] = {}
            std_simCounts[sample] = {}
            lower_CIs_refined[sample] = {}
            upper_CIs_refined[sample] = {}
            for chrom in chromosomes:
                (
                    y2s[sample][chrom],
                    bincenters2s[sample][chrom],
                    q_values[sample][chrom],
                    interval_lines[sample][chrom],
                    orig_mutations_samps[sample][chrom],
                    avg_simCounts[sample][chrom],
                    std_simCounts[sample][chrom],
                    imds[sample][chrom],
                    lower_CIs[sample][chrom],
                    upper_CIs[sample][chrom],
                    lower_CIs_refined[sample][chrom],
                    upper_CIs_refined[sample][chrom],
                    avg_bin_counts_samp[sample][chrom],
                ) = first_run(
                    overall_distances_all[chrom],
                    distances_orig_all_samps[chrom],
                    distances_orig_all[chrom],
                    vcf_path_clust,
                    vcf_path_nonClust,
                    sample,
                    original,
                    sim_count,
                    project,
                    genome,
                    clustering_vaf,
                    correctionData,
                    correction,
                    regions,
                    imds_corrected,
                    windowSize,
                    chromLengths,
                )
        else:
            if (
                len(overall_distances_all) == 0
                or len(distances_orig_all_samps) == 0
                or len(distances_orig_all) == 0
            ):
                continue
            (
                y2s[sample],
                bincenters2s[sample],
                q_values[sample],
                interval_lines[sample],
                orig_mutations_samps[sample],
                avg_simCounts[sample],
                std_simCounts[sample],
                imds[sample],
                lower_CIs[sample],
                upper_CIs[sample],
                lower_CIs_refined[sample],
                upper_CIs_refined[sample],
                avg_bin_counts_samp[sample],
            ) = first_run(
                overall_distances_all,
                distances_orig_all_samps,
                distances_orig_all,
                vcf_path_clust,
                vcf_path_nonClust,
                sample,
                original,
                sim_count,
                project,
                genome,
                clustering_vaf,
                correctionData,
                correction,
                regions,
                imds_corrected,
                windowSize,
                chromLengths,
            )

    return (
        imds,
        y2s,
        bincenters2s,
        q_values,
        interval_lines,
        orig_mutations_samps,
        lower_CIs,
        upper_CIs,
        avg_bin_counts_samp,
        avg_simCounts,
        std_simCounts,
        upper_CIs_refined,
        lower_CIs_refined,
        imds_corrected,
        y2s_corrected,
        bincenters2s_corrected,
        q_values_corrected,
        interval_lines_corrected,
        orig_mutations_samps_corrected,
        lower_CIs_corrected,
        upper_CIs_corrected,
        avg_bin_counts_samp_corrected,
        regions,
        regionsSamps,
        avg_simCounts_corrected,
        std_simCounts_corrected,
        upper_CIs_refined_corrected,
        lower_CIs_refined_corrected,
    )
    # return(imds)


def getResults(result):
    (
        imds,
        y2s,
        bincenters2s,
        q_values,
        interval_lines,
        orig_mutations_samps,
        lower_CIs,
        upper_CIs,
        avg_bin_counts_samp,
        avg_simCounts,
        std_simCounts,
        upper_CIs_refined,
        lower_CIs_refined,
        imds_corrected,
        y2s_corrected,
        bincenters2s_corrected,
        q_values_corrected,
        interval_lines_corrected,
        orig_mutations_samps_corrected,
        lower_CIs_corrected,
        upper_CIs_corrected,
        avg_bin_counts_samp_corrected,
        regions,
        regionsSamps,
        avg_simCounts_corrected,
        std_simCounts_corrected,
        upper_CIs_refined_corrected,
        lower_CIs_refined_corrected,
    ) = result
    imdsFinal.update(imds)
    y2sFinal.update(y2s)
    bincenters2sFinal.update(bincenters2s)
    q_valuesFinal.update(q_values)
    interval_linesFinal.update(interval_lines)
    orig_mutations_sampsFinal.update(orig_mutations_samps)
    lower_CIsFinal.update(lower_CIs)
    upper_CIsFinal.update(upper_CIs)
    avg_bin_counts_sampFinal.update(avg_bin_counts_samp)
    avg_simCountsFinal.update(avg_simCounts)
    std_simCountsFinal.update(std_simCounts)
    upper_CIs_refinedFinal.update(upper_CIs_refined)
    lower_CIs_refinedFinal.update(lower_CIs_refined)

    imds_correctedFinal.update(imds_corrected)
    y2s_correctedFinal.update(y2s_corrected)
    bincenters2s_correctedFinal.update(bincenters2s_corrected)
    q_values_correctedFinal.update(q_values_corrected)
    interval_lines_correctedFinal.update(interval_lines_corrected)
    orig_mutations_samps_correctedFinal.update(orig_mutations_samps_corrected)
    lower_CIs_correctedFinal.update(lower_CIs_corrected)
    upper_CIs_correctedFinal.update(upper_CIs_corrected)
    avg_bin_counts_samp_correctedFinal.update(avg_bin_counts_samp_corrected)
    # print(regions, regionsFinal)
    # regionsFinal += regions
    regionsSampsFinal.update(regionsSamps)
    avg_simCounts_correctedFinal.update(avg_simCounts_corrected)
    std_simCounts_correctedFinal.update(std_simCounts_corrected)
    upper_CIs_refined_correctedFinal.update(upper_CIs_refined_corrected)
    lower_CIs_refined_correctedFinal.update(lower_CIs_refined_corrected)

    # allResults.append(result)


def hotSpotAnalysis(
    project,
    genome,
    contexts,
    simContext,
    ref_dir,
    windowSize,
    processors,
    plotIMDfigure,
    exome=False,
    chromLengths=None,
    binsDensity=None,
    original=False,
    signature=False,
    percentage=False,
    firstRun=False,
    clustering_vaf=False,
    calculateIMD=True,
    chrom_based=False,
    correction=True,
):
    """
    The main, parent function to calculate sample-dependent IMDs across a data set. Generates output data structures and resulting plots.

    Parameters:
                     project	->	user provided project name (string)
                      genome	->	the reference genome used for the given analysis (string)
                    contexts	->	the contexts used to calculate IMDs (string; ie "96")
              simContext	->	the simulated context used for the background model (list of strings; ie ["6144"])
                     ref_dir	->	the directory for the given project (string)
              windowSize	->	the window size used to calculate the mutation densities across the genome (integer; default=None)
       plotIMDfigure	->	optional parameter that generates IMD and mutational spectra plots for each sample (boolean; default=True).
            chromLengths	->	a dictionary of the cumulative chromosome lengths for a given reference genome (dictionary)
             binsDensity	->	bins that are generated using a given window size (list)
                    original	->	option to plot original sample or not (boolean; default=True)
               signature	->	option to generate plots for signatures as percentages (boolean; default=False)
              percentage	->	option to generate normalized plots (boolean; default=False)
                    firstRun	->	optional parameter for debugging purposes. Always set to True by default (boolean)
      clustering_vaf	-> 	optional paramater to save the VAF from the original mutation files. Required if performing subclassification (default=False; switched to True if subclassify=True)
            calculateIMD	->	optional parameter to skip or calculate IMD. (boolean; default=True)
             chrom_based	->	option to generate IMDs per chromosome (boolean; default=False). This parameter is deprecated. Please DO NOT use.
              correction	->	optional parameter to perform a genome-wide mutational density correction (boolean; default=False)


    Returns:
            regionsSamps	->	a dictionary that contains all of the regions used for calculating corrected IMDs. If correction=False, then it returns an empty datastructure
      imds_corrected	->	a dictionary of all of the corrected IMDs. If correction=False, then it returns an empty datastructure
    """

    ################################################################
    # Organize path suffixes, paths, and several plotting parameters
    ################################################################
    path_suffix = ""
    if correction:
        path_suffix = "_corrected"
    height = 8
    width = 13
    scaled_width = (width / 1.75 * 0.95) / width
    scaled_height = (height / 4.5) / height

    if ref_dir[-1] != "/":
        ref_dir += "/"

    simContext = sorted(simContext, reverse=True)
    simContext = "_".join(simContext)
    file_context = contexts

    if contexts == "96":
        matrix_file_suffix = ".SBS96."
    elif contexts == "INDEL" or contexts == "ID":
        matrix_file_suffix = ".ID83."

    if exome:
        matrix_file_suffix_original = matrix_file_suffix + "exome"
        exomeSuffix = "_exome"
    else:
        matrix_file_suffix_original = matrix_file_suffix + "all"
        exomeSuffix = ""

    original_vcf = (
        ref_dir + "output/vcf_files" + path_suffix + "/single/" + project + "_all.txt"
    )
    directory = (
        ref_dir
        + "output/simulations/"
        + project
        + "_intradistance_"
        + genome
        + "_"
        + contexts
        + exomeSuffix
        + "/"
    )
    directory_out = ref_dir + "output/plots/"
    directory_orig = (
        ref_dir
        + "output/simulations/"
        + project
        + "_intradistance_original_"
        + genome
        + "_"
        + contexts
        + exomeSuffix
        + "/"
    )

    if file_context == "96":
        vcf_path_clust = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_clustered/SNV/"
        )
        vcf_path_nonClust = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_nonClustered/SNV/"
        )
        vcf_path_all_clust = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_all_clustered/SNV/"
        )
        vcf_path_all_nonClust = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_all_nonClustered/SNV/"
        )
        matrix_path = ref_dir + "output/SBS/" + project + matrix_file_suffix_original
        matrix_path_clustered = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_clustered/SNV/output/SBS/"
            + project
            + "_clustered"
            + matrix_file_suffix
            + "all"
        )
        matrix_path_nonClustered = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_nonClustered/SNV/output/SBS/"
            + project
            + "_nonClustered"
            + matrix_file_suffix
            + "all"
        )
        matrix_path_all_clustered = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_all_clustered/SNV/output/SBS"
            + project
            + "_all_clustered"
            + matrix_file_suffix
            + "all"
        )
        matrix_path_all_nonClustered = (
            ref_dir
            + "/references/matrix/"
            + project
            + "_all_nonClustered/"
            + project
            + "_all_nonClustered"
            + matrix_file_suffix
            + "all"
        )
        output_path = (
            directory_out
            + project
            + "_intradistance_plots_"
            + contexts
            + "_"
            + path_suffix
            + ".pdf"
        )

    else:
        vcf_path_clust = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_clustered/INDEL/"
        )
        vcf_path_nonClust = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_nonClustered/INDEL/"
        )
        vcf_path_all_clust = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_all_clustered/INDEL/"
        )
        vcf_path_all_nonClust = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_all_nonClustered/INDEL/"
        )
        matrix_path = ref_dir + "output/ID/" + project + matrix_file_suffix + "all"
        matrix_path_clustered = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_clustered/INDEL/output/ID/"
            + project
            + "_clustered"
            + matrix_file_suffix
            + "all"
        )
        matrix_path_nonClustered = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_nonClustered/INDEL/output/ID/"
            + project
            + "_nonClustered"
            + matrix_file_suffix
            + "all"
        )
        matrix_path_all_clustered = (
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_all_clustered/INDEL/output/ID"
            + project
            + "_all_clustered"
            + matrix_file_suffix
            + "all"
        )
        matrix_path_all_nonClustered = (
            ref_dir
            + "/references/matrix/"
            + project
            + "_all_nonClustered/"
            + project
            + "_all_nonClustered"
            + matrix_file_suffix
            + "all"
        )
        output_path = (
            directory_out
            + project
            + "_intradistance_plots_"
            + contexts
            + path_suffix
            + ".pdf"
        )

    if os.path.exists(directory_out) == False:
        os.mkdir(directory_out)

    if firstRun and calculateIMD:
        if os.path.exists(vcf_path_clust):
            shutil.rmtree(vcf_path_clust)

        if os.path.exists(vcf_path_nonClust):
            shutil.rmtree(vcf_path_nonClust)

    if not os.path.exists(vcf_path_clust):
        os.makedirs(vcf_path_clust)

    if not os.path.exists(vcf_path_nonClust):
        os.makedirs(vcf_path_nonClust)
    ################################################################

    # Instantiate the majority of the data structures
    global imdsFinal, y2sFinal, bincenters2sFinal, q_valuesFinal, interval_linesFinal, orig_mutations_sampsFinal, lower_CIsFinal, upper_CIsFinal, avg_bin_counts_sampFinal, avg_simCountsFinal, std_simCountsFinal, upper_CIs_refinedFinal, lower_CIs_refinedFinal, imds_correctedFinal, y2s_correctedFinal, bincenters2s_correctedFinal, q_values_correctedFinal, interval_lines_correctedFinal, orig_mutations_samps_correctedFinal, lower_CIs_correctedFinal, upper_CIs_correctedFinal, avg_bin_counts_samp_correctedFinal, regionsSampsFinal, avg_simCounts_correctedFinal, std_simCounts_correctedFinal, upper_CIs_refined_correctedFinal, lower_CIs_refined_correctedFinal
    # imds = {}
    # y2s = {}
    # bincenters2s = {}
    # q_values = {}
    # interval_lines = {}
    # orig_mutations_samps = {}
    # lower_CIs = {}
    # upper_CIs = {}
    # avg_bin_counts_samp = {}
    # avg_simCounts = {}
    # std_simCounts = {}
    # upper_CIs_refined = {}
    # lower_CIs_refined = {}

    # imds_corrected = {}
    # y2s_corrected = {}
    # bincenters2s_corrected = {}
    # q_values_corrected = {}
    # interval_lines_corrected = {}
    # orig_mutations_samps_corrected = {}
    # lower_CIs_corrected = {}
    # upper_CIs_corrected = {}
    # avg_bin_counts_samp_corrected = {}
    # regions = []
    # regionsSamps = {}
    # avg_simCounts_corrected = {}
    # std_simCounts_corrected = {}
    # upper_CIs_refined_corrected = {}
    # lower_CIs_refined_corrected = {}

    imdsFinal = {}
    y2sFinal = {}
    bincenters2sFinal = {}
    q_valuesFinal = {}
    interval_linesFinal = {}
    orig_mutations_sampsFinal = {}
    lower_CIsFinal = {}
    upper_CIsFinal = {}
    avg_bin_counts_sampFinal = {}
    avg_simCountsFinal = {}
    std_simCountsFinal = {}
    upper_CIs_refinedFinal = {}
    lower_CIs_refinedFinal = {}

    imds_correctedFinal = {}
    y2s_correctedFinal = {}
    bincenters2s_correctedFinal = {}
    q_values_correctedFinal = {}
    interval_lines_correctedFinal = {}
    orig_mutations_samps_correctedFinal = {}
    lower_CIs_correctedFinal = {}
    upper_CIs_correctedFinal = {}
    avg_bin_counts_samp_correctedFinal = {}
    # regionsFinal = []
    regionsSampsFinal = {}
    avg_simCounts_correctedFinal = {}
    std_simCounts_correctedFinal = {}
    upper_CIs_refined_correctedFinal = {}
    lower_CIs_refined_correctedFinal = {}

    folders = os.listdir(directory)

    ############################################################################################################
    # Start calculating the IMDs for each sample. This includes IMD corrections and calculations for both the
    # original sample and the simulated samples.
    ############################################################################################################
    if firstRun and calculateIMD:
        if os.path.exists(vcf_path_clust + project + "_clustered.txt"):
            os.remove(vcf_path_clust + project + "_clustered.txt")
        if os.path.exists(vcf_path_nonClust + project + "_nonClustered.txt"):
            os.remove(vcf_path_nonClust + project + "_nonClustered.txt")
        if os.path.exists(vcf_path_all_clust + project + "_all_clustered.txt"):
            os.remove(vcf_path_all_clust + project + "_all_clustered.txt")
        if os.path.exists(vcf_path_all_nonClust + project + "_all_nonClustered.txt"):
            os.remove(vcf_path_all_nonClust + project + "_all_nonClustered.txt")

        with open(vcf_path_clust + project + "_clustered.txt", "a") as clust:
            print(
                "\t".join(
                    [
                        x
                        for x in [
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
                            "probability",
                        ]
                    ]
                ),
                file=clust,
            )
        with open(vcf_path_nonClust + project + "_nonClustered.txt", "a") as nonclust:
            print(
                "\t".join(
                    [
                        x
                        for x in [
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
                            "probability",
                        ]
                    ]
                ),
                file=nonclust,
            )

        print(
            "Determining sample-dependent intermutational distance (IMD) cutoff...",
            end="",
            flush=True,
        )

        #####################################################################################################################
        # DONE - This section of code should be put into a function so that we can parallelize it
        #####################################################################################################################
        numSamples = len(folders)
        if numSamples < processors:
            max_seed = numSamples
        else:
            max_seed = processors
        pool = mp.Pool(max_seed)
        samples_parallel = [[] for i in range(max_seed)]
        pool_bin = 0
        for file in folders:
            if pool_bin == max_seed:
                pool_bin = 0
            samples_parallel[pool_bin].append(file)
            pool_bin += 1
        results = []
        global allResults
        allResults = []
        for i in range(0, len(samples_parallel), 1):
            r = pool.apply_async(
                calculateSampleIMDs,
                callback=getResults,
                args=(
                    project,
                    samples_parallel[i],
                    directory,
                    directory_orig,
                    vcf_path_clust,
                    vcf_path_nonClust,
                    original,
                    genome,
                    windowSize,
                    clustering_vaf,
                    contexts,
                    exomeSuffix,
                    chrom_based,
                    correction,
                    chromLengths,
                ),
            )
            results.append(r)
        pool.close()
        pool.join()
        for r in results:
            r.wait()
            if not r.successful():
                # Raises an error when not successful
                r.get()

        #####################################################################################################################
        #####################################################################################################################
        print("Completed!", flush=True)
        ############################################################################################################

        # print("YESSSSS")
        # Generate matrices for all clustered mutations and for all non-clustered mutations
        print("\nAnalyzing clustered mutations...", flush=True)
        matrices = matGen.SigProfilerMatrixGeneratorFunc(
            project + "_clustered", genome, vcf_path_clust, plot=False
        )
        ##To handle the none matrices
        if matrices.get("96") is None:
            print("No clustered mutations found")
            sys.exit()
        elif matrices["96"].sum().sum() == 0:
            print("No clustered mutations found")
            sys.exit()
        # if contexts == "96":
        #     if "96" in matrices:
        #         if matrices["96"].sum().sum() == 0:
        #             print("No clustered mutations found")
        #             sys.exit()
        #     else:
        #         print("No clustered mutations found")
        #         sys.exit()
        elif contexts == "ID":
            print()
            if "ID" in matrices:
                if matrices["ID"].sum().sum() == 0:
                    print("No clustered mutations found")
                    sys.exit()
            else:
                print("No clustered mutations found")
                sys.exit()
        print("\nAnalyzing non-clustered mutations...", flush=True)
        matGen.SigProfilerMatrixGeneratorFunc(
            project + "_nonClustered", genome, vcf_path_nonClust, plot=False
        )
        print(flush=True)

        ############################################################################################################
        # Save all of the relevant data including the imds, imds_corrected, the mutation counts per bin, etc.
        ############################################################################################################
        if not os.path.exists(ref_dir + "output/simulations/data/"):
            os.makedirs(ref_dir + "output/simulations/data/")
        pickleSuffix = ""
        if chrom_based:
            pickleSuffix = "_chrom"
        if correction:
            with open(
                ref_dir
                + "output/simulations/data/imds_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    imds_correctedFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
                )
            with open(
                ref_dir
                + "output/simulations/data/IMDBinHeights_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    y2s_correctedFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
                )
            with open(
                ref_dir
                + "output/simulations/data/IMDBins_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    bincenters2s_correctedFinal,
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            with open(
                ref_dir
                + "output/simulations/data/qvalues_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    q_values_correctedFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
                )
            with open(
                ref_dir
                + "output/simulations/data/interval_lines_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    interval_lines_correctedFinal,
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            with open(
                ref_dir
                + "output/simulations/data/orig_mutations_samps_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    orig_mutations_samps_correctedFinal,
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            with open(
                ref_dir
                + "output/simulations/data/lower_CIs_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    lower_CIs_correctedFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
                )
            with open(
                ref_dir
                + "output/simulations/data/upper_CIs_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    upper_CIs_correctedFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
                )
            with open(
                ref_dir
                + "output/simulations/data/avg_bin_counts_samp_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(
                    avg_bin_counts_samp_correctedFinal,
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )
            with open(
                ref_dir
                + "output/simulations/data/regions_corrected"
                + pickleSuffix
                + ".pickle",
                "wb",
            ) as handle:
                pickle.dump(regionsSampsFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)

        with open(
            ref_dir + "output/simulations/data/imds" + pickleSuffix + ".pickle", "wb"
        ) as handle:
            pickle.dump(imdsFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir
            + "output/simulations/data/IMDBinHeights"
            + pickleSuffix
            + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(y2sFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir + "output/simulations/data/IMDBins" + pickleSuffix + ".pickle", "wb"
        ) as handle:
            pickle.dump(bincenters2sFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir + "output/simulations/data/qvalues" + pickleSuffix + ".pickle", "wb"
        ) as handle:
            pickle.dump(q_valuesFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir
            + "output/simulations/data/interval_lines"
            + pickleSuffix
            + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(interval_linesFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir
            + "output/simulations/data/orig_mutations_samps"
            + pickleSuffix
            + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(
                orig_mutations_sampsFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
            )
        with open(
            ref_dir + "output/simulations/data/lower_CIs" + pickleSuffix + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(lower_CIsFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir + "output/simulations/data/upper_CIs" + pickleSuffix + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(upper_CIsFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir
            + "output/simulations/data/avg_bin_counts_samp"
            + pickleSuffix
            + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(
                avg_bin_counts_sampFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
            )
        with open(
            ref_dir
            + "output/simulations/data/avgSimCounts_samp"
            + pickleSuffix
            + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(avg_simCountsFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir
            + "output/simulations/data/stdSimCounts_samp"
            + pickleSuffix
            + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(std_simCountsFinal, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(
            ref_dir
            + "output/simulations/data/lower_CIs_refined"
            + pickleSuffix
            + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(
                lower_CIs_refinedFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
            )
        with open(
            ref_dir
            + "output/simulations/data/upper_CIs_refined"
            + pickleSuffix
            + ".pickle",
            "wb",
        ) as handle:
            pickle.dump(
                upper_CIs_refinedFinal, handle, protocol=pickle.HIGHEST_PROTOCOL
            )
    ############################################################################################################

    ############################################################################################################
    # If calculateIMD==False, you can skip the previous steps and simply load the saved data for the remainder of
    # the analysis.
    ############################################################################################################
    if not calculateIMD:
        pickleSuffix = ""
        if chrom_based:
            pickleSuffix = "_chrom"
        if correction:
            with open(
                ref_dir
                + "output/simulations/data/imds_corrected"
                + pickleSuffix
                + ".pickle",
                "rb",
            ) as handle:
                imds_correctedFinal = pickle.load(handle)
            with open(
                ref_dir
                + "output/simulations/data/regions_corrected"
                + pickleSuffix
                + ".pickle",
                "rb",
            ) as handle:
                regionsSampsFinal = pickle.load(handle)
        with open(
            ref_dir + "output/simulations/data/imds" + pickleSuffix + ".pickle", "rb"
        ) as handle:
            imdsFinal = pickle.load(handle)
        with open(
            ref_dir
            + "output/simulations/data/IMDBinHeights"
            + pickleSuffix
            + ".pickle",
            "rb",
        ) as handle:
            y2sFinal = pickle.load(handle)
        with open(
            ref_dir + "output/simulations/data/IMDBins" + pickleSuffix + ".pickle", "rb"
        ) as handle:
            bincenters2sFinal = pickle.load(handle)
        with open(
            ref_dir + "output/simulations/data/qvalues" + pickleSuffix + ".pickle", "rb"
        ) as handle:
            q_valuesFinal = pickle.load(handle)
        with open(
            ref_dir
            + "output/simulations/data/interval_lines"
            + pickleSuffix
            + ".pickle",
            "rb",
        ) as handle:
            interval_linesFinal = pickle.load(handle)
        with open(
            ref_dir
            + "output/simulations/data/orig_mutations_samps"
            + pickleSuffix
            + ".pickle",
            "rb",
        ) as handle:
            orig_mutations_sampsFinal = pickle.load(handle)
        with open(
            ref_dir + "output/simulations/data/lower_CIs" + pickleSuffix + ".pickle",
            "rb",
        ) as handle:
            lower_CIsFinal = pickle.load(handle)
        with open(
            ref_dir + "output/simulations/data/upper_CIs" + pickleSuffix + ".pickle",
            "rb",
        ) as handle:
            upper_CIsFinal = pickle.load(handle)
        with open(
            ref_dir
            + "output/simulations/data/avg_bin_counts_samp"
            + pickleSuffix
            + ".pickle",
            "rb",
        ) as handle:
            avg_bin_counts_sampFinal = pickle.load(handle)
        with open(
            ref_dir
            + "output/simulations/data/avgSimCounts_samp"
            + pickleSuffix
            + ".pickle",
            "rb",
        ) as handle:
            avg_simCountsFinal = pickle.load(handle)
        with open(
            ref_dir
            + "output/simulations/data/stdSimCounts_samp"
            + pickleSuffix
            + ".pickle",
            "rb",
        ) as handle:
            std_simCountsFinal = pickle.load(handle)
        with open(
            ref_dir
            + "output/simulations/data/lower_CIs_refined"
            + pickleSuffix
            + ".pickle",
            "rb",
        ) as handle:
            lower_CIs_refinedFinal = pickle.load(handle)
        with open(
            ref_dir
            + "output/simulations/data/upper_CIs_refined"
            + pickleSuffix
            + ".pickle",
            "rb",
        ) as handle:
            upper_CIs_refinedFinal = pickle.load(handle)
    ############################################################################################################

    # Recollect samples names
    if contexts == "96":
        with open(
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_clustered/SNV/output/SBS/"
            + project
            + "_clustered.SBS6.all"
        ) as f:
            first_line = f.readline()
            samples = first_line.strip().split()
            samples = samples[1:]
    elif contexts == "INDEL" or contexts == "ID":
        with open(
            ref_dir
            + "output/vcf_files"
            + path_suffix
            + "/"
            + project
            + "_clustered/INDEL/output/ID/"
            + project
            + "_clustered.ID83.all"
        ) as f:
            first_line = f.readline()
            samples = first_line.strip().split()
            samples = samples[1:]

    ################################
    # Generate the IMD/spectra plots
    ################################
    if plotIMDfigure:
        if exome:
            simContext += "_exome"
        pp = PdfPages(
            directory_out
            + project
            + "_intradistance_plots_"
            + simContext
            + path_suffix
            + ".pdf"
        )
        histo = True
        print("Plotting SigProfilerClusters Results...", end="", flush=True)
        for folder in folders:
            if folder == ".DS_Store_intradistance.txt" or folder == ".DS_Store":
                continue
            if folder not in samples:
                histo = False
            sample = folder
            files = os.listdir(directory + sample + "/")
            if not chrom_based:
                fig = plt.figure(figsize=(width, height))
                panel1 = plt.axes(
                    [0.075, 0.225 + scaled_height * 2, scaled_width, scaled_height]
                )
                panel2 = plt.axes(
                    [
                        0.125 + scaled_width,
                        0.225 + scaled_height * 2,
                        0.3,
                        scaled_height,
                    ]
                )
                panel3 = plt.axes(
                    [0.075, 0.15 + scaled_height, scaled_width, scaled_height]
                )
                panel4 = plt.axes(
                    [0.125 + scaled_width, 0.15 + scaled_height, 0.3, scaled_height]
                )
                panel5 = plt.axes([0.075, 0.075, scaled_width, scaled_height])
                panel6 = plt.axes([0.125 + scaled_width, 0.075, 0.3, scaled_height])

                if histo:
                    if not chrom_based:
                        clustered = plot_hist(
                            y2sFinal[sample],
                            bincenters2sFinal[sample],
                            q_valuesFinal[sample],
                            interval_linesFinal[sample],
                            orig_mutations_sampsFinal[sample],
                            avg_simCountsFinal[sample],
                            std_simCountsFinal[sample],
                            imdsFinal[sample],
                            lower_CIsFinal[sample],
                            upper_CIsFinal[sample],
                            lower_CIs_refinedFinal[sample],
                            upper_CIs_refinedFinal[sample],
                            avg_bin_counts_sampFinal[sample],
                            sample,
                            original,
                            panel2,
                            panel3,
                            panel4,
                            panel5,
                            panel6,
                        )
                        if clustered:
                            if file_context == "96":
                                plottingFunctions.plot96_same(
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
                                )
                            else:
                                plottingFunctions.plotINDEL_same(
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
                                )
                            pp.savefig()
                        plt.close()
                histo = True

        plt.close()
        pp.close()
    print("Completed!\n", flush=True)
    ################################

    return (regionsSampsFinal, imds_correctedFinal)
