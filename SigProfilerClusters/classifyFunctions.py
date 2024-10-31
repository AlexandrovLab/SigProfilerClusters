import os
import numpy as np
import os
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerClusters import hotspot
import shutil
from numpy import median
import pickle
import bisect
import re
import multiprocessing as mp
import sys


def processivitySubclassification(event, out2Y, out2K, out2S, out2N):
    """
    Separates Class 2 kataeigs events into different subclassifications of processivity.
    Specifically, 	Class2Y are when 75% of mutations in an event are processive pyrimidines (C and T or A and G)
                                    Class2K are when 75% of mutations in an event are processive keto bases (T and G or A and C)
                                    Class2S are when 75% of mutations in an event are processive strong bases (C and G or A and T)
                                    Class2N are for all remaining class 2 events that do not fall into one of the other groups.

    Parameters:
            event	->	the mutations present in a given event (list)
            out2Y	->	path to write the class2Y events in (string)
            out2K	->	path to write the class2K events in (string)
            out2S	->	path to write the class2S events in (string)
            out2N	->	path to write the class2N events in (string)

    Returns:
            None

    Outputs:
            Writes the subclassifed class2 events into a respective text file.
    """
    refs = "".join([x[8] for x in event])
    pyrProcessCG = max(
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join([x for x in ["1" if x == "C" else "0" for x in refs]])
                )
            )
        ),
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join([x for x in ["1" if x == "G" else "0" for x in refs]])
                )
            )
        ),
    )
    pyrProcessTA = max(
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join([x for x in ["1" if x == "T" else "0" for x in refs]])
                )
            )
        ),
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join([x for x in ["1" if x == "A" else "0" for x in refs]])
                )
            )
        ),
    )
    purityCG = max(refs.count("G") / len(refs), refs.count("C") / len(refs))
    purityTA = max(refs.count("T") / len(refs), refs.count("A") / len(refs))
    pyrProcess = max(
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join(
                        [
                            x
                            for x in [
                                "1" if x == "C" or x == "T" else "0" for x in refs
                            ]
                        ]
                    )
                )
            )
        ),
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join(
                        [
                            x
                            for x in [
                                "1" if x == "G" or x == "A" else "0" for x in refs
                            ]
                        ]
                    )
                )
            )
        ),
    )
    ketoProcess = max(
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join(
                        [
                            x
                            for x in [
                                "1" if x == "T" or x == "G" else "0" for x in refs
                            ]
                        ]
                    )
                )
            )
        ),
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join(
                        [
                            x
                            for x in [
                                "1" if x == "A" or x == "C" else "0" for x in refs
                            ]
                        ]
                    )
                )
            )
        ),
    )
    strongProcess = max(
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join(
                        [
                            x
                            for x in [
                                "1" if x == "C" or x == "G" else "0" for x in refs
                            ]
                        ]
                    )
                )
            )
        ),
        len(
            max(
                re.compile("(1+1)*").findall(
                    "".join(
                        [
                            x
                            for x in [
                                "1" if x == "A" or x == "T" else "0" for x in refs
                            ]
                        ]
                    )
                )
            )
        ),
    )

    maxProcessivity = max(pyrProcess, ketoProcess, strongProcess)
    processiveLenthRequirement = len(refs) * 0.75
    if (
        pyrProcessCG > processiveLenthRequirement
        or pyrProcessTA > processiveLenthRequirement
        or purityCG > 0.9
        or purityTA > 0.9
    ):
        for line in event:
            print("\t".join([x for x in line]), file=out2Y)
    else:
        if (
            maxProcessivity == pyrProcess
            and maxProcessivity >= processiveLenthRequirement
        ):
            for line in event:
                print("\t".join([x for x in line]), file=out2Y)
        elif (
            maxProcessivity == ketoProcess
            and maxProcessivity >= processiveLenthRequirement
        ):
            for line in event:
                print("\t".join([x for x in line]), file=out2K)
        elif (
            maxProcessivity == strongProcess
            and maxProcessivity >= processiveLenthRequirement
        ):
            for line in event:
                print("\t".join([x for x in line]), file=out2S)
        else:
            for line in event:
                print("\t".join([x for x in line]), file=out2N)


def pullCCF(project, project_path, correction=True):
    """
    Collects the CCFs from the original mutation files. Assumes that these are provided in the last column of each mutation line

    Parameters:
                     project	->	user provided project name (string)
            project_path	->	the directory for the given project (string)
              correction	->	optional parameter to perform a genome-wide mutational density correction (boolean; default=False)

    Returns:
            None

    Outputs:
            New clustered mutation files that contain the CCF for each mutation.
    """
    path_suffix = ""
    if correction:
        path_suffix += "_corrected"
    vcf_path = project_path + "input/"
    clusteredMutsPath = (
        project_path + "output/vcf_files" + path_suffix + "/" + project + "_clustered/"
    )
    clusteredMutsFile = (
        project_path
        + "output/vcf_files"
        + path_suffix
        + "/"
        + project
        + "_clustered"
        + "/SNV/"
        + project
        + "_clustered.txt"
    )
    vcf_files = [x for x in os.listdir(vcf_path) if x != ".DS_Store"]

    ccfs = {}
    for vcfFile in vcf_files:
        sample = vcfFile.split(".")[0]
        ccfs[sample] = {}
        with open(vcf_path + vcfFile) as f:
            for lines in f:
                if lines[0] == "#":
                    continue
                lines = lines.strip().split()
                chrom = lines[0]
                if len(chrom) > 3:
                    chrom = chrom[3:]
                pos = lines[1]
                ref = lines[3]
                alt = lines[4]
                try:
                    ccf = float(lines[-1])
                except:
                    ccf = -1.5
                    # print("There does not seem to be CCF scores in this input file.\n\t", vcfFile)
                    # break
                if len(ref) == len(alt) and len(ref) > 1:
                    for i in range(len(ref)):
                        keyLine = ":".join([chrom, str(int(pos) + i), ref[i], alt[i]])
                        ccfs[sample][keyLine] = ccf
                else:
                    keyLine = ":".join([chrom, pos, ref, alt])
                    ccfs[sample][keyLine] = ccf

    with open(clusteredMutsFile) as f, open(
        clusteredMutsPath + project + "_clustered_vaf.txt", "w"
    ) as out:
        next(f)
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
                        "VAF/CCF",
                    ]
                ]
            ),
            file=out,
        )
        for lines in f:
            lines = lines.strip().split()
            sample = lines[1]
            newKey = ":".join([lines[5], lines[6], lines[8], lines[9]])

            try:
                ccf = ccfs[sample][newKey]
            except:
                newKey = ":".join(["chr" + lines[5], lines[6], lines[8], lines[9]])
                ccf = ccfs[sample][newKey]
            lines.append(str(ccf))
            print("\t".join([x for x in lines]), file=out)


def pullVaf(project, project_path, variant_caller=None, correction=True):
    """
    Collects the VAFs from the original mutation files. Assumes that these are provided in the same
    format as Sanger or TCGA.

    Parameters:
                     project	->	user provided project name (string)
            project_path	->	the directory for the given project (string)
            variant_caller	->	optional parameter that informs the tool of what format the VAF scores are provided (boolean; default=None). This is required when subClassify=True. Currently, there are four supported formats:
                                -> sanger: If your VAF is recorded in the 11th column of your VCF as the last number of the colon delimited values, set variant_caller="sanger".
                                -> TCGA: If your VAF is recorded in the 8th column of your VCF as VCF=xx, set variant_caller="TCGA".
                                -> standardVC: If your VAF is recorded in the 10th column of your VCF as AF=xx, set variant_caller="standardVC".
                                -> mutect2: If your VAF is recorded in the 11th column of your VCF as AF=xx, set variant_caller="mutect2".
              correction	->	optional parameter to perform a genome-wide mutational density correction (boolean; default=False)

    Returns:
            None

    Outputs:
            New clustered mutation files that contain the VAF for each mutation.
    """
    path_suffix = ""
    if correction:
        path_suffix += "_corrected"
    vcf_path = project_path + "input/"
    clusteredMutsPath = (
        project_path + "output/vcf_files" + path_suffix + "/" + project + "_clustered/"
    )
    clusteredMutsFile = (
        project_path
        + "output/vcf_files"
        + path_suffix
        + "/"
        + project
        + "_clustered"
        + "/SNV/"
        + project
        + "_clustered.txt"
    )
    vcf_files = [x for x in os.listdir(vcf_path) if x != ".DS_Store"]

    # Dictionary for variant caller mapping
    variant_type_dict = {
        "sanger": "sanger",
        "tcga": "TCGA",
        "standardvc": "standardVC",
        "mutect2": "mutect2",
    }
    # Check if variant_caller is provided
    if variant_caller is None:
        raise ValueError(
            "Please specify your variant caller. Currently we are supporting four different variant callers: sanger, TCGA, standardVC and mutect2."
        )

    # Map variant_caller to standardized value using variant_type_dict
    if variant_caller.lower() in variant_type_dict:
        variant_caller = variant_type_dict[variant_caller.lower()]
    else:
        raise ValueError(
            f"Unsupported variant caller: {variant_caller}. Currently we are supporting four different variant callers: sanger, TCGA, standardVC and mutect2"
        )

    if variant_caller == "sanger":
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
                    if len(chrom) > 3:
                        chrom = chrom[3:]
                    pos = lines[1]
                    ref = lines[3]
                    alt = lines[4]
                    try:
                        vaf = float(lines[-1].split(":")[-1])
                    except:
                        print(
                            "There does not seem to be VAF scores in this Sanger-produced file.\n\t",
                            vcfFile,
                        )
                        break
                    if len(ref) == len(alt) and len(ref) > 1:
                        for i in range(len(ref)):
                            keyLine = ":".join(
                                [chrom, str(int(pos) + i), ref[i], alt[i]]
                            )
                            vafs[sample][keyLine] = vaf
                    else:
                        keyLine = ":".join([chrom, pos, ref, alt])
                        vafs[sample][keyLine] = vaf

    elif variant_caller == "TCGA":
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
                    if len(ref) == len(alt) and len(ref) > 1:
                        for i in range(len(ref)):
                            keyLine = ":".join(
                                [chrom, str(int(pos) + i), ref[i], alt[i]]
                            )
                            vafs[sample][keyLine] = vaf
                    else:
                        keyLine = ":".join([chrom, pos, ref, alt])
                        vafs[sample][keyLine] = vaf

    elif variant_caller == "standardVC":
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
                        vaf = float(lines[9].split("AF=")[1].split(";")[0])
                    except:
                        vaf = -1.5
                    if len(ref) == len(alt) and len(ref) > 1:
                        for i in range(len(ref)):
                            keyLine = ":".join(
                                [chrom, str(int(pos) + i), ref[i], alt[i]]
                            )
                            vafs[sample][keyLine] = vaf
                    else:
                        keyLine = ":".join([chrom, pos, ref, alt])
                        vafs[sample][keyLine] = vaf
    elif variant_caller == "mutect2":
        field = "AF"
        vcf_col = 11
        vafs = {}
        for vcfFile in vcf_files:
            sample = vcfFile.split(".")[0]
            vafs[sample] = {}
            with open(os.path.join(vcf_path, vcfFile)) as f:
                for lines in f:
                    if lines[0] == "#":
                        continue
                    lines = lines.strip().split()
                    chrom = lines[0]
                    if chrom.startswith("chr") or chrom.startswith("Chr"):
                        chrom = chrom[3:]
                    pos = lines[1]
                    ref = lines[3]
                    alt = lines[4]
                    try:
                        ## Column 9 is assumed to be FORMAT
                        fmt = lines[8].split(":")
                        ## Extract VAF index and use it to extract VAF from the provided column
                        vaf_ind = fmt.index(field)
                        ## Use vcf_col-1 because humans think lists as 1-indexed but python thinks 0-indexed
                        vaf = float(lines[vcf_col - 1].split(":")[vaf_ind])
                    except:
                        print(
                            "Provided VAF field does not match any field in VCF.\n\t",
                            vcfFile,
                        )
                        break
                    if len(ref) == len(alt) and len(ref) > 1:
                        for i in range(len(ref)):
                            keyLine = ":".join(
                                [chrom, str(int(pos) + i), ref[i], alt[i]]
                            )
                            vafs[sample][keyLine] = vaf
                    else:
                        keyLine = ":".join([chrom, pos, ref, alt])
                        vafs[sample][keyLine] = vaf

    with open(clusteredMutsFile) as f, open(
        clusteredMutsPath + project + "_clustered_vaf.txt", "w"
    ) as out:
        next(f)
        # print("HEADER", file=out)
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
                        "VAF/CCF",
                    ]
                ]
            ),
            file=out,
        )
        for lines in f:
            lines = lines.strip().split()
            sample = lines[1]
            newKey = ":".join([lines[5], lines[6], lines[8], lines[9]])

            try:
                vaf = vafs[sample][newKey]
            except:
                newKey = ":".join(["chr" + lines[5], lines[6], lines[8], lines[9]])
                vaf = vafs[sample][newKey]
            lines.append(str(vaf))
            print("\t".join([x for x in lines]), file=out)


def findClustersOfClusters(
    project,
    chrom_based,
    project_parent_path,
    windowSize,
    chromLengths,
    regions,
    log_out,
    genome,
    processors,
    imds,
    correction=True,
    includedCCFs=False,
):
    """
    Subclassifies the clustered mutations into different categories including DBS, MBS, omikli, kataegis, and other. This is only performed
    for simple single base substutions. Indels are not subclassified using this scheme.

    Parameters:
                            project	->	user provided project name (string)
                    chrom_based	->	option to generate IMDs per chromosome (boolean; default=False)
    project_parent_path	->	the directory for the given project (string)
                     windowSize	->	the window size used to calculate the mutation densities across the genome (integer; default=None)
               chromLengths	->	a dictionary of the cumulative chromosome lengths for a given reference genome (dictionary)
                            regions	->	a dictionary that contains all of the regions used for calculating corrected IMDs. If correction=False, then it returns an empty datastructure. (dictionary)
                             genome	->	the reference genome used for the given analysis (string)
                               imds	->	a dictionary of all of the corrected IMDs. If correction=False, then it returns an empty datastructure (dictionary)
                     correction	->	optional parameter to perform a genome-wide mutational density correction (boolean; default=False)

    Returns:
            None

    Outputs:
            Subclassification files for all clustered mutations.
    """

    ####################################################################################
    # Organize paths and file names
    ####################################################################################
    path_suffix = ""
    path_suffix2 = ""
    if chrom_based:
        path_suffix = "_chrom"
    if correction:
        path_suffix2 += "_corrected"

    project_path = (
        project_parent_path
        + "output/vcf_files"
        + path_suffix2
        + "/"
        + project
        + "_clustered/"
    )
    file = project_path + project + "_clustered_vaf.txt"
    out_file = project_path + project + "_clusters_of_clusters.txt"
    out_file2 = project_path + project + "_clusters_of_clusters_imd.txt"

    out_file3 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class1/"
        + project
        + "_clustered_class1.txt"
    )
    out_file4 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2/"
        + project
        + "_clustered_class2.txt"
    )
    out_file5 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class3/"
        + project
        + "_clustered_class3.txt"
    )
    out_file6 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class1a/"
        + project
        + "_clustered_class1a.txt"
    )
    out_file7 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class1b/"
        + project
        + "_clustered_class1b.txt"
    )
    out_file8 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class1c/"
        + project
        + "_clustered_class1c.txt"
    )
    out_file9 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2Y/"
        + project
        + "_clustered_class2Y.txt"
    )
    out_file10 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2K/"
        + project
        + "_clustered_class2K.txt"
    )
    out_file11 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2S/"
        + project
        + "_clustered_class2S.txt"
    )
    out_file12 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2N/"
        + project
        + "_clustered_class2N.txt"
    )

    if os.path.exists(project_path + "subclasses" + path_suffix + "/class1/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class1/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class1/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class3/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class3/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class3/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class1a/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class1a/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class1a/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class1b/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class1b/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class1b/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class1c/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class1c/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class1c/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2Y/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2Y/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2Y/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2K/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2K/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2K/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2S/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2S/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2S/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2N/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2N/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2N/")
    ####################################################################################

    # Hard-coded cutoffs, which can later be used as additional parameters to the tool.
    cutoff = 10001  # Min distance required between adjacenet events
    vaf_cut = 0.1  # Max difference in VAFs of adjacent mutations
    if includedCCFs:
        vaf_cut = 0.25  # Max difference in VAFs of adjacent mutations

    # Load the IMD data structures
    if chrom_based:
        with open(
            project_parent_path + "output/simulations/data/imds_chrom.pickle", "rb"
        ) as handle:
            imdsData = pickle.load(handle)
    else:
        with open(
            project_parent_path + "output/simulations/data/imds.pickle", "rb"
        ) as handle:
            imdsData = pickle.load(handle)
    if correction:
        with open(
            project_parent_path + "output/simulations/data/imds_corrected.pickle", "rb"
        ) as handle:
            imds_corrected = pickle.load(handle)

    first_pass = (
        True  # Deprecated option when originally developing. Can probably be deleted
    )
    total_muts = {}
    if first_pass:
        mnv_length = 0
        len_mnvs = {}
        total_mnvs = {}
        distances = []
        count = 1
        out = open(out_file, "w")
        print(
            "\t".join(
                [
                    x
                    for x in [
                        "clust_group",
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
                        "VAF/CCF",
                    ]
                ]
            ),
            file=out,
        )
        with open(file) as f:
            next(f)
            lines = [line.strip().split() for line in f]

        for i in range(1, len(lines), 1):
            prev_chrom = lines[i - 1][5]
            prev_pos = int(lines[i - 1][6])
            prev_samp = lines[i - 1][1]
            chrom = lines[i][5]
            pos = int(lines[i][6])
            samp = lines[i][1]
            if samp not in total_muts:
                total_muts[samp] = 0
            if prev_samp == samp:
                if prev_chrom == chrom:
                    if pos - prev_pos < cutoff or (
                        correction
                        and len(regions[samp]) > 0
                        and (
                            (
                                regions[samp][
                                    hotspot.catch(
                                        [".", ".", chrom, pos],
                                        regions[samp],
                                        chromLengths,
                                        genome,
                                    )
                                ]
                                - (pos + chromLengths[genome][chrom])
                                < windowSize
                            )
                            and hotspot.cutoffCatch(
                                [(pos - prev_pos), samp, chrom, pos],
                                imds_corrected,
                                regions[samp],
                                hotspot.catch(
                                    [".", ".", chrom, pos],
                                    regions[samp],
                                    chromLengths,
                                    genome,
                                ),
                                imdsData[samp],
                                chromLengths[genome],
                                windowSize,
                            )
                        )
                    ):  # (pos - prev_pos) < imds_corrected[samp][regions[samp][hotspot.catch([".",".",chrom, pos], regions[samp], chromLengths, genome, imds_corrected[samp])]])):
                        distances.append(pos - prev_pos)
                        mnv_length += 1
                        lines[i - 1] = [str(count)] + lines[i - 1]
                        print("\t".join([x for x in lines[i - 1]]), file=out)
                        total_muts[samp] += 1
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
                        lines[i - 1] = [str(count)] + lines[i - 1]
                        print("\t".join([x for x in lines[i - 1]]), file=out)
                        total_muts[samp] += 1
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

                    lines[i - 1] = [str(count)] + lines[i - 1]
                    print("\t".join([x for x in lines[i - 1]]), file=out)
                    total_muts[samp] += 1
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

                lines[i - 1] = [str(count)] + lines[i - 1]
                print("\t".join([x for x in lines[i - 1]]), file=out)
                total_muts[samp] += 1
                count += 1
                count = 1
                print("\n\n################ New Sample #################", file=out)
                print("\n\n", file=out)

        lines[i] = [str(count)] + lines[i]
        print("\t".join([x for x in lines[i]]), file=out)
        total_muts[samp] += 1
        out.close()
        with open(
            project_parent_path
            + "output/simulations/data/originalCounts"
            + path_suffix2
            + ".pickle",
            "wb",
        ) as f:
            pickle.dump(total_muts, f)

    if True:
        len_mnvs = {"I": {}, "II": {}, "III": {}, "Ia": {}, "Ib": {}, "Ic": {}}
        distances = []
        distances_mnv = {}
        lines = []
        count = 1
        subclassesHeader = "\t".join(
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
                    "group",
                    "IMD",
                    "VAF/CCF",
                    "subclass",
                ]
            ]
        )

        with open(out_file) as f, open(out_file2, "w") as out2, open(
            out_file3, "w"
        ) as out3, open(out_file4, "w") as out4, open(out_file5, "w") as out5, open(
            out_file6, "w"
        ) as out6, open(
            out_file7, "w"
        ) as out7, open(
            out_file8, "w"
        ) as out8, open(
            out_file9, "w"
        ) as out2Y, open(
            out_file10, "w"
        ) as out2K, open(
            out_file11, "w"
        ) as out2S, open(
            out_file12, "w"
        ) as out2N:
            print(
                "\t".join(
                    [
                        x
                        for x in [
                            "clust_group",
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
                            "VAF/CCF",
                        ]
                    ]
                ),
                file=out2,
            )
            print(subclassesHeader, file=out3)
            print(subclassesHeader, file=out4)
            print(subclassesHeader + "\treason", file=out5)
            print(subclassesHeader, file=out6)
            print(subclassesHeader, file=out7)
            print(subclassesHeader, file=out8)
            print(subclassesHeader, file=out2Y)
            print(subclassesHeader, file=out2K)
            print(subclassesHeader, file=out2S)
            print(subclassesHeader, file=out2N)

            next(f)
            for line in f:
                # line = line.strip().split()[1:]
                line = line.strip().split()
                if line != []:
                    line = line[1:-2] + [line[0]] + line[-2:]
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
                            if correction:
                                distancesLine = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if int(y[7]) - int(x[7]) > 1
                                        and (
                                            int(y[7]) - int(x[7]) <= imdsData[y[1]]
                                            or (
                                                len(regions[y[1]]) > 0
                                                and regions[y[1]][
                                                    hotspot.catch(
                                                        [".", ".", y[5], y[7]],
                                                        regions[y[1]],
                                                        chromLengths,
                                                        genome,
                                                    )
                                                ]
                                                - (
                                                    int(y[7])
                                                    + chromLengths[genome][y[5]]
                                                )
                                                < windowSize
                                                and hotspot.cutoffCatch(
                                                    [int(y[-2]), y[1], y[5], y[7]],
                                                    imds_corrected,
                                                    regions[y[1]],
                                                    hotspot.catch(
                                                        [".", ".", y[5], y[7]],
                                                        regions[y[1]],
                                                        chromLengths,
                                                        genome,
                                                    ),
                                                    imdsData[y[1]],
                                                    chromLengths[genome],
                                                    windowSize,
                                                )
                                            )
                                        )
                                    ]
                                )  #   int(y[-2]) < imds_corrected[y[1]][regions[y[1]][hotspot.catch([".",".",y[5], y[7]], regions[y[1]], chromLengths, genome, imds_corrected[y[1]])]])   )])
                                distancesFailed = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if (
                                            int(y[7]) - int(x[7]) > imdsData[y[1]]
                                            and (
                                                len(regions[y[1]]) > 0
                                                and regions[y[1]][
                                                    hotspot.catch(
                                                        [".", ".", y[5], y[7]],
                                                        regions[y[1]],
                                                        chromLengths,
                                                        genome,
                                                    )
                                                ]
                                                - (
                                                    int(y[7])
                                                    + chromLengths[genome][y[5]]
                                                )
                                                < windowSize
                                                and hotspot.cutoffCatch(
                                                    [int(y[-2]), y[1], y[5], y[7]],
                                                    imds_corrected[y[1]],
                                                    regions[y[1]],
                                                    hotspot.catch(
                                                        [".", ".", y[5], y[7]],
                                                        regions[y[1]],
                                                        chromLengths,
                                                        genome,
                                                    ),
                                                    imdsData[y[1]],
                                                    chromLengths[genome],
                                                    windowSize,
                                                )
                                            )
                                        )
                                    ]
                                )
                                zeroDistances = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if int(y[7]) - int(x[7]) > 1
                                    ]
                                )

                            else:
                                distancesLine = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if int(y[7]) - int(x[7]) > 1
                                        and (int(y[7]) - int(x[7]) <= imdsData[y[1]])
                                    ]
                                )
                                distancesFailed = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if (int(y[7]) - int(x[7]) > imdsData[y[1]])
                                    ]
                                )
                                zeroDistances = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if int(y[7]) - int(x[7]) > 1
                                    ]
                                )

                            vafs = len(
                                [
                                    abs(float(y[-1]) - float(x[-1]))
                                    for x, y in zip(lines, lines[1:])
                                    if round(abs(float(y[-1]) - float(x[-1])), 4)
                                    > vaf_cut
                                ]
                            )
                            unknownVafs = len(
                                [
                                    abs(float(y[-1]) - float(x[-1]))
                                    for x, y in zip(lines, lines[1:])
                                    if abs(float(y[-1]) - float(x[-1])) < 0
                                ]
                            )
                            distancesActual = [
                                int(y[7]) - int(x[7]) for x, y in zip(lines, lines[1:])
                            ]
                            vafsActual = [
                                abs(float(y[-1]) - float(x[-1]))
                                for x, y in zip(lines, lines[1:])
                            ]

                            if unknownVafs == 0 and vafs == 0 and distancesFailed == 0:
                                # Class I or II
                                if distancesLine > 1:
                                    # Class IC or II
                                    if len(lines) >= 4:
                                        writeClassII = True
                                    else:
                                        writeClassI = True
                                        writeClassIc = True
                                else:
                                    writeClassI = True
                                    if zeroDistances == 0:
                                        if len(lines) == 2:
                                            writeClassIa = True
                                        else:
                                            writeClassIb = True
                                    else:
                                        writeClassIc = True
                            else:
                                writeClassIII = True

                        if writeClassII:
                            processivitySubclassification(
                                lines, out2Y, out2K, out2S, out2N
                            )
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
                                        print(
                                            "\t".join([x for x in lines[i]]), file=out3
                                        )
                                except:
                                    print(lines)

                                if writeClassIc:
                                    # Writes Class Ic (Omiklis)
                                    try:
                                        for i in range(0, len(lines), 1):
                                            lines[i][-1] = "ClassIC"
                                            print(
                                                "\t".join([x for x in lines[i]]),
                                                file=out8,
                                            )
                                    except:
                                        print(lines)
                                elif writeClassIa:
                                    # Writes Class Ia (DBSs)
                                    try:
                                        for i in range(0, len(lines), 1):
                                            lines[i][-1] = "ClassIA"
                                            print(
                                                "\t".join([x for x in lines[i]]),
                                                file=out6,
                                            )
                                    except:
                                        print(lines)
                                elif writeClassIb:
                                    # Writes Class 1b (MBSs)
                                    try:
                                        for i in range(0, len(lines), 1):
                                            lines[i][-1] = "ClassIB"
                                            print(
                                                "\t".join([x for x in lines[i]]),
                                                file=out7,
                                            )
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
                                            if correction:
                                                if abs(
                                                    float(linesSubClass[i][-1])
                                                    - float(lineRef[-1])
                                                ) < vaf_cut and (
                                                    int(linesSubClass[i][7])
                                                    - int(lineRef[7])
                                                    <= imdsData[lineRef[1]][lineRef[5]]
                                                    or (
                                                        regions[lineRef[1]][
                                                            bisect.bisect_left(
                                                                regions[lineRef[1]],
                                                                (
                                                                    int(lineRef[7])
                                                                    + chromLengths[
                                                                        genome
                                                                    ][lineRef[5]]
                                                                ),
                                                            )
                                                        ]
                                                        - (
                                                            int(lineRef[7])
                                                            + chromLengths[genome][
                                                                lineRef[5]
                                                            ]
                                                        )
                                                        < windowSize
                                                        and int(lineRef[-2])
                                                        < imds_corrected[lineRef[1]][
                                                            regions[lineRef[1]][
                                                                bisect.bisect_left(
                                                                    regions[lineRef[1]],
                                                                    int(lineRef[7])
                                                                    + chromLengths[
                                                                        genome
                                                                    ][lineRef[5]],
                                                                )
                                                            ]
                                                        ]
                                                    )
                                                ):
                                                    saveNewEvent.append(
                                                        linesSubClass[i]
                                                    )
                                                    lineRef = linesSubClass[i]
                                            else:
                                                if (
                                                    abs(
                                                        float(linesSubClass[i][-1])
                                                        - float(lineRef[-1])
                                                    )
                                                    < vaf_cut
                                                    and int(linesSubClass[i][7])
                                                    - int(lineRef[7])
                                                    <= imdsData[lineRef[1]][lineRef[5]]
                                                ):
                                                    saveNewEvent.append(
                                                        linesSubClass[i]
                                                    )
                                                    lineRef = linesSubClass[i]
                                        else:
                                            if correction:
                                                if abs(
                                                    float(linesSubClass[i][-1])
                                                    - float(lineRef[-1])
                                                ) < vaf_cut and (
                                                    int(linesSubClass[i][7])
                                                    - int(lineRef[7])
                                                    <= imdsData[lineRef[1]]
                                                    or (
                                                        len(regions[lineRef[1]]) > 0
                                                        and regions[lineRef[1]][
                                                            hotspot.catch(
                                                                [
                                                                    ".",
                                                                    ".",
                                                                    lineRef[5],
                                                                    lineRef[7],
                                                                ],
                                                                regions[lineRef[1]],
                                                                chromLengths,
                                                                genome,
                                                            )
                                                        ]
                                                        - (
                                                            int(lineRef[7])
                                                            + chromLengths[genome][
                                                                lineRef[5]
                                                            ]
                                                        )
                                                        < windowSize
                                                        and hotspot.cutoffCatch(
                                                            [
                                                                int(lineRef[-2]),
                                                                lineRef[1],
                                                                lineRef[5],
                                                                lineRef[7],
                                                            ],
                                                            imds_corrected,
                                                            regions[lineRef[1]],
                                                            hotspot.catch(
                                                                [
                                                                    ".",
                                                                    ".",
                                                                    lineRef[5],
                                                                    lineRef[7],
                                                                ],
                                                                regions[lineRef[1]],
                                                                chromLengths,
                                                                genome,
                                                            ),
                                                            imdsData[lineRef[1]],
                                                            chromLengths[genome],
                                                            windowSize,
                                                        )
                                                    )
                                                ):
                                                    # if abs(float(linesSubClass[i][-1]) - float(lineRef[-1])) < vaf_cut and (int(linesSubClass[i][7])-int(lineRef[7]) <= imdsData[lineRef[1]] or (regions[lineRef[1]][hotspot.catch([".",".",lineRef[5], lineRef[7]], regions[lineRef[1]], chromLengths, genome, imds_corrected[lineRef[1]])] - (int(lineRef[7]) + chromLengths[genome][lineRef[5]]) < windowSize and int(lineRef[-2]) < imds_corrected[lineRef[1]][regions[lineRef[1]][hotspot.catch([".",".",lineRef[5], lineRef[7]], regions[lineRef[1]], chromLengths, genome, imds_corrected[lineRef[1]])]])):
                                                    # if abs(float(linesSubClass[i][-1]) - float(lineRef[-1])) < vaf_cut and (int(linesSubClass[i][7])-int(lineRef[7]) <= imdsData[lineRef[1]] or (regions[lineRef[1]][bisect.bisect_left(regions[lineRef[1]], (int(lineRef[7]) + chromLengths[genome][lineRef[5]]))] - (int(lineRef[7]) + chromLengths[genome][lineRef[5]]) < windowSize and int(lineRef[-2]) < imds_corrected[lineRef[1]][regions[lineRef[1]][bisect.bisect_left(regions[lineRef[1]], int(lineRef[7]) + chromLengths[genome][lineRef[5]])]])):
                                                    saveNewEvent.append(
                                                        linesSubClass[i]
                                                    )
                                                    lineRef = linesSubClass[i]
                                            else:
                                                if (
                                                    abs(
                                                        float(linesSubClass[i][-1])
                                                        - float(lineRef[-1])
                                                    )
                                                    < vaf_cut
                                                    and int(linesSubClass[i][7])
                                                    - int(lineRef[7])
                                                    <= imdsData[lineRef[1]]
                                                ):
                                                    saveNewEvent.append(
                                                        linesSubClass[i]
                                                    )
                                                    lineRef = linesSubClass[i]

                                    if len(saveNewEvent) > 1:
                                        if correction:
                                            # distancesLine = len([int(y[7])-int(x[7]) for x,y in zip(saveNewEvent, saveNewEvent[1:]) if int(y[7])-int(x[7]) > 1 and (int(y[7])-int(x[7]) <= imdsData[y[1]] or (regions[y[1]][bisect.bisect_left(regions[y[1]], (int(y[7]) + chromLengths[genome][y[5]]))] - (int(y[7]) + chromLengths[genome][y[5]]) < windowSize and int(y[-2]) < imds_corrected[y[1]][regions[y[1]][bisect.bisect_left(regions[y[1]], int(y[7]) + chromLengths[genome][y[5]])]]))])
                                            # distancesLine = len([int(y[7])-int(x[7]) for x,y in zip(saveNewEvent, saveNewEvent[1:]) if int(y[7])-int(x[7]) > 1 and (int(y[7])-int(x[7]) <= imdsData[y[1]] or (regions[y[1]][hotspot.catch([".",".",y[5], y[7]], regions[y[1]], chromLengths, genome, imds_corrected[y[1]])] - (int(y[7]) + chromLengths[genome][y[5]]) < windowSize and int(y[-2]) < imds_corrected[y[1]][regions[y[1]][hotspot.catch([".",".",y[5], y[7]], regions[y[1]], chromLengths, genome, imds_corrected[y[1]])]]))])
                                            distancesLine = len(
                                                [
                                                    int(y[7]) - int(x[7])
                                                    for x, y in zip(
                                                        saveNewEvent, saveNewEvent[1:]
                                                    )
                                                    if int(y[7]) - int(x[7]) > 1
                                                    and (
                                                        int(y[7]) - int(x[7])
                                                        <= imdsData[y[1]]
                                                        or (
                                                            len(regions[y[1]]) > 0
                                                            and regions[y[1]][
                                                                hotspot.catch(
                                                                    [
                                                                        ".",
                                                                        ".",
                                                                        y[5],
                                                                        y[7],
                                                                    ],
                                                                    regions[y[1]],
                                                                    chromLengths,
                                                                    genome,
                                                                )
                                                            ]
                                                            - (
                                                                int(y[7])
                                                                + chromLengths[genome][
                                                                    y[5]
                                                                ]
                                                            )
                                                            < windowSize
                                                            and hotspot.cutoffCatch(
                                                                [
                                                                    int(y[-2]),
                                                                    y[1],
                                                                    y[5],
                                                                    y[7],
                                                                ],
                                                                imds_corrected,
                                                                regions[y[1]],
                                                                hotspot.catch(
                                                                    [
                                                                        ".",
                                                                        ".",
                                                                        y[5],
                                                                        y[7],
                                                                    ],
                                                                    regions[y[1]],
                                                                    chromLengths,
                                                                    genome,
                                                                ),
                                                                imdsData[y[1]],
                                                                chromLengths[genome],
                                                                windowSize,
                                                            )
                                                        )
                                                    )
                                                ]
                                            )

                                        else:
                                            distancesLine = len(
                                                [
                                                    int(y[7]) - int(x[7])
                                                    for x, y in zip(
                                                        saveNewEvent, saveNewEvent[1:]
                                                    )
                                                    if int(y[7]) - int(x[7]) > 1
                                                    and (
                                                        int(y[7]) - int(x[7])
                                                        <= imdsData[y[1]]
                                                    )
                                                ]
                                            )
                                        vafs = len(
                                            [
                                                abs(float(y[-1]) - float(x[-1]))
                                                for x, y in zip(
                                                    saveNewEvent, saveNewEvent[1:]
                                                )
                                                if round(
                                                    abs(float(y[-1]) - float(x[-1])), 4
                                                )
                                                > vaf_cut
                                            ]
                                        )
                                        unknownVafs = len(
                                            [
                                                abs(float(y[-1]) - float(x[-1]))
                                                for x, y in zip(
                                                    saveNewEvent, saveNewEvent[1:]
                                                )
                                                if abs(float(y[-1]) - float(x[-1])) < 0
                                            ]
                                        )
                                        distancesActual = [
                                            int(y[7]) - int(x[7])
                                            for x, y in zip(
                                                saveNewEvent, saveNewEvent[1:]
                                            )
                                        ]
                                        vafsActual = [
                                            abs(float(y[-1]) - float(x[-1]))
                                            for x, y in zip(
                                                saveNewEvent, saveNewEvent[1:]
                                            )
                                        ]
                                        # Class I or II
                                        if distancesLine > 1:
                                            if len(saveNewEvent) >= 4:
                                                writeClassII = True
                                            else:
                                                writeClassI = True
                                                writeClassIc = True
                                        else:
                                            writeClassI = True
                                            if distancesLine == 0:
                                                if len(saveNewEvent) == 2:
                                                    writeClassIa = True
                                                else:
                                                    writeClassIb = True
                                            else:
                                                writeClassIc = True
                                    else:
                                        writeClassIII = True
                                        if unknownVafs > 0:
                                            category = "unknown"
                                        else:
                                            category = "vaf"
                                    for line in saveNewEvent:
                                        linesSubClass.remove(line)

                                    if writeClassII:
                                        processivitySubclassification(
                                            saveNewEvent, out2Y, out2K, out2S, out2N
                                        )
                                        for i in range(0, len(saveNewEvent), 1):
                                            saveNewEvent[i].append("ClassII")
                                            print(
                                                "\t".join([x for x in saveNewEvent[i]]),
                                                file=out4,
                                            )
                                            saveNewEvent[i] = [
                                                str(count)
                                            ] + saveNewEvent[i]
                                            print(
                                                "\t".join([x for x in saveNewEvent[i]]),
                                                file=out2,
                                            )
                                        count += 1
                                        print("\n\n", file=out2)

                                    else:
                                        if writeClassI:
                                            # Writes Class I (Single Events)
                                            try:
                                                for i in range(0, len(saveNewEvent), 1):
                                                    saveNewEvent[i].append("ClassI")
                                                    print(
                                                        "\t".join(
                                                            [x for x in saveNewEvent[i]]
                                                        ),
                                                        file=out3,
                                                    )
                                            except:
                                                print(saveNewEvent)

                                            if writeClassIc:
                                                # Writes Class Ic (extended MBSs)
                                                try:
                                                    for i in range(
                                                        0, len(saveNewEvent), 1
                                                    ):
                                                        saveNewEvent[i][-1] = "ClassIC"
                                                        print(
                                                            "\t".join(
                                                                [
                                                                    x
                                                                    for x in saveNewEvent[
                                                                        i
                                                                    ]
                                                                ]
                                                            ),
                                                            file=out8,
                                                        )
                                                except:
                                                    print(saveNewEvent)
                                            elif writeClassIa:
                                                # Writes Class Ia (DBSs)
                                                try:
                                                    for i in range(
                                                        0, len(saveNewEvent), 1
                                                    ):
                                                        saveNewEvent[i][-1] = "ClassIA"
                                                        print(
                                                            "\t".join(
                                                                [
                                                                    x
                                                                    for x in saveNewEvent[
                                                                        i
                                                                    ]
                                                                ]
                                                            ),
                                                            file=out6,
                                                        )
                                                except:
                                                    print(saveNewEvent)
                                            elif writeClassIb:
                                                # print("yes")
                                                # Writes Class 1b (MBSs)
                                                try:
                                                    for i in range(
                                                        0, len(saveNewEvent), 1
                                                    ):
                                                        saveNewEvent[i][-1] = "ClassIB"
                                                        # lines[i].append(category)
                                                        print(
                                                            "\t".join(
                                                                [
                                                                    x
                                                                    for x in saveNewEvent[
                                                                        i
                                                                    ]
                                                                ]
                                                            ),
                                                            file=out7,
                                                        )
                                                    # print("yes", saveNewEvent)
                                                except:
                                                    print(saveNewEvent)

                                        elif writeClassIII:
                                            # print("yes")
                                            try:
                                                for i in range(0, len(saveNewEvent), 1):
                                                    saveNewEvent[i].append("ClassIII")
                                                    saveNewEvent[i].append(category)
                                                    print(
                                                        "\t".join(
                                                            [x for x in saveNewEvent[i]]
                                                        ),
                                                        file=out5,
                                                    )
                                            except:
                                                print(saveNewEvent)
                                if len(linesSubClass) != 0:
                                    try:
                                        category = "vaf"
                                        for i in range(0, len(linesSubClass), 1):
                                            linesSubClass[i].append("ClassIII")
                                            linesSubClass[i].append(category)
                                            print(
                                                "\t".join(
                                                    [x for x in linesSubClass[i]]
                                                ),
                                                file=out5,
                                            )
                                    except:
                                        print(linesSubClass)

                    lines = []

        with open(project_path + "clusteredStats" + path_suffix + ".pickle", "wb") as f:
            pickle.dump(len_mnvs, f)

    # Parallelize and generate matrices for each subclass
    subclasses = [
        "class1",
        "class2",
        "class3",
        "class1a",
        "class1b",
        "class1c",
        "class2Y",
        "class2K",
        "class2S",
        "class2N",
    ]
    if processors > len(subclasses):
        max_seed = len(subclasses)
    else:
        max_seed = processors
    pool = mp.Pool(max_seed)

    subclasses_parallel = [[] for i in range(max_seed)]
    subclass_bin = 0
    for subclass in subclasses:
        if subclass_bin == max_seed:
            subclass_bin = 0
        subclasses_parallel[subclass_bin].append(subclass)
        subclass_bin += 1
    results = []
    for i in range(0, len(subclasses_parallel), 1):
        r = pool.apply_async(
            generateMatrices,
            args=(genome, log_out, project_path, path_suffix, subclasses_parallel[i]),
        )
        results.append(r)
    pool.close()
    pool.join()
    for r in results:
        r.wait()
        if not r.successful():
            # Raises an error when not successful
            r.get()

    # try:
    # 	print("Generating matrices for Class 1 mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class1", genome, project_path + 'subclasses'+ path_suffix + '/class1/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2 mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2", genome, project_path + 'subclasses'+ path_suffix +'/class2/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 3 mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class3", genome, project_path + 'subclasses'+ path_suffix +'/class3/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 1a mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class1a", genome, project_path + 'subclasses'+ path_suffix +'/class1a/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 1b mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class1b", genome, project_path + 'subclasses'+ path_suffix +'/class1b/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 1c mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class1c", genome, project_path + 'subclasses'+ path_suffix +'/class1c/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2Y mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2Y", genome, project_path + 'subclasses'+ path_suffix +'/class2Y/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2K mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2K", genome, project_path + 'subclasses'+ path_suffix +'/class2K/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2S mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2S", genome, project_path + 'subclasses'+ path_suffix +'/class2S/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2N mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2N", genome, project_path + 'subclasses'+ path_suffix +'/class2N/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass


def findClustersOfClusters_noVAF(
    project,
    chrom_based,
    project_parent_path,
    windowSize,
    chromLengths,
    regions,
    log_out,
    genome,
    processors,
    imds,
    correction=True,
):
    """
    Subclassifies the clustered mutations into different categories including DBS, MBS, omikli, and kataegis. This is only performed
    for simple single base substutions when no VAFs are present. Indels are not subclassified using this scheme.

    Parameters:
                            project	->	user provided project name (string)
                    chrom_based	->	option to generate IMDs per chromosome (boolean; default=False)
    project_parent_path	->	the directory for the given project (string)
                     windowSize	->	the window size used to calculate the mutation densities across the genome (integer; default=None)
               chromLengths	->	a dictionary of the cumulative chromosome lengths for a given reference genome (dictionary)
                            regions	->	a dictionary that contains all of the regions used for calculating corrected IMDs. If correction=False, then it returns an empty datastructure. (dictionary)
                             genome	->	the reference genome used for the given analysis (string)
                               imds	->	a dictionary of all of the corrected IMDs. If correction=False, then it returns an empty datastructure (dictionary)
                     correction	->	optional parameter to perform a genome-wide mutational density correction (boolean; default=False)

    Returns:
            None

    Outputs:
            Subclassification files for all clustered mutations.
    """

    ####################################################################################
    # Organize paths and file names
    ####################################################################################
    path_suffix = ""
    path_suffix2 = ""
    if chrom_based:
        path_suffix = "_chrom"
    if correction:
        path_suffix2 += "_corrected"

    project_path = (
        project_parent_path
        + "output/vcf_files"
        + path_suffix2
        + "/"
        + project
        + "_clustered/"
    )
    file = project_path + "/SNV/" + project + "_clustered.txt"
    out_file = project_path + project + "_clusters_of_clusters.txt"
    out_file2 = project_path + project + "_clusters_of_clusters_imd.txt"

    out_file3 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class1/"
        + project
        + "_clustered_class1.txt"
    )
    out_file4 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2/"
        + project
        + "_clustered_class2.txt"
    )
    out_file5 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class3/"
        + project
        + "_clustered_class3.txt"
    )
    out_file6 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class1a/"
        + project
        + "_clustered_class1a.txt"
    )
    out_file7 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class1b/"
        + project
        + "_clustered_class1b.txt"
    )
    out_file8 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class1c/"
        + project
        + "_clustered_class1c.txt"
    )
    out_file9 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2Y/"
        + project
        + "_clustered_class2Y.txt"
    )
    out_file10 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2K/"
        + project
        + "_clustered_class2K.txt"
    )
    out_file11 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2S/"
        + project
        + "_clustered_class2S.txt"
    )
    out_file12 = (
        project_path
        + "subclasses"
        + path_suffix
        + "/class2N/"
        + project
        + "_clustered_class2N.txt"
    )

    if os.path.exists(project_path + "subclasses" + path_suffix + "/class1/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class1/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class1/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class3/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class3/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class3/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class1a/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class1a/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class1a/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class1b/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class1b/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class1b/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class1c/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class1c/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class1c/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2Y/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2Y/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2Y/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2K/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2K/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2K/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2S/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2S/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2S/")
    if os.path.exists(project_path + "subclasses" + path_suffix + "/class2N/"):
        shutil.rmtree(project_path + "subclasses" + path_suffix + "/class2N/")
    os.makedirs(project_path + "subclasses" + path_suffix + "/class2N/")
    ####################################################################################

    # Hard-coded cutoffs, which can later be used as additional parameters to the tool.
    cutoff = 10001  # Min distance required between adjacenet events

    # Load the IMD data structures
    if chrom_based:
        with open(
            project_parent_path + "output/simulations/data/imds_chrom.pickle", "rb"
        ) as handle:
            imdsData = pickle.load(handle)
    else:
        with open(
            project_parent_path + "output/simulations/data/imds.pickle", "rb"
        ) as handle:
            imdsData = pickle.load(handle)
    if correction:
        with open(
            project_parent_path + "output/simulations/data/imds_corrected.pickle", "rb"
        ) as handle:
            imds_corrected = pickle.load(handle)

    first_pass = (
        True  # Deprecated option when originally developing. Can probably be deleted
    )
    total_muts = {}
    if first_pass:
        mnv_length = 0
        len_mnvs = {}
        total_mnvs = {}
        distances = []
        count = 1
        out = open(out_file, "w")
        with open(file) as f:
            next(f)
            lines = [line.strip().split() for line in f]

        for i in range(1, len(lines), 1):
            prev_chrom = lines[i - 1][5]
            prev_pos = int(lines[i - 1][6])
            prev_samp = lines[i - 1][1]
            chrom = lines[i][5]
            pos = int(lines[i][6])
            samp = lines[i][1]
            if samp not in total_muts:
                total_muts[samp] = 0
            if prev_samp == samp:
                if prev_chrom == chrom:
                    if pos - prev_pos < cutoff or (
                        correction
                        and len(regions[samp]) > 0
                        and (
                            (
                                regions[samp][
                                    hotspot.catch(
                                        [".", ".", chrom, pos],
                                        regions[samp],
                                        chromLengths,
                                        genome,
                                    )
                                ]
                                - (pos + chromLengths[genome][chrom])
                                < windowSize
                            )
                            and hotspot.cutoffCatch(
                                [(pos - prev_pos), samp, chrom, pos],
                                imds_corrected,
                                regions[samp],
                                hotspot.catch(
                                    [".", ".", chrom, pos],
                                    regions[samp],
                                    chromLengths,
                                    genome,
                                ),
                                imdsData[samp],
                                chromLengths[genome],
                                windowSize,
                            )
                        )
                    ):  # (pos - prev_pos) < imds_corrected[samp][regions[samp][hotspot.catch([".",".",chrom, pos], regions[samp], chromLengths, genome, imds_corrected[samp])]])):
                        distances.append(pos - prev_pos)
                        mnv_length += 1
                        lines[i - 1] = [str(count)] + lines[i - 1]
                        print("\t".join([x for x in lines[i - 1]]), file=out)
                        total_muts[samp] += 1
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
                        lines[i - 1] = [str(count)] + lines[i - 1]
                        print("\t".join([x for x in lines[i - 1]]), file=out)
                        total_muts[samp] += 1
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

                    lines[i - 1] = [str(count)] + lines[i - 1]
                    print("\t".join([x for x in lines[i - 1]]), file=out)
                    total_muts[samp] += 1
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

                lines[i - 1] = [str(count)] + lines[i - 1]
                print("\t".join([x for x in lines[i - 1]]), file=out)
                total_muts[samp] += 1
                count += 1
                count = 1
                print("\n\n################ New Sample #################", file=out)
                print("\n\n", file=out)

        lines[i] = [str(count)] + lines[i]
        print("\t".join([x for x in lines[i]]), file=out)
        total_muts[samp] += 1
        out.close()
        with open(
            project_parent_path
            + "output/simulations/data/originalCounts"
            + path_suffix2
            + ".pickle",
            "wb",
        ) as f:
            pickle.dump(total_muts, f)

    if True:
        len_mnvs = {"I": {}, "II": {}, "III": {}, "Ia": {}, "Ib": {}, "Ic": {}}
        distances = []
        distances_mnv = {}
        lines = []
        count = 1
        subclassesHeader = "\t".join(
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
                    "group",
                    "IMD",
                    "subclass",
                ]
            ]
        )
        with open(out_file) as f, open(out_file2, "w") as out2, open(
            out_file3, "w"
        ) as out3, open(out_file4, "w") as out4, open(out_file5, "w") as out5, open(
            out_file6, "w"
        ) as out6, open(
            out_file7, "w"
        ) as out7, open(
            out_file8, "w"
        ) as out8, open(
            out_file9, "w"
        ) as out2Y, open(
            out_file10, "w"
        ) as out2K, open(
            out_file11, "w"
        ) as out2S, open(
            out_file12, "w"
        ) as out2N:
            print(
                "\t".join(
                    [
                        x
                        for x in [
                            "clust_group",
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
                        ]
                    ]
                ),
                file=out2,
            )
            print(subclassesHeader, file=out3)
            print(subclassesHeader, file=out4)
            print(subclassesHeader + "\treason", file=out5)
            print(subclassesHeader, file=out6)
            print(subclassesHeader, file=out7)
            print(subclassesHeader, file=out8)
            print(subclassesHeader, file=out2Y)
            print(subclassesHeader, file=out2K)
            print(subclassesHeader, file=out2S)
            print(subclassesHeader, file=out2N)

            for line in f:
                # line = line.strip().split()[1:]
                line = line.strip().split()
                if line != []:
                    line = line[1:-2] + [line[0]] + line[-2:]
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
                            if correction:
                                distancesLine = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if int(y[7]) - int(x[7]) > 1
                                        and (
                                            int(y[7]) - int(x[7]) <= imdsData[y[1]]
                                            or (
                                                len(regions[y[1]]) > 0
                                                and regions[y[1]][
                                                    hotspot.catch(
                                                        [".", ".", y[5], y[7]],
                                                        regions[y[1]],
                                                        chromLengths,
                                                        genome,
                                                    )
                                                ]
                                                - (
                                                    int(y[7])
                                                    + chromLengths[genome][y[5]]
                                                )
                                                < windowSize
                                                and hotspot.cutoffCatch(
                                                    [int(y[-2]), y[1], y[5], y[7]],
                                                    imds_corrected,
                                                    regions[y[1]],
                                                    hotspot.catch(
                                                        [".", ".", y[5], y[7]],
                                                        regions[y[1]],
                                                        chromLengths,
                                                        genome,
                                                    ),
                                                    imdsData[y[1]],
                                                    chromLengths[genome],
                                                    windowSize,
                                                )
                                            )
                                        )
                                    ]
                                )  #   int(y[-2]) < imds_corrected[y[1]][regions[y[1]][hotspot.catch([".",".",y[5], y[7]], regions[y[1]], chromLengths, genome, imds_corrected[y[1]])]])   )])
                                distancesFailed = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if (
                                            int(y[7]) - int(x[7]) > imdsData[y[1]]
                                            and (
                                                len(regions[y[1]]) > 0
                                                and regions[y[1]][
                                                    hotspot.catch(
                                                        [".", ".", y[5], y[7]],
                                                        regions[y[1]],
                                                        chromLengths,
                                                        genome,
                                                    )
                                                ]
                                                - (
                                                    int(y[7])
                                                    + chromLengths[genome][y[5]]
                                                )
                                                < windowSize
                                                and hotspot.cutoffCatch(
                                                    [int(y[-2]), y[1], y[5], y[7]],
                                                    imds_corrected[y[1]],
                                                    regions[y[1]],
                                                    hotspot.catch(
                                                        [".", ".", y[5], y[7]],
                                                        regions[y[1]],
                                                        chromLengths,
                                                        genome,
                                                    ),
                                                    imdsData[y[1]],
                                                    chromLengths[genome],
                                                    windowSize,
                                                )
                                            )
                                        )
                                    ]
                                )
                                zeroDistances = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if int(y[7]) - int(x[7]) > 1
                                    ]
                                )

                            else:
                                distancesLine = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if int(y[7]) - int(x[7]) > 1
                                        and (int(y[7]) - int(x[7]) <= imdsData[y[1]])
                                    ]
                                )
                                distancesFailed = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if (int(y[7]) - int(x[7]) > imdsData[y[1]])
                                    ]
                                )
                                zeroDistances = len(
                                    [
                                        int(y[7]) - int(x[7])
                                        for x, y in zip(lines, lines[1:])
                                        if int(y[7]) - int(x[7]) > 1
                                    ]
                                )

                            if distancesFailed == 0:
                                # Class I or II
                                if distancesLine > 1:
                                    # Class IC or II
                                    if len(lines) >= 4:
                                        writeClassII = True
                                    else:
                                        writeClassI = True
                                        writeClassIc = True
                                else:
                                    writeClassI = True
                                    if zeroDistances == 0:
                                        if len(lines) == 2:
                                            writeClassIa = True
                                        else:
                                            writeClassIb = True
                                    else:
                                        writeClassIc = True
                            else:
                                writeClassIII = True

                        if writeClassII:
                            processivitySubclassification(
                                lines, out2Y, out2K, out2S, out2N
                            )
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
                                        print(
                                            "\t".join([x for x in lines[i]]), file=out3
                                        )
                                except:
                                    print(lines)

                                if writeClassIc:
                                    # Writes Class Ic (Omiklis)
                                    try:
                                        for i in range(0, len(lines), 1):
                                            lines[i][-1] = "ClassIC"
                                            print(
                                                "\t".join([x for x in lines[i]]),
                                                file=out8,
                                            )
                                    except:
                                        print(lines)
                                elif writeClassIa:
                                    # Writes Class Ia (DBSs)
                                    try:
                                        for i in range(0, len(lines), 1):
                                            lines[i][-1] = "ClassIA"
                                            print(
                                                "\t".join([x for x in lines[i]]),
                                                file=out6,
                                            )
                                    except:
                                        print(lines)
                                elif writeClassIb:
                                    # Writes Class 1b (MBSs)
                                    try:
                                        for i in range(0, len(lines), 1):
                                            lines[i][-1] = "ClassIB"
                                            print(
                                                "\t".join([x for x in lines[i]]),
                                                file=out7,
                                            )
                                    except:
                                        print(lines)

                    lines = []

            write_out = False
            category = None
            writeClassI = False
            writeClassII = False
            writeClassIII = False
            writeClassII = False
            writeClassIb = False
            writeClassIc = False
            writeClassIa = False

            if correction:
                distancesLine = len(
                    [
                        int(y[7]) - int(x[7])
                        for x, y in zip(lines, lines[1:])
                        if int(y[7]) - int(x[7]) > 1
                        and (
                            int(y[7]) - int(x[7]) <= imdsData[y[1]]
                            or (
                                len(regions[y[1]]) > 0
                                and regions[y[1]][
                                    hotspot.catch(
                                        [".", ".", y[5], y[7]],
                                        regions[y[1]],
                                        chromLengths,
                                        genome,
                                    )
                                ]
                                - (int(y[7]) + chromLengths[genome][y[5]])
                                < windowSize
                                and hotspot.cutoffCatch(
                                    [int(y[-2]), y[1], y[5], y[7]],
                                    imds_corrected,
                                    regions[y[1]],
                                    hotspot.catch(
                                        [".", ".", y[5], y[7]],
                                        regions[y[1]],
                                        chromLengths,
                                        genome,
                                    ),
                                    imdsData[y[1]],
                                    chromLengths[genome],
                                    windowSize,
                                )
                            )
                        )
                    ]
                )  #   int(y[-2]) < imds_corrected[y[1]][regions[y[1]][hotspot.catch([".",".",y[5], y[7]], regions[y[1]], chromLengths, genome, imds_corrected[y[1]])]])   )])
                distancesFailed = len(
                    [
                        int(y[7]) - int(x[7])
                        for x, y in zip(lines, lines[1:])
                        if (
                            int(y[7]) - int(x[7]) > imdsData[y[1]]
                            and (
                                len(regions[y[1]]) > 0
                                and regions[y[1]][
                                    hotspot.catch(
                                        [".", ".", y[5], y[7]],
                                        regions[y[1]],
                                        chromLengths,
                                        genome,
                                    )
                                ]
                                - (int(y[7]) + chromLengths[genome][y[5]])
                                < windowSize
                                and hotspot.cutoffCatch(
                                    [int(y[-2]), y[1], y[5], y[7]],
                                    imds_corrected[y[1]],
                                    regions[y[1]],
                                    hotspot.catch(
                                        [".", ".", y[5], y[7]],
                                        regions[y[1]],
                                        chromLengths,
                                        genome,
                                    ),
                                    imdsData[y[1]],
                                    chromLengths[genome],
                                    windowSize,
                                )
                            )
                        )
                    ]
                )
                zeroDistances = len(
                    [
                        int(y[7]) - int(x[7])
                        for x, y in zip(lines, lines[1:])
                        if int(y[7]) - int(x[7]) > 1
                    ]
                )

            else:
                distancesLine = len(
                    [
                        int(y[7]) - int(x[7])
                        for x, y in zip(lines, lines[1:])
                        if int(y[7]) - int(x[7]) > 1
                        and (int(y[7]) - int(x[7]) <= imdsData[y[1]])
                    ]
                )
                distancesFailed = len(
                    [
                        int(y[7]) - int(x[7])
                        for x, y in zip(lines, lines[1:])
                        if (int(y[7]) - int(x[7]) > imdsData[y[1]])
                    ]
                )
                zeroDistances = len(
                    [
                        int(y[7]) - int(x[7])
                        for x, y in zip(lines, lines[1:])
                        if int(y[7]) - int(x[7]) > 1
                    ]
                )

            if distancesFailed == 0:
                # Class I or II
                if distancesLine > 1:
                    # Class IC or II
                    if len(lines) >= 4:
                        writeClassII = True
                    else:
                        writeClassI = True
                        writeClassIc = True
                else:
                    writeClassI = True
                    if zeroDistances == 0:
                        if len(lines) == 2:
                            writeClassIa = True
                        else:
                            writeClassIb = True
                    else:
                        writeClassIc = True
            else:
                writeClassIII = True

            if writeClassII:
                processivitySubclassification(lines, out2Y, out2K, out2S, out2N)
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
                        # Writes Class Ic (Omiklis)
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
                        # Writes Class 1b (MBSs)
                        try:
                            for i in range(0, len(lines), 1):
                                lines[i][-1] = "ClassIB"
                                print("\t".join([x for x in lines[i]]), file=out7)
                        except:
                            print(lines)

        with open(project_path + "clusteredStats" + path_suffix + ".pickle", "wb") as f:
            pickle.dump(len_mnvs, f)

    # Parallelize and generate matrices for each subclass
    subclasses = [
        "class1",
        "class2",
        "class1a",
        "class1b",
        "class1c",
        "class2Y",
        "class2K",
        "class2S",
        "class2N",
    ]
    if processors > len(subclasses):
        max_seed = len(subclasses)
    else:
        max_seed = processors
    pool = mp.Pool(max_seed)

    subclasses_parallel = [[] for i in range(max_seed)]
    subclass_bin = 0
    for subclass in subclasses:
        if subclass_bin == max_seed:
            subclass_bin = 0
        subclasses_parallel[subclass_bin].append(subclass)
        subclass_bin += 1
    results = []
    for i in range(0, len(subclasses_parallel), 1):
        r = pool.apply_async(
            generateMatrices,
            args=(genome, log_out, project_path, path_suffix, subclasses_parallel[i]),
        )
        results.append(r)
    pool.close()
    pool.join()
    for r in results:
        r.wait()
        if not r.successful():
            # Raises an error when not successful
            r.get()

    # try:
    # 	print("Generating matrices for Class 1 mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class1", genome, project_path + 'subclasses'+ path_suffix + '/class1/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2 mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2", genome, project_path + 'subclasses'+ path_suffix +'/class2/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 1a mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class1a", genome, project_path + 'subclasses'+ path_suffix +'/class1a/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 1b mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class1b", genome, project_path + 'subclasses'+ path_suffix +'/class1b/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 1c mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class1c", genome, project_path + 'subclasses'+ path_suffix +'/class1c/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2Y mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2Y", genome, project_path + 'subclasses'+ path_suffix +'/class2Y/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2K mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2K", genome, project_path + 'subclasses'+ path_suffix +'/class2K/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2S mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2S", genome, project_path + 'subclasses'+ path_suffix +'/class2S/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass
    # try:
    # 	print("Generating matrices for Class 2N mutations:")
    # 	matrices = matGen.SigProfilerMatrixGeneratorFunc("class2N", genome, project_path + 'subclasses'+ path_suffix +'/class2N/', seqInfo=True)#, plot=True)
    # 	print()
    # except:
    # 	pass


def generateMatrices(genome, log_out, project_path, path_suffix, subclassesProcessor):
    temp = sys.stdout
    sys.stdout = open(log_out, "a")
    try:
        for subclass in subclassesProcessor:
            matrices = matGen.SigProfilerMatrixGeneratorFunc(
                subclass,
                genome,
                project_path + "subclasses" + path_suffix + "/" + subclass + "/",
                seqInfo=True,
            )  # , plot=True)
        sys.stdout.close()
    except:
        pass
