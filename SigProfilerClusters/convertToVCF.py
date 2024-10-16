import os
import shutil


def checkPaths(currentPath, suffix):
    inputPath = os.path.join(currentPath, suffix)
    if os.path.exists(inputPath):
        shutil.rmtree(inputPath)
    os.makedirs(inputPath)


def generateAllPaths(input_path, contexts):
    """
    Generates all necessary output paths to save the clustered output as VCF files.

    Parameters:
            input_path	->	the directory for the given project. This should contain the output from SigProfilerMatrixGenerator/Simulator (string)

    Returns:
            None

    Outputs:
            Clustered output paths found under : [project_path]/output/clustered/
            Non-clustered output paths found under: [project_path]/output/nonClustered/
    """

    # Clustered and non-clustered paths
    clusteredOutputPath = os.path.join(input_path, "output", "clustered")
    nonClusteredOutputPath = os.path.join(input_path, "output", "nonClustered")

    if not os.path.exists(clusteredOutputPath):
        os.makedirs(clusteredOutputPath)
    if not os.path.exists(nonClusteredOutputPath):
        os.makedirs(nonClusteredOutputPath)

    if contexts == "ID":
        checkPaths(clusteredOutputPath, "ID")
        checkPaths(nonClusteredOutputPath, "ID")
    else:
        checkPaths(clusteredOutputPath, "DBS")
        checkPaths(clusteredOutputPath, "MBS")
        checkPaths(clusteredOutputPath, "omikli")
        checkPaths(clusteredOutputPath, "kataegis")
        checkPaths(clusteredOutputPath, "other")
        checkPaths(nonClusteredOutputPath, "SBS")


def pullLineInfo(line, collectVAF=True):
    line = line.strip().split()
    sample = line[1]
    chrom = line[5]
    pos = line[6]
    ref = line[8]
    alt = line[9]

    if collectVAF:
        return (sample, chrom, pos, ref, alt, line[14], line[12])
    else:
        return (sample, chrom, pos, ref, alt)


def convertFiles(input_path, contexts, project, correction):
    subclasses = {
        "class1a": "DBS",
        "class1b": "MBS",
        "class1c": "omikli",
        "class2": "kataegis",
        "class3": "other",
    }
    if contexts == "ID":
        if correction == True:
            idPath = os.path.join(
                input_path,
                "output",
                "vcf_files_corrected",
                project + "_clustered",
                "INDEL",
                project + "_clustered.txt",
            )
        else:
            idPath = os.path.join(
                input_path,
                "output",
                "vcf_files",
                project + "_clustered",
                "INDEL",
                project + "_clustered.txt",
            )
        outputFile = os.path.join(input_path, "output", "clustered", "ID")
        with open(idPath) as f:
            next(f)
            allLines = f.readlines()
            if len(allLines) == 0:
                pass
            (
                currentSample,
                chrom,
                pos,
                ref,
                alt,
            ) = pullLineInfo(allLines[0], False)
            out = open(
                os.path.join(
                    input_path, "output", "clustered", "ID", currentSample + ".vcf"
                ),
                "a",
            )
            print(
                "\t".join(
                    ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                ),
                file=out,
            )
            print(
                "\t".join(
                    [chrom, pos, currentSample, ref, alt, ".", ".", "ClusteredAnalysis"]
                ),
                file=out,
            )
            for lines in allLines[1:]:
                (
                    sample,
                    chrom,
                    pos,
                    ref,
                    alt,
                ) = pullLineInfo(lines, False)
                if sample != currentSample:
                    out.close()
                    currentSample = sample
                    out = open(
                        os.path.join(
                            input_path,
                            "output",
                            "clustered",
                            "ID",
                            currentSample + ".vcf",
                        ),
                        "a",
                    )
                    print(
                        "\t".join(
                            [
                                "#CHROM",
                                "POS",
                                "ID",
                                "REF",
                                "ALT",
                                "QUAL",
                                "FILTER",
                                "INFO",
                            ]
                        ),
                        file=out,
                    )
                print(
                    "\t".join(
                        [chrom, pos, sample, ref, alt, ".", ".", "ClusteredAnalysis"]
                    ),
                    file=out,
                )
        out.close()

        if correction == True:
            nonClusteredFile = os.path.join(
                input_path,
                "output",
                "vcf_files_corrected",
                project + "_nonClustered",
                "INDEL",
                project + "_nonClustered.txt",
            )
        else:
            nonClusteredFile = os.path.join(
                input_path,
                "output",
                "vcf_files",
                project + "_nonClustered",
                "INDEL",
                project + "_nonClustered.txt",
            )
        with open(nonClusteredFile) as f:
            next(f)
            allLines = f.readlines()
            if len(allLines) == 0:
                pass
            currentSample, chrom, pos, ref, alt = pullLineInfo(allLines[0], False)
            out = open(
                os.path.join(
                    input_path, "output", "nonClustered", "ID", currentSample + ".vcf"
                ),
                "a",
            )
            print(
                "\t".join(
                    ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                ),
                file=out,
            )
            print(
                "\t".join(
                    [chrom, pos, currentSample, ref, alt, ".", ".", "NonClustered"]
                ),
                file=out,
            )
            for lines in allLines[1:]:
                sample, chrom, pos, ref, alt = pullLineInfo(lines, False)
                if sample != currentSample:
                    out.close()
                    currentSample = sample
                    out = open(
                        os.path.join(
                            input_path,
                            "output",
                            "nonClustered",
                            "ID",
                            currentSample + ".vcf",
                        ),
                        "a",
                    )
                    print(
                        "\t".join(
                            [
                                "#CHROM",
                                "POS",
                                "ID",
                                "REF",
                                "ALT",
                                "QUAL",
                                "FILTER",
                                "INFO",
                            ]
                        ),
                        file=out,
                    )
                print(
                    "\t".join([chrom, pos, sample, ref, alt, ".", ".", "NonClustered"]),
                    file=out,
                )
        out.close()

    else:
        if correction == True:
            subclassPath = os.path.join(
                input_path,
                "output",
                "vcf_files_corrected",
                project + "_clustered",
                "subclasses",
            )
        else:
            subclassPath = os.path.join(
                input_path, "output", "vcf_files", project + "_clustered", "subclasses"
            )
        for subclass in subclasses:
            subclassFile = os.path.join(
                subclassPath, subclass, project + "_clustered_" + subclass + ".txt"
            )
            outputFile = os.path.join(
                input_path, "output", "clustered", subclasses[subclass]
            )
            with open(subclassFile) as f:
                next(f)
                allLines = f.readlines()
                if len(allLines) == 0:
                    continue
                currentSample, chrom, pos, ref, alt, vaf, group = pullLineInfo(
                    allLines[0]
                )
                out = open(
                    os.path.join(
                        input_path,
                        "output",
                        "clustered",
                        subclasses[subclass],
                        currentSample + ".vcf",
                    ),
                    "a",
                )
                print(
                    "\t".join(
                        ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                    ),
                    file=out,
                )
                print(
                    "\t".join(
                        [
                            chrom,
                            pos,
                            currentSample,
                            ref,
                            alt,
                            ".",
                            ".",
                            "ClusteredAnalysis;VAF=" + vaf + ";groupNumber" + group,
                        ]
                    ),
                    file=out,
                )
                for lines in allLines[1:]:
                    sample, chrom, pos, ref, alt, vaf, group = pullLineInfo(lines)
                    if sample != currentSample:
                        out.close()
                        currentSample = sample
                        out = open(
                            os.path.join(
                                input_path,
                                "output",
                                "clustered",
                                subclasses[subclass],
                                currentSample + ".vcf",
                            ),
                            "a",
                        )
                        print(
                            "\t".join(
                                [
                                    "#CHROM",
                                    "POS",
                                    "ID",
                                    "REF",
                                    "ALT",
                                    "QUAL",
                                    "FILTER",
                                    "INFO",
                                ]
                            ),
                            file=out,
                        )
                    print(
                        "\t".join(
                            [
                                chrom,
                                pos,
                                sample,
                                ref,
                                alt,
                                ".",
                                ".",
                                "ClusteredAnalysis;VAF=" + vaf + ";groupNumber" + group,
                            ]
                        ),
                        file=out,
                    )
            out.close()

        if correction == True:
            nonClusteredFile = os.path.join(
                input_path,
                "output",
                "vcf_files_corrected",
                project + "_nonClustered",
                "SNV",
                project + "_nonClustered.txt",
            )
        else:
            nonClusteredFile = os.path.join(
                input_path,
                "output",
                "vcf_files",
                project + "_nonClustered",
                "SNV",
                project + "_nonClustered.txt",
            )
        with open(nonClusteredFile) as f:
            next(f)
            allLines = f.readlines()
            if len(allLines) == 0:
                pass
            currentSample, chrom, pos, ref, alt = pullLineInfo(allLines[0], False)
            out = open(
                os.path.join(
                    input_path, "output", "nonClustered", "SBS", currentSample + ".vcf"
                ),
                "a",
            )
            print(
                "\t".join(
                    ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                ),
                file=out,
            )
            print(
                "\t".join(
                    [chrom, pos, currentSample, ref, alt, ".", ".", "NonClustered"]
                ),
                file=out,
            )
            for lines in allLines[1:]:
                sample, chrom, pos, ref, alt = pullLineInfo(lines, False)
                if sample != currentSample:
                    out.close()
                    currentSample = sample
                    out = open(
                        os.path.join(
                            input_path,
                            "output",
                            "nonClustered",
                            "SBS",
                            currentSample + ".vcf",
                        ),
                        "a",
                    )
                    print(
                        "\t".join(
                            [
                                "#CHROM",
                                "POS",
                                "ID",
                                "REF",
                                "ALT",
                                "QUAL",
                                "FILTER",
                                "INFO",
                            ]
                        ),
                        file=out,
                    )
                print(
                    "\t".join([chrom, pos, sample, ref, alt, ".", ".", "NonClustered"]),
                    file=out,
                )
        out.close()
