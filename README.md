[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/qpmzw/wiki/home/) [![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause) [![Build Status](https://app.travis-ci.com/AlexandrovLab/SigProfilerClusters.svg?token=Fyqun3zxxDr3YzDRaKKm&branch=master)](https://app.travis-ci.com/AlexandrovLab/SigProfilerClusters)


# SigProfilerClusters
Tool for analyzing the inter-mutational distances between SNV-SNV and INDEL-INDEL mutations. Tool separates mutations into clustered and non-clustered groups on a sample-dependent basis and subclassifies all SNVs into a set category of clustered event: i) DBS; ii) MBS; iii) omikli; and iv) kataegis. Indels are not subclassifed. 

This tool was previously under the project name of SigProfilerHotSpots, but has been renamed to SigProfilerClusters. For all instructions below, SigProfilerClusters may be interchanged with SigProfilerHotSpots if the older version of the tools is being used.

**INTRODUCTION**

The purpose of this document is to provide a guide for using the SigProfilerClusters framework. An extensive Wiki page detailing the usage of this tool can be found at https://osf.io/qpmzw/wiki/home/.


**PREREQUISITES**

The framework is written in PYTHON, and uses additional SigProfiler packages:

  * PYTHON          version 3.4 or newer
  * SigProfilerMatrixGenerator (https://github.com/AlexandrovLab/SigProfilerMatrixGenerator)
  * SigProfilerSimulator (https://github.com/AlexandrovLab/SigProfilerSimulator)

Please visit their respective GitHub pages for detailed installation and usage instructions.

**QUICK START GUIDE**

This section will guide you through the minimum steps required to perform clustered analysis:

1a. Install the python package using pip (current package):
                          pip install SigProfilerClusters

1b. Install the python package using pip (deprecated version):
                          pip install SigProfilerHotSpots
                          
                          
Install your desired reference genome from the command line/terminal as follows (available reference genomes are: GRCh37, GRCh38, mm9, and mm10):
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37', rsync=False, bash=True)
```
This will install the human 37 assembly as a reference genome. You may install as many genomes as you wish. If you have a firewall on your server, you may need to install rsync and use the rsync=True parameter. Similarly, if you do not have bash, 
use bash=False.

2. Place your vcf files in your desired output folder. It is recommended that you name this folder based on your project's name. Before you can analyze clustered mutations, you need to generate a background model for each of your samples. To do this, generate a minimum of 100 simulations for your project (see SigProfilerSimulator for a detailed list of parameters):
```
>>from SigProfilerSimulator import SigProfilerSimulator as sigSim
>>sigSim.SigProfilerSimulator(project, project_path, genome, contexts=["96"], simulations=100, chrom_based=True)
```
3. Now the original mutations can be partitioned into clustered and non-clustered sets using the required parameters below:
```
>> from SigProfilerClusters import SigProfilerClusters as hp
>> hp.analysis(project, genome, contexts, simContext, input_path)
```
See below for a detailed list of available parameters

4. The partitioned vcf files are placed under [project_path]/ouput/vcf_files/[project]_clustered/ and  [project_path]/ouput/vcf_files/[project]_nonClustered/. You can visualize the results by looking at the IMD plots available under [project_path]/ouput/simulations/[project]_simulations_[genome]_[context]_intradistance_plots/.

**AVAILABLE PARAMETERS**

	Required:
            project:			[string] Unique name for the given project
            genome:			[string] Reference genome to use. Must be installed using SigProfilerMatrixGenerator
            contexts:			[string] Mutation context for measuring IMD (e.g. "6", "96", "1536", etc,)
            simContext: 		[list of strings] Mutations context that was used for generating the background model (e.g ["6144"] or ["96"])
            input_path:			[string] Path to the given project
    
    	Optional:
            analysis:	 		[string] Desired analysis pipeline. By default output_type='all'. Other options include "subClassify" and "hotspot". 
            sortSims:			[boolean] Option to sort the simulated files if they have already been sorted. By default sortSims=True to ensure accurate results. The files must be sorted for accurate results. 
            interdistance:			[string] The mutation types to calculate IMDs between - Use only when performing analysis of indels (default='ID').
            calculateIMD:		[boolean] Parameter to calculate the IMDs. This will save time if you need to rerun the subclassification step only (default=True).
            max_cpu:			[integer] Change the number of allocated CPUs. By default all CPUs are used
            subClassify:		[boolean] Subclassify the clustered mutations. Requires that VAF scores are available in TCGA or Sanger format. By default subClassify=False 
            plotIMDfigure:	[boolean] Parameter that generates IMD and mutational spectra plots for each sample (default=True).
            plotRainfall		[boolean] Parameter that generates rainfall plots for each sample using the subclassification of clustered events (default=True).
            
            The following parameters are used if the subClassify argument is True:
            includedVAFs:	[boolean] Parameter that informs the tool of the inclusion of VAFs in the dataset (default=True)
            sanger:			[boolean] The input files are from Sanger. By default sanger=True
            TCGA:			[boolean] The input files are from TCGA. By default TCGA=False
            windowSize:		[integer] Window size for calculating mutation density in the rainfall plots. By default windowSize=10000000
            correction		[boolean] Optional parameter to perform a genome-wide mutational density correction (boolean; default=False)


**LOG FILES**

All errors and progress checkpoints are saved into SigProfilerClusters_[project]_[genome].err and SigProfilerClusters_[project]_[genome].out, respectively. For all errors, please email the error and progress log files to the primary contact under CONTACT INFORMATION.

CITATIONS

Erik N Bergstrom, Mousumy Kundu, Noura Tbeileh, Ludmil B Alexandrov. bioRxiv 2022.02.11.480117; doi: https://doi.org/10.1101/2022.02.11.480117

COPYRIGHT

Copyright (c) 2022, Erik Bergstrom [Alexandrov Lab] All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

CONTACT INFORMATION

Please address any queries or bug reports to Erik Bergstrom at ebergstr@eng.ucsd.edu
