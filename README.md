# SigProfilerHotSpots
Tool for analyzing the inter-mutational distances between SNV-SNV and INDEL-INDEL mutations. Tool separates mutations into clustered and non-clustered groups on a sample-dependent basis. 


**INTRODUCTION**

The purpose of this document is to provide a guide for using the SigProfilerHotSpots framework.
To be completed...

**PREREQUISITES**

The framework is written in PYTHON, and uses additional SigProfiler packages:

SigProfilerMatrixGenerator
SigProfilerSimulator

**QUICK START GUIDE**

This section will guide you through the minimum steps required to perform a hot-spot analysis:

Install the python package using pip:
                          pip install SigProfilerHotSpots
                          
Install your desired reference genome from the command line/terminal as follows (available reference genomes are: GRCh37, GRCh38, mm9, and mm10):
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37', rsync=False, bash=True)
This will install the human 37 assembly as a reference genome. You may install as many genomes as you wish. If you have a firewall on your server, you may need to install rsync and use the rsync=True parameter. Similarly, if you do not have bash, 
use bash=False.
Place your vcf files in your desired output folder. It is recommended that you name this folder based on your project's name

To be completed...

Output File Structure

To be completed...


**LOG FILES**

All errors and progress checkpoints are saved into sigProfilerMatrixGenerator_[project]_[genome].err and sigProfilerMatrixGenerator_[project]_[genome].out, respectively. For all errors, please email the error and progress log files to the primary contact under CONTACT INFORMATION.

CITATION

E.N. Bergstrom, M.N. Huang, U. Mahto, M. Barnes, M.R. Stratton, S.G. Rozen, and L.B. Alexandrov (2019) SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. https://www.biorxiv.org/content/10.1101/653097v1

COPYRIGHT

Copyright (c) 2019, Erik Bergstrom [Alexandrov Lab] All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

CONTACT INFORMATION

Please address any queries or bug reports to Erik Bergstrom at ebergstr@eng.ucsd.edu
