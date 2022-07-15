# diCal2

`diCal2` is a software program for inferring parameters of demographic histories from whole genome sequence data.  To this end, it combines several demography-aware conditional sampling distributions in a composite likelihood framework. This composite likelihood is then used in an EM algorithm, to find the demographic parameters maximizing the likelihood. The program also implements a genetic algorithm on top of the EM to allow for a more thorough navigation of the parameter space.

diCal2 is an implementation of the method described in the paper:

- Matthias Steinrücken, John A. Kamm, Yun S. Song: [Inference of complex population histories using whole-genome sequences from multiple populations](https://doi.org/10.1073/pnas.1905060116), PNAS 116 (34), pp. 17115-17120 (2019).  [Preprint: [here](http://dx.doi.org/10.1101/026591)]

It has been applied for demographic inference in the paper:

- Maanasa Raghavan, Matthias Steinrücken, Kelley Harris, Stephan Schiffels, Simon Rasmussen, Michael DeGiorgio, Anders Albrechtsen, Cristina Valdiosera, María C. Ávila-Arcos, Anna-Sapfo Malaspinas, ..., Yun S Song, Rasmus Nielsen, Eske Willerslev. [Genomic evidence for the Pleistocene and recent population history of Native Americans](https://doi.org/10.1126/science.aab3884).  Science 349, aab3884 (2015).


### Table of Contents
* [Quick Start](#quick-start)
* [EXAMPLES](#EXAMPLES)
* [BUILD FROM SOURCE](#BUILD-FROM-SOURCE)
* [LICENSES](#LICENSES)
* [EXTERNAL LIBRARIES LICENSES](#EXTERNAL-LIBRARIES-LICENSES)
* [CONTACT](#CONTACT)
* [VERSION HISTORY](#VERSION-HISTORY)

### Quick Start

##### REQUIREMENTS:

`diCal2` requires JRE 1.8 or higher to run.  To compile the sourcecode, JDK 1.8 or higher is required.

##### DOWNLOAD:

Download the java-executable file [diCal2.jar](https://github.com/popgenmethods/diCal2/raw/main/diCal2.jar) in the top level directory of this repositiry.

##### USAGE:

For additional usage, options, and examples, please consult the [manual](https://github.com/popgenmethods/diCal2/blob/main/manual/manual.pdf).

To run `diCal2`, in the main directory, execute

```
java -jar diCal2.jar <arguments>
```

To print the usage and see the available command line arguments, execute the command

```
java -jar diCal2.jar --help
```

For some example calls, see the section 6. EXAMPLES.

Note: Depending on the command line arguments, the program might require a substantial amount of memory. If the program fails due to insufficient memory, you can provide the '-Xmx' flag to have the JVM use more memory. For example, the following command will allow the program to use 10 GB:

```
java -Xmx10g -jar diCal2.jar <arguments>
```

Further, if the program produces unexpected output, it can be helpful to enable assertions (more checking), using the command

```
java -ea -jar diCal2.jar <arguments>
```

### EXAMPLES:

Note that all population genetic parameters are UNSCALED and the population size is given in number of diploid individuals. Furthermore, the sequence data used in the following examples consists only of a small number of sites, for illustrative purposes. Depending on the demographic scenario, it is necessary to use hundreds of Mbp to get accurate estimates.

#### EXAMPLE 1

 The following example estimates demographic parameters in a scenario of exponential growth:

```
java -jar diCal2.jar --paramFile examples/fromReadme/test.param --demoFile examples/fromReadme/exp.demo --ratesFile examples/fromReadme/exp.rates --vcfFile examples/fromReadme/test.vcf --vcfFilterPassString '.' --vcfReferenceFile examples/fromReadme/test.fa --configFile examples/fromReadme/exp.config --metaStartFile examples/fromReadme/exp.rand --seed 541816302422 --lociPerHmmStep 3 --compositeLikelihood lol --metaNumIterations 2 --metaKeepBest 1 --metaNumPoints 3 --numberIterationsEM 2 --numberIterationsMstep 2 --disableCoordinateWiseMStep --intervalType logUniform --intervalParams '11,0.01,4' --bounds '1.00001,1000;0.01,0.06;0.01,0.23;0.02,2;0.5,4' --metaParallelEmSteps 2 --parallel 4
```

The file "examples/test.param" containing mutation and recombination parameters is specified via the argument "--paramFile". The demographic history and the exponential rates are given in the files "examples/exp.demo" and "examples/exp.rates" (parameters "--demoFile" and "--ratesFile").  The wildcards in these files ('?' followed by a number) indicate which parameters should be estimated (and their order).  Using the same wildcard twice to indicate an identical parameter is permitted.  The parameters "--vcfFile examples/test.vcf --vcfFilterPassString '.'" indicate that the sequence data to analyze is given in VCF format in the file "examples/test.vcf", and the string to required in the FILTER column for a site to be considered valid is a dot '.'. "--configFile" specifies the config file "examples/exp.config" that indicates the sequence length, number of alleles, and number of subpopulation, followed by a list that specifies the multiplicites of the sequences from the sequence file in the different subpopulations.  The length of this list has to equal the number of sequences provided. The file "examples/exp.rand" provided using "--metaStartFile" contains lines of tab-separated starting points.

"--seed" specifies the seed, and "--lociPerHmmStep 3" indicates that 3 loci should be grouped into a meta locus. "--compositeLikelihood lol" indicates that the leave-one-out composite likelihood should be used.  

"--metaNumIterations 2 --metaKeepBest 1 --metaNumPoints 3" are parameters for the genetic algorithm.  Here, 2 generations are supposed to be computed, where we choose 1 parent for each, the generation should consist of 3 parameter sets. "--numberIterationsEM 2 --numberIterationsMstep 2 --disableCoordinateWiseMStep --printEmPath" indicates that for each particle, 2 EM steps should be done, and 2 iterations in each M step. "--disableCoordinateWiseMStep" indicates that the multi-dimensional Nelder-Mead algorithm should be used in the M-Step, rather than the performing independent NM-steps in each coordinate, while holding the others fixed (the default). Further, information about the EM path should be printed. The parameters "--intervalType logUniform --intervalParams '11,0.01,4'" specify the 11 discretization intervals for the HMM should be uniform in a logarithmic space. The lowest point is 0.01 and the highest 4 (coalescent-scaled time).  "--bounds" specifies lower and upper bounds for all 5 parameters.  Finally, the program should use 4 threads ("--parallel 4") and 2 meta steps should be computed in parallel ("--metaParallelEmSteps 2").

Every line of the output that starts with a '#' denotes logging information. The relevant results are not preceded by a '#' and has eight tab-separated values on each line: The log-likelihood of current estimate, time in ms for one step, current estimates for five, and an id string. The id string is consists of three integers separated by '_': The first is the meta iteration, the second the EM step, and the third is the index of the particle in the current generation. To find the maximum, you should choose the one with the highest likelihood from the last meta iteration after the maximal number of EM steps has been performed.


#### EXAMPLE 2

The following example estimates demographic parameter in a scenario of isolation with migration:

```
java -jar diCal2.jar --paramFile examples/fromReadme/test.param --demoFile examples/fromReadme/IM.demo --vcfFile examples/fromReadme/test.vcf --vcfFilterPassString '.' --vcfReferenceFile examples/fromReadme/test.fa --configFile examples/fromReadme/IM.config --metaStartFile examples/fromReadme/IM.rand --seed 60643714832 --lociPerHmmStep 4 --compositeLikelihood pcl --metaNumIterations 2 --metaKeepBest 1 --metaNumPoints 3 --numberIterationsEM 2 --numberIterationsMstep 2 --intervalType logUniform --intervalParams '11,0.01,4' --bounds '0.01,0.32;0.05,1.0001;0.05,5;0.05,5;0.02,2;0.9,5;0.1,500'
```

The file "examples/test.param" file containing mutation and recombination parameters is specified via the argument "--paramFile". The demographic history is given in the file "examples/IM.demo" (specified using "--demoFile").  The wildcards in this file ('?' followed by a number) indicates which parameters should be estimated (and their order).  Using the same wildcard twice to indicate an identical parameter is permitted.  The parameter "--sequenceFile examples/test_ref.seq" indicates that the sequence data to analyze is given in the file "examples/test_ref.seq" in the diCal2 specific format.  "--configFile" specifies the config file "examples/IM.config" that indicates the sequence length, number of alleles, and number of subpopulation, followed by a list that specifies the multiplicites of the sequences from the sequence file in the different subpopulations. The length of this list has to equal the number of sequences provided. The file "examples/IM.rand" provided using "--metaStartFile" contains lines of tab-separated starting points.

"--seed" specifies the seed, and "--lociPerHmmStep 4" indicates that 4 loci should be grouped into a meta locus. "--compositeLikelihood pcl" indicates that the pairwise composite likelihood should be used. "--metaNumIterations 2 --metaKeepBest 1 --metaNumPoints 3" are parameters for the genetic algorithm.  Here, 2 generations are supposed to be computed, where we choose 1 parent for each, the generation should consist of 3 parameter sets.  "--numberIterationsEM 2 --numberIterationsMstep 2 --printEmPath" indicates that for each particle, 2 EM steps should be done, and 2 iterations in each M step. Further, information about the EM path should be printed. The parameters "--intervalType logUniform --intervalParams '11,0.01,4'" specify the 11 discretization intervals for the HMM should be uniform in a logarithmic space. The lowest point is 0.01 and the highest 4 (coalescent-scaled time).  "--bounds" specifies lower and upper bounds for all 7 parameters.

Every line of the output that starts with a '#' denotes logging information. The relevant results are not preceded by a '#' and has eight tab-separated values on each line: The log-likelihood of current estimate, time in ms for one step, current estimates for seven, and an id string. The id string is consists of three integers separated by '_': The first is the meta iteration, the second the EM step, and the third is the index of the particle in the current generation. To find the maximum, you should choose the one with the highest likelihood from the last meta iteration after the maximal number of EM steps.

#### EXAMPLE 4

The following example estimates demographic parameter in a scenario of isolation with three extant populations:

```
java -jar diCal2.jar --paramFile examples/fromReadme/test.param --demoFile examples/fromReadme/three.demo --vcfFile examples/fromReadme/test.vcf --vcfFilterPassString '.' --vcfReferenceFile examples/fromReadme/test.fa --configFile examples/fromReadme/three.config --seed 241438375231 --lociPerHmmStep 3 --compositeLikelihood lol --numberIterationsEM 5 --numberIterationsMstep 4 --startPoint '0.1,0.2,0.1,0.1,0.1,0.1,0.2' --intervalType logUniform --intervalParams '11,0.01,4' --bounds '0.005,0.2;0.05,1.0001;0.05,5;0.05,5;0.05,5;0.03,3;0.1,5'
```

The file "examples/test.param" file containing mutation and recombination parameters is specified via the argument "--paramFile". The demographic history is given in the file "examples/three.demo" (specified using "--demoFile").  The wildcards in this file ('?' followed by a number) indicates which parameters should be estimated (and their order).  Using the same wildcard twice to indicate an identical parameter is permitted. The parameter "--sequenceFile examples/test_sites.seq" indicates that the sequence data to analyze is given in the file "examples/test_sites.seq" in the another diCal2 specific format.  "--configFile" specifies the config file "examples/three.config" that indicates the sequence length, number of alleles, and number of subpopulation, followed by a list that specifies the multiplicites of the sequences from the sequence file in the different subpopulations.

"--seed" specifies the seed, and "--lociPerHmmStep 3" indicates that 3 loci should be grouped into a meta locus. "--compositeLikelihood lol" indicates that the leave-one-out composite likelihood should be used. "--numberIterationsEM 5 --numberIterationsMstep 4 --printEmPath --coordinateWiseMStep --nmFraction 0.2 --startPoint '0.1,0.2,0.1,0.1,0.1,0.1,0.2'" that 5 EM steps should be done, and 4 iterations in each M step. Further, information about the EM path should be printed. A start point for the EM is also specified. The parameters "--intervalType logUniform --intervalParams '11,0.01,4'" specify the 11 discretization intervals for the HMM should be uniform in a logarithmic space. The lowest point is 0.01 and the highest 4 (coalescent-scaled time).  "--bounds" specifies lower and upper bounds for all 7 parameters.

Every line of the output that starts with a '#' denotes logging information. The relevant results are not preceded by a '#' and has eight tab-separated values on each line: The log-likelihood of current estimate, time in ms for one step, current estimates for seven, and an id string. The id string is consists of three integers separated by '_': The first is the meta iteration, the second the EM step, and the third is the index of the particle in the current generation. To find the maximum, you should choose the one with the highest likelihood from the last meta iteration after the maximal number of EM steps.


### BUILD FROM SOURCE:

Download the file diCal2_2_0_2.tar.gz and unpack it.  The external jars are included in this release and can be found in the subdirectory diCal2_lib/.  If they should be missing, consult section 7. EXTERNAL LIBRARIES LICENSES for download information.

Compile the source code by executing the command (in the main directory):

```
javac -cp ./src:./diCal2_lib/arpack_combined_all.jar:./diCal2_lib/commons-math3-3.0.jar:./diCal2_lib/Jama-1.0.2.jar:./diCal2_lib/Jampack.jar:./diCal2_lib/JSAP-2.1.jar:./diCal2_lib/lapack_simple.jar:./diCal2_lib/trove-3.0.0.jar src/edu/berkeley/diCal2/maximum_likelihood/StructureEstimationEM.java 
```

To run the software, use the command:

```
java -cp ./src:./diCal2_lib/arpack_combined_all.jar:./diCal2_lib/commons-math3-3.0.jar:./diCal2_lib/Jama-1.0.2.jar:./diCal2_lib/Jampack.jar:./diCal2_lib/JSAP-2.1.jar:./diCal2_lib/lapack_simple.jar:./diCal2_lib/trove-3.0.0.jar edu.berkeley.diCal2.maximum_likelihood.StructureEstimationEM
```

For example, to see the usage information, run:

```
java -cp ./src:./diCal2_lib/arpack_combined_all.jar:./diCal2_lib/commons-math3-3.0.jar:./diCal2_lib/Jama-1.0.2.jar:./diCal2_lib/Jampack.jar:./diCal2_lib/JSAP-2.1.jar:./diCal2_lib/lapack_simple.jar:./diCal2_lib/trove-3.0.0.jar edu.berkeley.diCal2.maximum_likelihood.StructureEstimationEM --help
```

### LICENSES:

The source code is released under the GNU General Public License, version 3.  The full text of the license can be found in LICENSE_GPLv3.txt, which should have been included with this README.txt.

This software comes packaged with several external libraries. See Section 7. EXTERNAL LIBRARIES of the REAMDE for the respective license information and download links.

### EXTERNAL LIBRARIES LICENSES

The following external are included in this release of diCal2:


The following libraries (jar-files) are required. Put the jar-files (or symbolic links) into the subdirectory diCal2_lib/

- JSAP-2.1.jar
	Download from http://sourceforge.net/projects/jsap/files/jsap/2.1/ (alt: http://www.martiansoftware.com/jsap/).
	License: GNU Library or Lesser General Public License version 2.0 (LGPLv2) (http://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

- arpack_combined_all.jar
	Download from http://en.sourceforge.jp/projects/sfnet_f2j/releases/
	License: BSD License

- Jampack.jar
	Download Jampack.zip from ftp://math.nist.gov/pub/Jampack/Jampack/AboutJampack.html and unpack. This creates the subdirectory 'Jampack'. Execute the command 'jar cf Jampack.jar Jampack/*.class' to create the jar-file 'jampack.jar'
	License: Consult ftp://math.nist.gov/pub/Jampack/Jampack/AboutJampack.html

- Jama-1.0.2.jar
	Download from http://math.nist.gov/javanumerics/jama/
	License: public domain

- commons-math3-3.0.jar
	Download commons-math3-3.0-bin.tar.gz from http://archive.apache.org/dist/commons/math/binaries/ and unpack. commons-math3-3.0.jar can be found in the subdirectory 'commons-math3-3.0'.
	License: Apache v2.0 (http://www.apache.org/licenses/LICENSE-2.0)

- trove-3.0.0.jar
	Download trove-3.0.0.zip from https://sourceforge.net/projects/trove4j/files/trove/3.0.0/ and unpack (put it into a new subdirectory before unpacking). trove-3.0.0.jar can be found in the subdirectory 'lib'.
	License: GNU Library or Lesser General Public License version 2.0 (LGPLv2) (http://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

- lapack_simple.jar
	Download jlapack-0.8.tgz from http://www.netlib.org/java/f2j/ and unpack. lapack_simple.jar can be found in this archive.
	License: Custom License:
	Copyright � 2015 The University of Tennessee. All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: 
	- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
	- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer listed in this license in the documentation and/or other materials provided with the distribution. 
	- Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

	This software is provided by the copyright holders and contributors "as is" and any express or implied warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose are disclaimed. in no event shall the copyright owner or contributors be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.


### CONTACT:

Please contact steinrue@uchicago.edu with bugs, comments, or questions regarding the software.

### VERSION HISTORY:

2.0.0: Initial release.

2.0.1: Some bugfixes.

2.0.2: Some bugfixes. Commmand line interface changed to take less parameters and instead use more default values.

2.0.3: Now reads haploid vcf-files.

2.0.4: Now with manual & some bugfixes.

2.0.5: Better vcf-reading & memory management.
