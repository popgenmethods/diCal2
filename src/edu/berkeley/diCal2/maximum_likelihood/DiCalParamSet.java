 /*
    This file is part of diCal2.

    diCal2 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    diCal2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with diCal2.  If not, see <http://www.gnu.org/licenses/>.
  */

package edu.berkeley.diCal2.maximum_likelihood;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.SimpleJSAP;
import com.martiansoftware.jsap.Switch;

import Jama.Matrix;
import edu.berkeley.diCal2.csd.DemoConfiguration;
import edu.berkeley.diCal2.csd.DemoConfiguration.ConfigInfo;
import edu.berkeley.diCal2.csd.DemoState;
import edu.berkeley.diCal2.csd.DemoState.DemoStateFactory;
import edu.berkeley.diCal2.csd.EigenParamSet;
import edu.berkeley.diCal2.csd.HmmStepHandler;
import edu.berkeley.diCal2.csd.HmmStepHandler.SimpleLocusSkipper;
import edu.berkeley.diCal2.csd.HmmStepHandler.SimpleXLocusSkipper;
import edu.berkeley.diCal2.csd.IntervalFactory;
import edu.berkeley.diCal2.csd.IntervalFactory.CustomFixedIntervalFactory;
import edu.berkeley.diCal2.csd.IntervalFactory.LogUniformFactory;
import edu.berkeley.diCal2.csd.TrunkProcess;
import edu.berkeley.diCal2.csd.TrunkProcess.CakeStyle;
import edu.berkeley.diCal2.csd.TrunkProcess.TrunkProcessFactory;
import edu.berkeley.diCal2.csd.UberDemographyCore;
import edu.berkeley.diCal2.demography.DemographyFactory;
import edu.berkeley.diCal2.haplotype.FSAParamSet;
import edu.berkeley.diCal2.haplotype.FSAParamSetReader;
import edu.berkeley.diCal2.haplotype.GeneticType;
import edu.berkeley.diCal2.haplotype.HFSAXFullHaplotype;
import edu.berkeley.diCal2.haplotype.HFSAXFullHaplotype.HFSAXFullHaplotypeShell;
import edu.berkeley.diCal2.haplotype.ReadSequences;
import edu.berkeley.diCal2.maximum_likelihood.StructureEstimationEM.CompositeLikelihoodType;
import edu.berkeley.diCal2.maximum_likelihood.StructureEstimationEM.ConditionalObjectiveFunctionType;
import edu.berkeley.diCal2.utility.StationaryDistribution;

public class DiCalParamSet {

	public static class ExtendedConfigInfo {
		public final DemoConfiguration<HFSAXFullHaplotypeShell>  structuredConfig;
		public final HmmStepHandler<HFSAXFullHaplotypeShell> fancyTransitionMap;
		public final List<DemeHapPair<HFSAXFullHaplotypeShell>> demeHaplotypeList;
		public final EigenParamSet pSet;
		public final TreeSet<Integer> printLoci;

		public final ConfigInfo configInfo;
		public final HFSAXFullHaplotypeShell additionalHap;
		public final List<List<Integer>> csdList;
		
		public ExtendedConfigInfo (List<HFSAXFullHaplotype> preSequenceList, ConfigInfo preConfigInfo, int[] numPerDeme, double theta, double rho, Matrix mutMatrix, boolean useLocusSkipping, int lociPerHmmStep, TreeSet<Integer> printLoci, boolean useStationaryForPartially, Integer additionalHapIdx, List<List<Integer>> csdList) throws IOException {
			assert (preConfigInfo == null || numPerDeme == null);
			assert (preConfigInfo != null || numPerDeme != null);

			HFSAXFullHaplotype preAdditionalHap = null;
			List<HFSAXFullHaplotypeShell> sequenceList = null;
			// just remember the csd list
			this.csdList = csdList;
			
			// read in the sequences correctly
			if (preConfigInfo != null) {
				// sieve sequences and config too (if there is one)
				assert (preConfigInfo.multiplicities.size() == preSequenceList.size());

				List<int[]> mult = new ArrayList<int[]>();

				// prepare the sequencelist to add the real sequences
				sequenceList = new ArrayList<HFSAXFullHaplotypeShell>();
				
				// go through and sieve
				for (int i=0; i<preSequenceList.size(); i++) {
					
					// see how many we have
					int sum = 0;
					int[] thisMult = preConfigInfo.multiplicities.get(i);
					// no more debugging
//					System.out.println(Arrays.toString(thisMult));
//					System.out.println(preSequenceList.get(i));
					for (int j : thisMult) {
						assert (j >= 0);
						sum += j;
					}
					
					// add mult if there is at least one
					if (sum > 0) {
						mult.add(thisMult);
					}
					
					// add the sequence if there is a valid one
					// add make sure that the multiplicities are there
					if (preSequenceList.get(i) != null) {
						if (sum > 0) {
							sequenceList.add (new HFSAXFullHaplotypeShell(preSequenceList.get(i)));
						}
						else {
							assert (i == additionalHapIdx);
							// we fake it for the moment (note, here we can be sure to have mutliplicities)
							assert (preAdditionalHap == null);
							preAdditionalHap = preSequenceList.get(i);
						}
					}
					else {
						// don't add it
						assert (sum == 0);
					}
				}
				
				// finally the config
				this.configInfo = new ConfigInfo(mult, preConfigInfo.numDemes, preConfigInfo.numLoci, preConfigInfo.numAlleles);
			}
			else {
				// no config file given, only a numPerDeme list
				assert (numPerDeme != null);
				// so just use the sequence list as given
				sequenceList = new ArrayList<HFSAXFullHaplotypeShell>();
				for (HFSAXFullHaplotype hap : preSequenceList) sequenceList.add (new HFSAXFullHaplotypeShell (hap));
				// create some fake mults
				List<int[]> mult = new ArrayList<int[]>();
				for (int deme = 0; deme < numPerDeme.length; deme++) {
					for (int i = 0; i < numPerDeme[deme]; i++) {
						// fake it
						int[] thisDeme = new int[numPerDeme.length];
						thisDeme[deme] = 1;
						mult.add(thisDeme);
					}
				}
				assert (mult.size() == sequenceList.size());
				// get the others from somwhere
				int numDemes = numPerDeme.length;
				int numLoci = sequenceList.get(0).getNumLoci();
				assert (mutMatrix.getColumnDimension() == mutMatrix.getRowDimension());
				int numAlleles = mutMatrix.getColumnDimension();
				// no config, so create a fake config info
				this.configInfo = new ConfigInfo (mult, numDemes, numLoci, numAlleles);
			}

			// if pre hap is null, this works too
			if (preAdditionalHap == null) {
				this.additionalHap = null;
			}
			else {
				this.additionalHap = new HFSAXFullHaplotypeShell (preAdditionalHap);
			}
			
			// some tests
			if (this.configInfo != null) {
				// make sure it agrees
				if (this.configInfo.multiplicities.size() != sequenceList.size()) {
					throw new IOException("Number of lines in configfile don't agree with number of haplotpyes given.");
				}
				assert (this.configInfo.numLoci == sequenceList.get(0).getNumLoci());
			}
			else {
				assert (this.additionalHap == null);
			}
			
			
			// finally make the list and the config
			
			// make the ordered list first
			this.demeHaplotypeList = new ArrayList<DemeHapPair<HFSAXFullHaplotypeShell>>();

			// go through all the multiplicites
			for (int m=0; m < this.configInfo.multiplicities.size(); m++) {
				// go through all demes
				for (int g=0; g<this.configInfo.multiplicities.get(m).length; g++) {
					// now add as many as there are
					// and with deme information
					for (int i=0; i<this.configInfo.multiplicities.get(m)[g]; i++) {
						this.demeHaplotypeList.add(new DemeHapPair<HFSAXFullHaplotypeShell>(g, sequenceList.get(m)));
					}
				}
			}

			// and then make the unordered config
			this.structuredConfig = new DemoConfiguration<HFSAXFullHaplotypeShell> (this.configInfo.numDemes, this.configInfo.numLoci, this.configInfo.numAlleles, this.demeHaplotypeList);
			
			this.pSet = new EigenParamSet(structuredConfig.numLoci(), this.configInfo.numAlleles, theta, mutMatrix, rho);			
			
			// get time
			long oldTime = System.currentTimeMillis();
			
			// create an appropriate step handler
			if (useLocusSkipping){
				assert (lociPerHmmStep == 1);
				if (useStationaryForPartially) throw new IOException ("Partially Missing only allowed with multilocus hmm");
				// we need it for the segregating sites
				List<HFSAXFullHaplotypeShell> validSequences = sequenceList;
				if (this.additionalHap != null) {
					// don't want to modify the list of sequences
					validSequences = new ArrayList<HFSAXFullHaplotypeShell>(sequenceList);
					validSequences.add(this.additionalHap);
				}
				this.fancyTransitionMap = new HmmStepHandler.SimpleXLocusSkipper (validSequences, pSet, printLoci);
			} else {
				if (lociPerHmmStep == 1) {
					if (useStationaryForPartially) throw new IOException ("Partially Missing only allowed with multilocus hmm");
					this.fancyTransitionMap = new HmmStepHandler.SingleStepLocusTransitionMap<HFSAXFullHaplotypeShell> (pSet);			
				} else {
					// need some copies
					List<HFSAXFullHaplotypeShell> trunkSequences = sequenceList;
					List<HFSAXFullHaplotypeShell> additionalSequences = null;
					
					// a single additional haplotype would have to go in additional
					if (this.additionalHap != null) {
						additionalSequences = new ArrayList<HFSAXFullHaplotypeShell>();
						additionalSequences.add(this.additionalHap);
					}
					// if we read from a file
					if (this.csdList != null) {
						// we need to partition as well
						// some sets
						Set<Integer> trunkSet = new TreeSet<Integer>();
						Set<Integer> additionalSet = new TreeSet<Integer>();
						
						// fill them
						for (List<Integer> csd : this.csdList) {
							for (int i=0; i<csd.size(); i++) {
								if (i == 0) {
									// first one always additional
									additionalSet.add (csd.get(i));
								}
								else {
									// rest is trunk
									trunkSet.add (csd.get(i));
								}
							}
						}
						
						// them come up with some lists
						trunkSequences = new ArrayList<HFSAXFullHaplotypeShell>();
						for (Integer trunkIdx : trunkSet) {
							trunkSequences.add (sequenceList.get(trunkIdx));
						}
						additionalSequences = new ArrayList<HFSAXFullHaplotypeShell>();
						for (Integer addIdx : additionalSet) {
							additionalSequences.add (sequenceList.get(addIdx));
						}
					}
					
					// now see about the step handler
					if (!useStationaryForPartially) {
						this.fancyTransitionMap = new HmmStepHandler.MultiLocusStepHandler<HFSAXFullHaplotypeShell>(trunkSequences, additionalSequences, pSet, lociPerHmmStep);
					}
					else {
						this.fancyTransitionMap = new HmmStepHandler.MultiLocusStepHandlerSuboptimalMissing<HFSAXFullHaplotypeShell>(trunkSequences, additionalSequences, pSet, lociPerHmmStep);
					}
					
					// after we constructed the multiLocus StepHandler, we can clear the old haplotypes to free some memory
					clearHaplotypes  (trunkSequences, additionalSequences);
					// suggest cleaning to the garbage collector
					System.gc();
				}
			}
			
			// timing
			long currTime = System.currentTimeMillis();
			long readTime = currTime - oldTime;
			System.out.println ("# " + readTime + " milliseconds needed to build the step handler.");

			
			// some more things
			this.printLoci = printLoci;
		}

		private void clearHaplotypes (List<HFSAXFullHaplotypeShell> trunkSequences, List<HFSAXFullHaplotypeShell> additionalSequences) {
			// go trough trunk sequences and clear them
			for (HFSAXFullHaplotypeShell hap : trunkSequences) {
				hap.clear();
			}
			// same with additional ones
			if (additionalSequences != null) {
				for (HFSAXFullHaplotypeShell hap : additionalSequences) {
					hap.clear();
				}
			}
		}
	}
	
	
	public final List<ExtendedConfigInfo> extendedConfigInfoList;
	
	public final DemographyFactory demoFactory;
	public final IntervalFactory intervalFactory;
	public final TrunkProcessFactory trunkFactory;
	public final DemoStateFactory demoStateFactory;
	public final JSAPResult jsapParams;

	public List<List<Integer>> loadedCsdList;

	public final Double nelderMeadFraction;

	public ConditionalObjectiveFunctionType conditionalObjectiveFunction;

	public final String compositeLikelihoodStringLow;
	public final CompositeLikelihoodType compositeLikelihood;

	private final String filterPassString;

	final public boolean useEigenCore;

	public static class DemeHapPair<H extends GeneticType> {
		public final Integer demeNumber;
		public final H haplotype;

		public DemeHapPair (Integer demeNumber, H haplotype) {
			this.demeNumber = demeNumber;
			this.haplotype = haplotype;
		}
	}
	
	public DiCalParamSet(String args[], Parameter[] additionalParams, Parameter[] additionalHiddenParams, PrintStream outStream) throws JSAPException, IOException {
		this(args, additionalParams, additionalHiddenParams, outStream, false);
	}
	
	public DiCalParamSet(String args[], Parameter[] additionalParams, Parameter[] additionalHiddenParams, PrintStream outStream, boolean requireOrderedConfig) throws JSAPException, IOException {
		
		Parameter[] defaultParams = new Parameter[] {
                
	            new FlaggedOption( "configFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'c', "configFile", 
	                "The input configuration file. Tells you meta parameter and which sequences go into which demes." ),
	            new FlaggedOption( "numPerDeme", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'n', "numPerDeme",
	            		"The number of individuals per deme. e.g. 5,3,2 means the first 5 haplotypes are in deme 0, the next 3 are in deme 1, and the next 2 are in deme 2."),
	                
	            new FlaggedOption( "vcfFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "vcfFile", 
	                "The sequence file. Just a bunch of sequences. In vcf format." ),
	            
		            new FlaggedOption("compositeLikelihood", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "compositeLikelihood",
	    			"The composite likelihood type (LOL, PAC, SuperPAC, PCL, PCLOL, or FILE)."),
	            
	            new FlaggedOption( "paramFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 'p', "paramFile", 
	                "The input parameter file." ),
	            // this is specified by the given mutation matrix anyway
//	            new Switch ("vcfDiAllelic", JSAP.NO_SHORTFLAG, "vcfDiAllelic",
//		        		"If this switch is set, the vcf file is represented with two alleles internally. Sites with more then two alleles are ignored."),
//		        new Switch ("vcfTriAllelic", JSAP.NO_SHORTFLAG, "vcfTriAllelic",
//		        		"If this switch is set, the vcf file is represented with three alleles internally. Sites with more then three alleles are ignored."),

	    		new FlaggedOption( "vcfFilterPassString", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "vcfFilterPassString", 
						"The String in the FILTER column of a vcf-file that marks a snp as having passed the filters. ALL OTHER SITES ARE IGNORED!"),

				new FlaggedOption( "vcfOffset", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "vcfOffset", 
						"Offset(s) to shift the positions given in the vcf file(s)."),
		
	    		new FlaggedOption( "vcfReferenceFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "vcfReferenceFile", 
						"Provide a file containing the reference sequence to be used with the vcf-file. Needs to be in fasta-format. Specifying this option overrides whatever is provided in the vcf-file."),
//	            new Switch ("vcfFixedReference", JSAP.NO_SHORTFLAG, "useExternalReference",
//                		"If this switch is set, all non-segregating sites are assumed to have the same arbitrary (fixed) allele. Specifying this overrides whatever is provided in the vcf-file."),
//	            new Switch ("vcfRandomReference", JSAP.NO_SHORTFLAG, "vcfRandomReference",
//                		"If this switch is set, all non-segregating sites are assumed to have the same random (uniformly distributed, fixed) allele. Specifying this overrides whatever is provided in the vcf-file."),
			
	            new FlaggedOption( "demoFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.REQUIRED, JSAP.NO_SHORTFLAG, "demoFile",
	            		"The demography parameter file."),
	            new FlaggedOption ("bounds", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "bounds", "A string of semi-colon separated list of pairs of doubles, each pair spearated by a comma. This list gives the bounds for the parameters during the EM-estimation."),
	                
	            new FlaggedOption( "ratesFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "ratesFile",
	               		"The file with the exponential growth rates. Has to kind of match the demogrpahy file."),

	            new FlaggedOption ("lociPerHmmStep", JSAP.INTEGER_PARSER, "1", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "lociPerHmmStep",
	            		"If this flag is used with k > 1, then each HMM-step consists of k loci, and we assume there are no recombinations within blocks of k loci. Not compatible with locus skipping."),
	            			                		
	            new FlaggedOption("intervalType", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "intervalType",
	        			"The method for creating the intervals. [single,simple,meanrate,mixtureuniform,old]"),
	        			
	        	new FlaggedOption ("intervalParams", JSAP.STRING_PARSER, "", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "intervalParams",
	        			"The parameters for generating the intervals. [single=<none>, simple=<none>, meanrate=<numIntervals>, mixtureuniform=<numIntervals>, old=<numIntervals>]"),
	        			
	            new Switch ("printIntervals", JSAP.NO_SHORTFLAG, "printIntervals",
	                    		"If this switch is set, the intervals used are printed."),

	            new Switch ("hidden", JSAP.NO_SHORTFLAG, "hidden",
                		"Show all command line parameter (including hidden ones). WARNING: Experts only."),
	        			
	            new FlaggedOption ("bedFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "bedFile",
	            		"Include *.bed file(s) that lists all regions of the chromosome that should be EXCLUDED from analysis")
	                    
	        };
		
		
		
		Parameter[] defaultHiddenParams = new Parameter[] {
	            new FlaggedOption( "sequenceFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, 's', "sequenceFile", 
		                "The sequence file. Just a bunch of sequences." ),
	            new Switch("unCompressedSequence", JSAP.NO_SHORTFLAG, "unCompressedSequence", "If this switch is set, the sequence file is uncompressed and contains full haplotypes."),
	            new Switch("stationaryForPartiallyMissing", JSAP.NO_SHORTFLAG, "stationaryForPartiallyMissing", "If this switch is set, the partially missing sites are not ignored, but the stationary distribution is used for them. Kinda only works with the multilocus hmm version."),
	            new Switch ("useRawMutationMatrix", JSAP.NO_SHORTFLAG, "useRawMutationMatrix",
                		"If this switch is set, the mutation matrix is used as is. Default mode is to renormalize the mutation matrix so rows sum to 1."),
    		new Switch ("makeUniformMutationMatrix", JSAP.NO_SHORTFLAG, "makeUniformMutationMatrix",
            		"If this switch is set, we use a uniform mutation matrix instead of the given one, and adjust the mutation rate accordingly."),   
            
    		new FlaggedOption( "alleles", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "alleles", 
    				"A comma separated string of alleles, e.g. 'A,C,G,T'."),
			new FlaggedOption( "missingAlleles", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "missingAlleles", 
    				"A comma separated string of the alleles that are interpreted as missing, e.g. 'N,T,S'."),
            new FlaggedOption( "minimalMigration", JSAP.DOUBLE_PARSER, "0.005", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "minimalMigration",
                    "If this factor is given, you can always migrate anywhere with a minial rate of min(theta,rho)*<factor>. (Except the rates you estimate)"),
            new Switch ("useLocusSkipping", JSAP.NO_SHORTFLAG, "useLocusSkipping",
            		"If this switch is set, we use the locus-skipping speedup (requires uniform PIM mutations)"),
			new FlaggedOption("trunkStyle", JSAP.STRING_PARSER, "migratingEthan", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "trunkStyle",
        			"The style of trunk. [simple,oldCake,meanCake,multiCake,multiCakeUpdating,migratingMultiCake,migMultiCakeUpdating,recursive,exactCake,migratingEthan]"),
        	  
        	new FlaggedOption("cakeStyle", JSAP.STRING_PARSER, "average", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "cakeStyle",
        		"Where in the interval to get the cake factor [beginning,middle,end,average]"),		
        	
        	new FlaggedOption("printLoci", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "printLoci",
        			"If this flag is specified with an integer k, additionally decodes at every k loci, so that we can print at them."),
        	new Switch ("printSnps", JSAP.NO_SHORTFLAG, "printSnps", "Print at the SNPs as well."),	
        		
        	new Switch ("acceptUnphasedAsMissing", JSAP.NO_SHORTFLAG, "acceptUnphasedAsMissing", "If this switch is set, we do accept unphased genotypes as missing alleles."),	

        	new Switch ("vcfIgnoreDoubleEntries", JSAP.NO_SHORTFLAG, "vcfIgnoreDoubleEntries", "If this switch is set, multiple entries in a vcf for same sequence position are igored, and only the first is used."),
        	
        	new FlaggedOption ("addTrunkIntervals", JSAP.INTEGER_PARSER, "0", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "addTrunkIntervals",
              		"Number of trunk intervals to add."),
              		
            new Switch ("oldCakeRenormalizeSettings", JSAP.NO_SHORTFLAG, "oldCakeRenormalizeSettings", "Use old cake intervals and don't renormalize"),
            
            new Switch ("oldCompressedDataFormat", JSAP.NO_SHORTFLAG, "oldCompressedDataFormat", 
            		"In old data format, the compressed sequence file didn't know how many total loci there were, and this info was gotten from the config file."),

            new FlaggedOption ("nmFraction", JSAP.DOUBLE_PARSER, "0.2", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "nmFraction",
                    "The initial stepsize of the nelder mead optimizer is this fraction of the starting point."),
		new Switch ("condOnTransitionType", JSAP.NO_SHORTFLAG, "condOnTransitionType", "If this switch is set, we condition on the type of transition (recombination or no recombination)."),
		new Switch ("marginalKL", JSAP.NO_SHORTFLAG, "marginalKL", "If this switch is set, instead of using EM objective function, we use KL divergence to average marginal posterior decoding."),
    	new Switch("ancientDemeStates", JSAP.NO_SHORTFLAG, "ancientDemeStates",
    			"If this switch is set, the diCal hidden states are a triplet (epoch, ancientDeme, haplotype)."),
		new Switch ("useEigenCore", JSAP.NO_SHORTFLAG, "useEigenCore", "If this switch is set, instead of using the default ODECore (transitions base on ODEs), use the EigenCore (transitions based on spectral decomposition of migMatrix)."),

//        new Switch ("vcfMissingReference", JSAP.NO_SHORTFLAG, "vcfMissingReference",
//        		"If this switch is set, all non-segregating sites are assumed missing. Specifying this overrides whatever is provided in the vcf-file.")

	        };
		
		
		// only the visible parameters
		List<Parameter> paramVisibleList = new ArrayList<Parameter>();
		for (Parameter param : defaultParams) paramVisibleList.add(param);
		for (Parameter param : additionalParams) paramVisibleList.add(param);
			
		// all parameters		
		List<Parameter> paramCombinedList = new ArrayList<Parameter>();
		for (Parameter param : defaultParams) paramCombinedList.add(param);
		for (Parameter param : additionalParams) paramCombinedList.add(param);
		for (Parameter param : defaultHiddenParams) paramCombinedList.add(param);
		for (Parameter param : additionalHiddenParams) paramCombinedList.add(param);
		

		// first get the visible thing
		// sort the parameters, by ID is not necessarily by parameter name, but it's best we can do for now
		paramVisibleList.sort(Comparator.comparing(Parameter::getID));
//		for (Parameter param : paramVisibleList) System.out.println (param.getID());
		// and create fake object
		SimpleJSAP visibleJsap = new SimpleJSAP( 
				"DiCalParamSet", 
	            "Parses data and params for DiCal applications",
	            paramVisibleList.toArray(new Parameter[0])
	            	);
		
		// get help and usage from this
		String visibleHelp = visibleJsap.getHelp();
		String visibleUsage = visibleJsap.getUsage();
		
		// now build the real parsing object
		paramCombinedList.sort(Comparator.comparing(Parameter::getID));
		SimpleJSAP jsap = new SimpleJSAP( 
				"DiCalParamSet", 
	            "Parses data and params for diCal2",
	            paramCombinedList.toArray(new Parameter[0])
	            	);
		
		String hiddenHelp = jsap.getHelp();
		
		// set help and usage manually
		jsap.setHelp(visibleHelp);
		jsap.setUsage(visibleUsage);
	
		// we have to do this here manually, because otherwise the parsing complains about missing arguments
		for (String thisArg : args) {
			if (thisArg.equals("--hidden")) {
				System.err.println(hiddenHelp);
				System.exit(1);
			}
		}
		
		// and now stuff should just work
		if (args.length > 0) {
			this.jsapParams = jsap.parse(args);
		}
		else {
			this.jsapParams = null;
			System.err.println ("Use --help for more options.");
			System.exit(1);
		}
		
		if (jsap.messagePrinted()) { System.exit(1); }
		
		if (jsapParams.getBoolean("oldCakeRenormalizeSettings")) {
			UberDemographyCore.RENORMALIZE_PROBS = false;
			DemoState.DemoStateFactory.OLD_CAKE_INTERVALS = true;
		}
		
		outStream.println("# Parameter values:");

		// something about a fraction for the nelder mead stepsize
		if (this.jsapParams.contains("nmFraction")) {
			this.nelderMeadFraction = this.jsapParams.getDouble("nmFraction");
		}
		else {
			this.nelderMeadFraction = null;
		}

		if (this.jsapParams.contains("vcfFilterPassString")) {
			this.filterPassString = this.jsapParams.getString ("vcfFilterPassString").trim();
		}
		else {
			this.filterPassString = null;
		}

		
		
		// likelihood things must go here
		if(jsapParams.contains("compositeLikelihood")){
			String compositeLikelihoodString = jsapParams.getString("compositeLikelihood");
			this.compositeLikelihoodStringLow = compositeLikelihoodString.toLowerCase();


			// get the composite likelihood type, and initialize the list of CSDConfigs
			// first get the type
			this.loadedCsdList = null;
			if (this.compositeLikelihoodStringLow.equals("lol")) {
				this.compositeLikelihood = CompositeLikelihoodType.LOL;
			} else if (this.compositeLikelihoodStringLow.equals("pac")) {
				this.compositeLikelihood = CompositeLikelihoodType.PAC;
			} else if (this.compositeLikelihoodStringLow.equals("oldpac")) {
				this.compositeLikelihood = CompositeLikelihoodType.OLD_PAC;
			} else if (this.compositeLikelihoodStringLow.equals("superpac")) {
				this.compositeLikelihood = CompositeLikelihoodType.SuperPAC;
			} else if (this.compositeLikelihoodStringLow.equals("pcl")) {
				// pairwise composite likelihood
				this.compositeLikelihood = CompositeLikelihoodType.PCL;
			} else if (this.compositeLikelihoodStringLow.equals("pclol")) {
				this.compositeLikelihood = CompositeLikelihoodType.PCLOL;
			} else {
				// other than the stuff before has to be a file that's given
				this.compositeLikelihood = CompositeLikelihoodType.FILE;

				// load the CSDs from the file
				this.loadedCsdList = loadCSDsFromFile (compositeLikelihoodString);
			}
		} else {
			this.compositeLikelihoodStringLow = null;
			this.compositeLikelihood = null;
		}
		
		// read in recombination and mutation parameters
		// possibly a list to do chromosome specific
		if (!jsapParams.contains("paramFile")) {
			System.err.println("Error: Parameter 'paramFile' is required.");
			System.exit(1);
		}
		String[] parameterFileNameList = jsapParams.getString("paramFile").split(",");
		if (parameterFileNameList.length < 1) throw new IOException("No parameter file given.");
		List<FSAParamSet> twoLociParamSets = new ArrayList<FSAParamSet>();
		// all param sets should have the same numAlleles
		int consensusNumAlleles = -1;
		double minParamRate = Double.POSITIVE_INFINITY;
		
		for (int i=0; i<parameterFileNameList.length; i++) {
			String parameterFile = parameterFileNameList[i];
			
			// store it in a temporary paramSet for now, because we don't know how many loci there will be until we read in the sequences
			// also, we might modify the mutation stuff
			boolean renormMutMatrix = !jsapParams.getBoolean("useRawMutationMatrix");
			FSAParamSet tmpParamSet = FSAParamSetReader.readParams (new FileReader(parameterFile), 2, UberDemographyCore.ONELOCUS_EPSILON, renormMutMatrix);
			double mutRate = tmpParamSet.getMutationRate(0);
			double recoRate = tmpParamSet.getRecombinationRate(0, 1);
			Matrix mutMatrix = tmpParamSet.getMutationMatrix(0);
			int numAlleles = tmpParamSet.numAlleles();
			// only first time
			if (consensusNumAlleles < 0) {
				consensusNumAlleles = numAlleles;
			}
			if (consensusNumAlleles != numAlleles) throw new IOException ("Parameter file " + i + " has wrong number of alleles.");
			if (numAlleles < 2) throw new IOException ("Must have at least 2 alleles.");
			if ((numAlleles != mutMatrix.getColumnDimension()) || (numAlleles != mutMatrix.getRowDimension())) {
				throw new IOException("Number of alleles given in config file incompatible with dimensions of mutation matrix.");
			}
				
			// if your matrix is not uniform (josh's symmetric), then get the rate and matrix that produces the same number of non-self-mutations
			if (jsapParams.getBoolean("makeUniformMutationMatrix")) {
		//			tmpParamSet = new FSAParamSet(2, tmpParamSet.numAlleles(), getAvgMutationRate(tmpParamSet.getMutationRate(0), tmpParamSet.getMutationMatrix(0)), new Matrix(getUniformMutMatrix(tmpParamSet.numAlleles())), tmpParamSet.getRecombinationRate(0, 1));
				mutRate = getAvgMutationRate(mutRate, mutMatrix);
				mutMatrix = new Matrix(getUniformMutMatrix(numAlleles));
			}
			
			// remember the minimum
			minParamRate = Math.min(minParamRate, mutRate);
			minParamRate = Math.min(minParamRate, recoRate);
			
			// and remember it
			twoLociParamSets.add(new FSAParamSet(2, numAlleles, mutRate, mutMatrix, recoRate));
		}

		
		// read in the structured configuration
		if (jsapParams.contains("configFile") == jsapParams.contains("numPerDeme")) {
			throw new IOException ("Exactly one of configFile or numPerDeme must be specified.");
		}
		
		ConfigInfo configInfo = null;
		int[] numPerDeme = null;
		if (jsapParams.contains("configFile")) {
			
			if (requireOrderedConfig) {
				throw new IOException ("An ordered configuration sample is required, and --configFile specifies an unordered configuration.");
			}
			
			String configFile = jsapParams.getString("configFile");
			
			Reader configReader = new FileReader (configFile);
			configInfo = DemoConfiguration.readConfigInfo(new BufferedReader(configReader));
			
			// only in vcf mode can we be ok with no good number given
			assert ((configInfo.numLoci != null) || jsapParams.contains("vcfFile"));
						
			// in case of vcf, we have to reset it
			if (jsapParams.contains("vcfFile")) {
				configInfo.numLoci = null;
			} 
			
			
			//safety first!
			if(jsapParams.contains("bedFile")){
				if(!jsapParams.contains("vcfFile")){
					throw new IOException("Using BED files for masking has only been implemented for the vcf reader, so you must also specify vcfFile");
				}
			}
			
			
			configReader.close ();

			if (consensusNumAlleles != configInfo.numAlleles) throw new IOException("Parameter and configuration file have differing number of alleles");
			
			outStream.println("# configFile = " + configFile);
		} else {
			String numPerDemeString = jsapParams.getString("numPerDeme");
			
			String[] countsStrings = numPerDemeString.split(",");
			
			numPerDeme = new int[countsStrings.length];
			for (int i = 0; i < numPerDeme.length; i++){
				numPerDeme[i] = Integer.parseInt(countsStrings[i]);
			}
			
			outStream.println("# numPerDeme = " + numPerDemeString);
		}
		
		// here is a good place
		List<Boolean> sequencesToRead = new ArrayList<Boolean>();

		Integer additionalHap = null;
		if (jsapParams.contains("conditionalHap")) {
			additionalHap = jsapParams.getInt("conditionalHap");
		}

		if (configInfo != null) {
			
			// go through the multiplicities and see which sequences we need
			for (int multIdx=0; multIdx<configInfo.multiplicities.size(); multIdx++) {
				int[] currMult = configInfo.multiplicities.get(multIdx);
				int sum = 0;
				for (int i : currMult) {
					assert (i >= 0);
					sum += i;
				}
				// decide
				if (((additionalHap != null) && (additionalHap == multIdx)) || (sum > 0)) {
					// we want it
					sequencesToRead.add (true);
				}
				else {
					// no
					sequencesToRead.add (false);
				}
			}
		}
		else {
			assert (numPerDeme != null);
			assert (!jsapParams.contains("conditionalHap"));
			
			// just read all the sequences
			for (int thisMult : numPerDeme) {
				for (int i=0; i<thisMult; i++) {
					sequencesToRead.add(true);
				}
			}
		}

		
		// factory for the trunk
		String cakeName = jsapParams.getString("cakeStyle");
		CakeStyle cakeStyle = TrunkProcess.getCakeStyle(cakeName);
		String trunkStyle = jsapParams.getString("trunkStyle");
		if ((!(trunkStyle.equals("simple") || trunkStyle.equals("exactCake")) && cakeStyle == null) || (trunkStyle.equals("recursive") && cakeStyle == CakeStyle.MIDDLE) || (trunkStyle.equals("migratingEthan") && cakeStyle == CakeStyle.MIDDLE)) {
			throw new IOException("Invalid cake style for " + trunkStyle);
		}
		
		this.trunkFactory = new TrunkProcessFactory(trunkStyle, cakeStyle);
		if (trunkFactory.trunkStyle == null) {
			throw new IOException("Invalid trunk style");
		}
		
		
		// see whether we have bounds for the demography parameters
		double[][] bounds = null;
		if (jsapParams.contains("bounds")) {
			String boundString = jsapParams.getString("bounds");
			// pair are separated by semi-colons
			String[] boundPairs = boundString.split(";");
			bounds = new double[boundPairs.length][2];
			for (int i=0; i<boundPairs.length; i++) {
				String[] thisPair = boundPairs[i].split(",");
				if (thisPair.length != 2) throw new IOException("Each pair in bounds has to be two doubles separated by a comma");
				bounds[i][0] = Double.parseDouble(thisPair[0]);
				bounds[i][1] = Double.parseDouble(thisPair[1]);
			}
		}

		
		// read in the core type
		this.useEigenCore = jsapParams.getBoolean("useEigenCore");
		if (this.useEigenCore && jsapParams.contains("ratesFile")) {
			throw new IOException ("Cannot use EigenCore with exponential growth. Use ODECore (default) instead.");
		}
		

		
		// read in demography factory
		String demoParamsFile = null;
		String growthRatesFile = null;

		// get the demo filename
		demoParamsFile = jsapParams.getString("demoFile");
		
		// do we have a rates file?
		FileReader ratesReader = null;
		if (jsapParams.contains("ratesFile")) {
			// get the file
			growthRatesFile = jsapParams.getString ("ratesFile");
			ratesReader = new FileReader (growthRatesFile);
		}
		
		// and now build the right factory
		this.demoFactory = new DemographyFactory.ParamFileDemoFactory (new FileReader(demoParamsFile), ratesReader, UberDemographyCore.DEMO_FACTORY_EPSILON, this.trunkFactory.halfMigrationRate(), bounds, this.useEigenCore);
		
		if ((bounds != null) && (this.demoFactory.getBounds().length != bounds.length)) throw new IOException ("Bounding for reco estimation or mutation estimation not implemented yet.");
		
		// adjust the migration, if requested
		if (jsapParams.contains("minimalMigration")) {
			double factor = jsapParams.getDouble("minimalMigration");
			if (factor <= 0) throw new IOException("Negative minimal migration factor makes no sense");
			// now adjust it
			double minimalRate = factor * Math.min (minParamRate, this.demoFactory.getMinNonZeroMigrationRate());
			System.out.println("# minimal migration rate: " + minimalRate);
			this.demoFactory.addMinimalMigration (minimalRate);
		}
		
		// get a fancy locus transition map
		
		boolean useLocusSkipping = jsapParams.getBoolean("useLocusSkipping");
		int lociPerHmmStep = jsapParams.getInt("lociPerHmmStep");
		if (lociPerHmmStep < 1) {
			throw new IOException ("lociPerHmmStep must be positive.");
		}
		else {
			if (lociPerHmmStep > 32765) {
				throw new IOException ("lociPerHmmStep bigger than what fits into a short int.");
			}
		}
		if (lociPerHmmStep > 1 && useLocusSkipping) {
			throw new IOException ("Locus skipping not compatible with multiple loci per HMM step.");
		}
		
		// get the factory for the demoStates
		int addTrunkIntervals = jsapParams.getInt("addTrunkIntervals");
		String intervalFactoryParamsInput = jsapParams.getString("intervalParams").replaceAll("\\s","");
		String demoStatesMethod = "";
		String demoStatesParams = null;
		
		if (jsapParams.contains("intervalType")) {
			if (jsapParams.getBoolean("ancientDemeStates")) {
				throw new IOException("Must specify exactly one of intervalType OR ancientDemeStates");
			}
			
			boolean printIntervals = jsapParams.getBoolean ("printIntervals");
			
			// get the interval factory
			String intervalTypeName = jsapParams.getString("intervalType").replaceAll("\\s","").toLowerCase();
			IntervalFactory intervalFactory = IntervalFactory.getIntervalFactory (intervalTypeName, intervalFactoryParamsInput, printIntervals);
			// everything all right?
			if (intervalFactory == null) {
				throw new IOException ("Could not initialize interval Factory!");
			}
			
			//mesh some dudes if estimating theta
			if(this.jsapParams.contains("estimateTheta") && this.jsapParams.getBoolean("estimateTheta")){
				List<IntervalFactory> intFacts = new ArrayList<IntervalFactory>();
				
				double lastPoint = 0d;
				String factoryParams = this.jsapParams.getString("thetaIntervalParams");
				String[] intervalFactoryParams = factoryParams.split(",");
				double[] points = new double[intervalFactoryParams.length];
				if(intervalFactoryParams.length < 1) throw new IOException ("You messed up the thetaIntervalParams");
				for (int i = 0; i < intervalFactoryParams.length; i++) {
					double nextPoint = Double.parseDouble(intervalFactoryParams[i]);
					if (nextPoint <= lastPoint) throw new IOException("Something in your thetaIntervalParams is messed up");
					points[i] = nextPoint;
					lastPoint = nextPoint;
				}
				intFacts.add(new IntervalFactory.CustomFixedIntervalFactory (points, false));
				intFacts.add(intervalFactory);
				intervalFactory = new IntervalFactory.CombinedIntervalFactory(intFacts, UberDemographyCore.ONELOCUS_EPSILON, printIntervals);
			}
			
			this.intervalFactory = intervalFactory;
			// for some special factories, we want the intervals now, regardless of printIntervals or not
			if ((intervalFactory instanceof LogUniformFactory) || intervalFactory instanceof CustomFixedIntervalFactory) {
				System.out.println ("# realIntervals:\t" + Arrays.toString (intervalFactory.getNontrivialIntervals (null, null, -1)));
			}
			
			this.demoStateFactory = new DemoStateFactory(intervalFactory, addTrunkIntervals, UberDemographyCore.ONELOCUS_EPSILON);
			demoStatesMethod = "intervalsOnly-" + intervalFactory.factoryType();
			demoStatesParams = intervalFactory.factoryParams();
		} else {
			this.intervalFactory = null;
			if (!jsapParams.getBoolean("ancientDemeStates")) {
				throw new IOException("Must either specify an intervalType, or specify to use ancientDemeStates.");
			}
			
			if (!intervalFactoryParamsInput.equals("")) {
				throw new IOException("Cannot specifiy intervalParams when using ancientDemeStates.");
			}
			
			// get the demoStateFactory that create demoStates for each (epoch, ancientDeme, presentDeme)
			this.demoStateFactory = new DemoStateFactory(null, addTrunkIntervals, UberDemographyCore.ONELOCUS_EPSILON);
			
			demoStatesMethod = "ancientDemeStates";
			demoStatesParams = null;
		}

		
		// get time
		long oldTime = System.currentTimeMillis();
		System.out.println ("# Reading input sequences.");

		
		// need to know this
		// do we want to print every K?
		Integer printEveryK = null;
		if (this.jsapParams.contains("printLoci")) {
			// check it
			if (lociPerHmmStep > 1) {
				throw new IOException ("Cannot specify print loci when using the multilocus-blocking speedup.");
			}
			// get value
			printEveryK = this.jsapParams.getInt("printLoci");
		}
		// what about the snps
		boolean printSnps = jsapParams.getBoolean("printSnps");

		
		// get the alleles
		char[] alleles = null;
		if (jsapParams.contains("alleles")) {
			alleles = getAlleles (jsapParams.getString("alleles"), consensusNumAlleles);
		}
		
		// the comment chars
		char[] commentChars = new char[] {'#', '>'};
		
		// and the missing
		char[] missingChars = null;
		if (jsapParams.contains("missingAlleles")) {
			missingChars = getAlleles (jsapParams.getString("missingAlleles"), null);
		}
		else {
			// default ones
			missingChars = new char[] {'n', 'N', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V'};
		}
		
		// to vcf or not to vcf
		String[] sequenceFileNameList = null;
		Boolean compressedDicalSequence = null;
		List<Integer> vcfOffsets = null;
		boolean readDical = false, readVCF = false;
		if (jsapParams.contains("sequenceFile")) {
			assert ((configInfo == null) || (configInfo.numLoci != null));

//			if (jsapParams.getBoolean("vcfDiAllelic")) throw new IOException ("vcfDiAllelic only allowed for in VCF input.");
//			if (jsapParams.getBoolean("vcfTriAllelic")) throw new IOException ("vcfTriAllelic only allowed for in VCF input.");
			if (jsapParams.contains("vcfOffset")) throw new IOException ("vcfOffset only allowed for in VCF input.");

			
			readDical = true;
			
			// get them names
			sequenceFileNameList = jsapParams.getString("sequenceFile").split(",");
			compressedDicalSequence = !jsapParams.getBoolean("unCompressedSequence");
		}
		else if (jsapParams.contains("vcfFile")) {
			assert ((configInfo == null) || (configInfo.numLoci == null));
			
			if (jsapParams.contains("alleles")) throw new IOException ("Specifying alleles is not allowed in VCF format.");
			if (jsapParams.contains("missingAlleles")) throw new IOException ("Specifying missingAlleles is not allowed in VCF format.");

//			if (jsapParams.getBoolean("vcfDiAllelic")) {
//				if (consensusNumAlleles != 2) throw new IOException("If vcfDiAllelic is set, the number of alleles given in the config file has to be two.");
//				System.out.println ("# [WARNING] The vcf file will be represented internally as diallelic. Incompatible sites will be ignored. Allele numbers are as in the vcf file");
//			}
//			if (jsapParams.getBoolean("vcfTriAllelic")) {
//				if (consensusNumAlleles != 3) throw new IOException("If vcfTriAllelic is set, the number of alleles given in the config file has to be three.");
//				System.out.println ("# [WARNING] The vcf file will be represented internally as triallelic. Incompatible sites will be ignored. Allele numbers are as in the vcf file");
//			}
//			if (jsapParams.getBoolean("vcfDiAllelic") && jsapParams.getBoolean("vcfTriAllelic")) {
//				throw new IOException("Only one of vcfDiAllelic or vcfTriAllelic allowed.");
//			}
			if (consensusNumAlleles < 4) {
				// not sure why this is necessary
				boolean allUniform = true;
				for (FSAParamSet thisSet : twoLociParamSets) {
					allUniform &= SimpleLocusSkipper.isUniformMutationMatrix(thisSet.getMutationMatrix(0).getArray(), UberDemographyCore.ONELOCUS_EPSILON);
				}
				if (!allUniform) {
					throw new IOException ("If you have less than 4 alleles for the vcf file, you should have a uniform mutation matrix.");
				}
			}
			if (consensusNumAlleles > 4) throw new IOException ("Maximally 4 alleles allowed for vcf.");
			
			assert (consensusNumAlleles == configInfo.numAlleles);
			
			// see about offset for vcf files
			if (jsapParams.contains("vcfOffset")) {
				String[] offsets = jsapParams.getString("vcfOffset").split(",");
				vcfOffsets = new ArrayList<Integer>();
				for (String offset : offsets) vcfOffsets.add(Integer.parseInt(offset));
			}

			
			
//			if (lociPerHmmStep <= 1) throw new IOException ("For now we only use VCF with the MultiLocusStepHandler.");
			if (jsapParams.getBoolean("unCompressedSequence")) throw new IOException ("Un-compressed sequence doesn't make sense for vcf.");
						
			readVCF = true;
			
			if (this.filterPassString == null) {
				throw new IOException ("You did not proved a --vcfFilterPassString");
			}
			
			
			// get them names
			sequenceFileNameList = jsapParams.getString("vcfFile").split(",");

			// this should be if we do vcf
			if (configInfo == null) throw new IOException("VCF input format needs a config file.");
		}
		else {
			// tertium non datur
			System.err.println ("Provide either sequence files or vcf files.");
			System.exit (-1);
		}

		// if we have more than one paramset, better have as many as sequences
		if (twoLociParamSets.size() > 1) {
			if (twoLociParamSets.size() != sequenceFileNameList.length) throw new IOException ("Number of parameter files is not equal to number of sequence files.");
		}
		else {
			// just make as many copies of paramSets as their are sequences, for convencience
			for (int j=1; j < sequenceFileNameList.length; j++) {
				twoLociParamSets.add(twoLociParamSets.get(0));
			}
		}
		// better be true now
		assert (twoLociParamSets.size() == sequenceFileNameList.length);
		
		
		// how about the reference file
		// in vcf or not, how to deal with it
		// also, should be as many files as sequenceFileNameList
		boolean referenceFileOnCommandLine = jsapParams.contains("vcfReferenceFile");
//		boolean wantFixedReference = jsapParams.getBoolean("vcfFixedReference");
//		boolean wantRandomReference = jsapParams.getBoolean("vcfRandomReference");
//		boolean wantMissingReference = jsapParams.getBoolean("vcfMissingReference");

		// init
		String[] vcfReferenceFilenameList = null;
//		ReferenceMode referenceMode = ReferenceMode.NONE;

		// and we only care if we want to read vcfs
		if (!readVCF) {
			if (referenceFileOnCommandLine) {
//			if (referenceFileOnCommandLine || wantFixedReference || wantRandomReference || wantMissingReference) {
//				String vcfClashString = vcfIncompatibleString (referenceFileOnCommandLine, wantFixedReference, wantRandomReference, wantMissingReference); 
//				throw new IOException ("The parameter(s) " + vcfClashString + " can only be used if the input sequences are given in a vcfFile.");
//			}
				throw new IOException ("The parameter --vcfReferenceFile can only be used if the input sequences are given in a vcfFile.");
			}
		}
		else {
			assert (readVCF);
			// now we wants vcf
			// but only at most one option
			
//			boolean compatible = true;
			
			// go through
			if (referenceFileOnCommandLine) {
//				if (wantFixedReference || wantRandomReference || wantMissingReference) {
//					compatible = false;
//				}

				vcfReferenceFilenameList = jsapParams.getString("vcfReferenceFile").split(",");
				// must be as many as we have sequence files
				if (vcfReferenceFilenameList.length != sequenceFileNameList.length) {
					throw new IOException("Must have as many reference-files as vcf-files (one for each).");
				}
				// all good
			}
//			else if (wantFixedReference) {
//				if (referenceFileOnCommandLine || wantRandomReference || wantMissingReference) {
//					compatible = false;
//				}
//				
//				referenceMode = ReferenceMode.FIXED;
//			}
//			else if (wantRandomReference) {
//				if (referenceFileOnCommandLine || wantFixedReference || wantMissingReference) {
//					compatible = false;
//				}
//				
//				referenceMode = ReferenceMode.RANDOM;
//			}
//			else if (wantMissingReference) {
//				if (referenceFileOnCommandLine || wantFixedReference || wantRandomReference) {
//					compatible = false;
//				}
//				
//				referenceMode = ReferenceMode.MISSING;
//			}
//			
//			//  are we gucci?
//			if (!compatible) {
//				String vcfClashString = vcfIncompatibleString (referenceFileOnCommandLine, wantFixedReference, wantRandomReference, wantMissingReference); 
//				throw new IOException ("The parameter(s) " + vcfClashString + " are incompatible. Use only one of them.");
//			}
		}
		
		
		
		
		// bed stuff
		//
		
		//safety first!
		if(jsapParams.contains("bedFile")){
			if(!jsapParams.contains("vcfFile")){
				throw new IOException("Using BED files for masking has only been implemented for the vcf reader, so you must also specify vcfFile");
			}
		}

		// get the bed file names
		String[] bedFileNameList = null;
		
		if(jsapParams.contains("bedFile")){
			bedFileNameList = jsapParams.getString("bedFile").split(",");
			// must be as many as we have sequence files
			if (bedFileNameList.length != sequenceFileNameList.length) {
				throw new IOException("Must have as many bed-files as vcf-files (one for each).");
			}
		}

		
		// read in the sequences, and store here
		boolean useStationaryForPartially = jsapParams.getBoolean("stationaryForPartiallyMissing");
		this.extendedConfigInfoList = new ArrayList<ExtendedConfigInfo>();
		if (vcfOffsets != null) {
			if ((vcfOffsets.size() != 1) && (vcfOffsets.size() != sequenceFileNameList.length)) {
				throw new IOException ("Either provide one offset for all files, or one for each.");
			}
		}
		for (int j=0; j<sequenceFileNameList.length; j++) {
			String sequenceFile = sequenceFileNameList[j];
			// open the stream
			Reader sequenceReader = new FileReader (sequenceFile);
			// get the directory
			System.out.println ("# sequenceFile: " + sequenceFile);
			// more robust
			String pathToSequenceFile = (new File ((new File(sequenceFile)).getAbsolutePath()).getParentFile().getAbsolutePath());
						// here to store the guys
			List<HFSAXFullHaplotype> sequenceList = null;
			
			// which format to read?
			if (readDical) {
				
				if (!compressedDicalSequence) {
					sequenceList = ReadSequences.xReadFullSequenceFormat(sequenceReader, sequencesToRead, consensusNumAlleles, alleles, commentChars, missingChars);
				} else {
					if (this.jsapParams.getBoolean("oldCompressedDataFormat")) {
						if (!this.jsapParams.contains("configFile")) {
							sequenceReader.close();
							throw new IOException ("Old compressed data format required a config file to also be specified.");
						}
						
						ConfigInfo tmpConfigInfo = DemoConfiguration.readConfigInfo (new BufferedReader (new FileReader(this.jsapParams.getString("configFile"))));
						sequenceList = ReadSequences.xReadOldCompressedSequenceFormat(sequenceReader, sequencesToRead, consensusNumAlleles, alleles, tmpConfigInfo.numLoci, commentChars, missingChars);
					} else {
						sequenceList = ReadSequences.xReadCompressedSequenceFormat(sequenceReader, sequencesToRead, consensusNumAlleles, alleles, commentChars, missingChars);
					}
				}
				// should be good
			}
			else {
				assert (readVCF);
				
				int offset = 0;
				if (vcfOffsets != null) {
					offset = vcfOffsets.get(0);
					if (vcfOffsets.size() > 1) offset = vcfOffsets.get(j);
				}				
				// do bed stuff
				Reader bedReader = null;
				if (bedFileNameList != null) {
					String bedFile = bedFileNameList[j];
					System.out.println ("# bedFile: " + bedFile);
					bedReader = new FileReader (bedFile);
				}
				
				
				// now read them vcf files
				// we want them phased
				String vcfReferenceFilename = null;
				if (vcfReferenceFilenameList != null) {
					assert (vcfReferenceFilenameList.length > j);
					vcfReferenceFilename =  vcfReferenceFilenameList[j];
				}
				sequenceList = ReadSequences.readVcf (sequenceReader, sequencesToRead, consensusNumAlleles, commentChars, this.jsapParams.getBoolean("acceptUnphasedAsMissing"), pathToSequenceFile, offset, bedReader, true, this.filterPassString, vcfReferenceFilename, this.jsapParams.getBoolean("vcfIgnoreDoubleEntries"));			
			}
			assert (sequenceList != null);
			// this should hold for all
//			if (!readVCF) {
			if (sequenceList.size() != sequencesToRead.size()) throw new IOException("Number of sequences in seqenceFile doesn't match number given in config file.");
//			}

			// check the number of loci given in the config file
			if ((configInfo != null) && (configInfo.numLoci != null)) {
				for (int i=0; i<sequenceList.size(); i++) {
					// some sequences are not read and thus point to null
					if ((sequenceList.get(i) != null) && (sequenceList.get(i).getNumLoci() != configInfo.numLoci)) {
						throw new IOException("Sequence length (number of loci) given in config file [" + configInfo.numLoci + "] differs from length of haplotype in sequence file [" + sequenceList.get(i).getNumLoci() + "].");
					}
				}
			}

//			System.out.println(Arrays.toString(((SimpleFSARef)sequenceList.get(0).getReference()).getAlleleConfig()));
//			for (int i=0; i<sequenceList.size(); i++) {
//				System.out.println(sequenceList.get(i));
//			}
			
			// be nice, close the stream
			sequenceReader.close ();

			// print every K loci maybe
			TreeSet<Integer> printLoci = null;
			if (printEveryK != null) {
				printLoci = makePrintLoci(sequenceList.get(0).getNumLoci(), printEveryK);
			}
			// and maybe the SNPs as well
			if (printSnps) {
				HashSet<HFSAXFullHaplotypeShell> fakeSet = new HashSet<HFSAXFullHaplotypeShell>();
				for (HFSAXFullHaplotype hap : sequenceList) fakeSet.add (new HFSAXFullHaplotypeShell(hap));
				TreeSet<Integer> snpLoci = SimpleXLocusSkipper.getSNPs(fakeSet);
				if (printLoci == null) {
					printLoci = snpLoci;
				} else {
					printLoci.addAll(snpLoci);
				}
			}
			
			
			// check the loci
			// make a local config info
			ConfigInfo thisConfigInfo = null;
			if (configInfo != null) {
				if (configInfo.numLoci == null) {
					assert (readVCF);
					// catch up
					int firstRealHapIdx = 0;
					while (sequenceList.get(firstRealHapIdx) == null) {
						firstRealHapIdx++;
						assert (firstRealHapIdx != sequenceList.size());
					}
					// make a config
					thisConfigInfo = new ConfigInfo (configInfo.multiplicities,configInfo.numDemes, sequenceList.get(firstRealHapIdx).getNumLoci(), configInfo.numAlleles);
				}
				else {
					assert (readDical);
					// have to adjust
					int i=0;
					while (sequenceList.get(i) == null) {
						i++;
						assert (i != sequenceList.size());
					}
					
					if (sequenceList.get(i).getNumLoci() != configInfo.numLoci) {
						throw new IOException ("The number of loci given in the config file doesn't match the sequence files");
					}
					thisConfigInfo = new ConfigInfo (configInfo.multiplicities,configInfo.numDemes, configInfo.numLoci, configInfo.numAlleles);
				}
			}
			
			// we should check somewhere, but still allow people to use the --numPerDeme flag, no???
			if (!jsapParams.contains("numPerDeme") && thisConfigInfo.numDemes != this.demoFactory.getNumPresentDemes()) throw new IOException("Demography has to have same number of demes at present as config file");
			double mutRate = twoLociParamSets.get(j).getMutationRate(0);
			double recoRate = twoLociParamSets.get(j).getRecombinationRate(0, 1);
			Matrix mutMatrix = twoLociParamSets.get(j).getMutationMatrix(0);
//			System.out.println ("(" + mutRate + ", " + recoRate + ")");
			ExtendedConfigInfo thisConfig = new ExtendedConfigInfo(sequenceList, thisConfigInfo, numPerDeme, mutRate, recoRate, mutMatrix, useLocusSkipping, lociPerHmmStep, printLoci, useStationaryForPartially, additionalHap, loadedCsdList);
			this.extendedConfigInfoList.add (thisConfig);
			// check whether they are ok with a potential csd file
			if (this.loadedCsdList != null) {
				assert (this.compositeLikelihood == CompositeLikelihoodType.FILE);
				// load the CSDs from the file
				checkNumHapsInCSDs (this.loadedCsdList, thisConfig.demeHaplotypeList.size());
			}
		}
		// sequences should be fine now
		
		
		// timing
		long currTime = System.currentTimeMillis();
		long readTime = currTime - oldTime;
		System.out.println ("# " + readTime + " milliseconds needed to read the sequence files (including building step handler).");
		
		
		// go on with objective function stuff
		this.conditionalObjectiveFunction = ConditionalObjectiveFunctionType.ConditionLineage;
		boolean condSwitch = jsapParams.getBoolean("condOnTransitionType");
		boolean klSwitch = jsapParams.getBoolean("marginalKL");
		if (condSwitch && klSwitch) {
			throw new IOException("Can only set one of switches, condOnTransitionType and marginalKL");
		}
		if (condSwitch) {
			this.conditionalObjectiveFunction = ConditionalObjectiveFunctionType.ConditionLineageTransitionType;
		}
		if (klSwitch) {
			this.conditionalObjectiveFunction = ConditionalObjectiveFunctionType.MarginalKLDivergence;
		}

		
		
		// print out the command line arguments
		assert (outStream != null);
		
		DateFormat dateFormat = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss");
		Date date = new Date();
		System.out.println("# Date/time: " + dateFormat.format(date));
		

		assert ((parameterFileNameList.length == 1) || (parameterFileNameList.length == twoLociParamSets.size()));
		for (int m=0; m<twoLociParamSets.size(); m++) {
			if (parameterFileNameList.length == 1) {
				outStream.println("# parameterFile = " + parameterFileNameList[0]);
			}
			else {
				outStream.println("# parameterFile = " + parameterFileNameList[m]);
			}
			outStream.println ("# [PARAM_START]");
			PrintDicalBasicParamSet (outStream, twoLociParamSets.get(m));
			outStream.println ("# [PARAM_END]");
		}
		outStream.println("# demographyParameterFile = " + demoParamsFile);
		outStream.println("# growthRatesFile = " + growthRatesFile);
		outStream.println ("# [DEMO_START]");		
		this.demoFactory.dump (outStream, true);
		outStream.println ("# [DEMO_END]");		
		outStream.println("# locus skipping = " + useLocusSkipping);
		outStream.println("# trunkStyle = " + trunkStyle);
		outStream.println("# additionalTrunkIntervals = " + addTrunkIntervals);
		outStream.println("# cakeStyle = " + cakeName);
		outStream.println("# printEveryKLoci = " + printEveryK);
		outStream.println("# printSnps = " + printSnps);
		if (this.conditionalObjectiveFunction == ConditionalObjectiveFunctionType.MarginalKLDivergence) outStream.println("# Using marginal KL divergence");
		outStream.println("# demoStatesMethod = " + demoStatesMethod);
		outStream.println("# demoStatesParams = " + demoStatesParams);
		outStream.println("# partiallyMissing = " + useStationaryForPartially);
		
		for (int i = 0; i < this.extendedConfigInfoList.size(); i++) {
			// say something about config for first one
			if (i==0) {
				// we actually only care about the multiplicities
				int[] total = new int[this.extendedConfigInfoList.get(0).configInfo.multiplicities.get(0).length];
				for (int[] thisMult : this.extendedConfigInfoList.get(0).configInfo.multiplicities) {
					for (int j=0; j<total.length; j++) {
						total[j] += thisMult[j];
					}
				}
				outStream.println ("# [CONFIG_MULT] " + Arrays.toString(total));
			}
			outStream.println("# sequenceFile " + i + " = " + sequenceFileNameList[i]);
	    	// and show some locus transition stuff
			this.extendedConfigInfoList.get(i).fancyTransitionMap.printInfo (outStream);
			outStream.println("#");
		}
	}
		
	private String vcfIncompatibleString (boolean referenceFileOnCommandLine, boolean wantFixedReference, boolean wantRandomReference, boolean wantMissingReference) {
		List<String> stringList = new ArrayList<String>();
		if (referenceFileOnCommandLine) stringList.add ("--vcfReferenceFile");
		if (wantFixedReference) stringList.add ("--vcfFixedReference");
		if (wantRandomReference) stringList.add ("--vcfRandomReference");
		if (wantMissingReference) stringList.add ("--vcfMissingReference");
		return String.join(", ", stringList);
	}

	private void PrintDicalBasicParamSet (PrintStream outStream, FSAParamSet tmpParamSet) {
		outStream.println ("# mut = " + tmpParamSet.getMutationRate(0));
		outStream.println ("# rec = " + tmpParamSet.getRecombinationRate(0, 1));
		outStream.println ("# mutation matrix = ");
		double[][] daMatrix = tmpParamSet.getMutationMatrix(0).getArray();
		for (int i=0; i<daMatrix.length; i++) {
			outStream.println ("# " + Arrays.toString(daMatrix[i])); 
		}
	}

	public static char[] getAlleles (String alleleString, Integer numAlleles) throws IOException {
		String[] preAlleles = alleleString.trim().split(",");
		if ((numAlleles != null) && (preAlleles.length != numAlleles)) {
			throw new IOException ("Inputted alleles don't agree with number of alleles in paramFile.");
		}
		
		char[] alleles = new char[preAlleles.length];
		for (int i = 0; i < alleles.length; i++) {
			assert (preAlleles[i].length() == 1);
			alleles[i] = preAlleles[i].charAt(0);
		}
		
		return alleles;
	}

	private static TreeSet<Integer> makePrintLoci(int numLoci, int everyKbases) {
		TreeSet<Integer> toReturn = new TreeSet<Integer>();
		
		int currLocus = 0;
		while (currLocus < numLoci) {
			toReturn.add(currLocus);
			currLocus += everyKbases;
		}
		
		return toReturn;
	}
	
	public static double getAvgMutationRate(double mutRate, Matrix probMatrix){
		// get average per generation mutation rate
		double[] stationaryDistn = StationaryDistribution.getStationaryDistribution(probMatrix.getArray());
		double avgMutRate = 0d;
		for (int i = 0; i < stationaryDistn.length; i++){
			avgMutRate += mutRate * (1 - probMatrix.get(i, i)) * stationaryDistn[i];
		}
		return avgMutRate;
	}
	
	public static double[][] getUniformMutMatrix(int numAlleles) {
		double[][] toReturn = new double[numAlleles][numAlleles];
		for (int i = 0; i < numAlleles; i++) {
			for (int j = 0; j < numAlleles; j++) {
				if (i == j) continue;
				toReturn[i][j] = 1d / (numAlleles - 1);
			}
		}
		return toReturn;
	}
	
	private static void checkNumHapsInCSDs (List<List<Integer>> csdList, int numHaps) throws IOException {
		for (List<Integer> thisCSD : csdList) {
			for (Integer hapIdx : thisCSD) {
				if (hapIdx >= numHaps) throw new IOException("Invalid index in CSD");
			}
		}
	}
	
	@SuppressWarnings("resource")
	private static List<List<Integer>> loadCSDsFromFile (String fileName) throws IOException {
		
		Set<Character> commentCharacterSet = new TreeSet<Character>();
		for (char c : new char[] {'#'}) commentCharacterSet.add(c);

		// get a reader for that file
		BufferedReader bufferedCsdReader = new BufferedReader(new FileReader(fileName));
				
		List<List<Integer>> toReturn = new ArrayList<List<Integer>>();
		// go through the file
		String line = null;
		while ((line = bufferedCsdReader.readLine()) != null) {
			// always ignore comment and empty lines
			if (line.trim().equals("") || commentCharacterSet.contains(line.charAt(0))) continue;
			
			// now we should have a list of integers
			String[] fields = line.split("\\s+");
			List<Integer> thisList = new ArrayList<Integer>();
			toReturn.add(thisList);
			for (String value : fields) {
				Integer thisInt = new Integer(value);
				if (thisInt < 0) throw new IOException("Invalid index in CSD");
				thisList.add(thisInt);
			}
			
			// make sure about uniqueness of number in each row
			Set<Integer> set = new TreeSet<Integer>(thisList);
			if (set.size() != thisList.size()) throw new IOException("Duplicate entry in CSD");
		}
		// now we should have our list		
		
		// be nice
		bufferedCsdReader.close();
		
		// and return things
		return toReturn;
	}

}
