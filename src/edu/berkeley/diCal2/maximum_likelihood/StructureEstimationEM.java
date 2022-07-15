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

import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optimization.PointValuePair;
import org.apache.commons.math3.optimization.direct.NelderMeadSimplex;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.Switch;

import edu.berkeley.diCal2.csd.DemoConfiguration;
import edu.berkeley.diCal2.csd.DemoState.DemoStateCollection;
import edu.berkeley.diCal2.csd.DemoState.DemoStateFactory;
import edu.berkeley.diCal2.csd.EigenParamSet;
import edu.berkeley.diCal2.csd.HmmStepHandler;
import edu.berkeley.diCal2.csd.HmmStepHandler.SingleStepMissingAllelesTransitionMap.CoreCache;
import edu.berkeley.diCal2.csd.IntervalFactory;
import edu.berkeley.diCal2.csd.SMCSDemo;
import edu.berkeley.diCal2.csd.SMCSDemoEMObjectiveFunction;
import edu.berkeley.diCal2.csd.SMCSDemoEMObjectiveFunction.EmptyConditionalConfigOF;
import edu.berkeley.diCal2.csd.TrunkProcess;
import edu.berkeley.diCal2.csd.TrunkProcess.TrunkProcessFactory;
import edu.berkeley.diCal2.csd.UberDemographyCore;
import edu.berkeley.diCal2.csd.auxiliary.DelayedRandom;
import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.demography.DemographyFactory;
import edu.berkeley.diCal2.haplotype.FSAHaplotype;
import edu.berkeley.diCal2.haplotype.GeneticType.GeneticTypeMultiplicity;
import edu.berkeley.diCal2.haplotype.HFSAXFullHaplotype.HFSAXFullHaplotypeShell;
import edu.berkeley.diCal2.maximum_likelihood.DiCalParamSet.DemeHapPair;
import edu.berkeley.diCal2.maximum_likelihood.MetaOptimization.MetaOptimizationArguments;
import edu.berkeley.diCal2.utility.LogSum;
import edu.berkeley.diCal2.utility.Pair;
import edu.berkeley.diCal2.utility.RealPartition.Interval;
import edu.berkeley.diCal2.utility.SimplePair;
import edu.berkeley.diCal2.utility.VisitPermuations;
import edu.berkeley.diCal2.utility.VisitPermuations.Permutation;
import edu.berkeley.diCal2.utility.VisitPermuations.PermutationVisitor;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

public class StructureEstimationEM {

	// enums to choose between different types of analyses
	public static enum CompositeLikelihoodType {PAC, SuperPAC, LOL, PCL, OLD_PAC, FILE, PCLOL};
	public static enum ConditionalObjectiveFunctionType {ConditionLineageTransitionType, ConditionHaplotype, ConditionLineage, MarginalKLDivergence}

	private static final double NM_MIMUM_STEPSIZE = 1e-8;

	private static IntervalFactory thetaIntervalFactory;
	private static IntervalFactory otherIntervalFactory;
	
	public static void main (String[] args) throws JSAPException, IOException {
		
		long programStartTime = System.currentTimeMillis();

		// print out the command line arguments
		PrintStream outStream = System.out;
		outStream.print("# Command-line arguments: ");
		for (String arg: args)	{
			outStream.print(arg + " ");
		}
		outStream.print("\n");
				
        Parameter[] additionalParams = new Parameter[] {
        		new Switch ("verbose", 'v', "verbose", "Verbose for the EM procedure."),
       		
            	new FlaggedOption( "seed", JSAP.LONG_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "seed", 
                    "The seed to initialize the randomness." ),
                    
                new FlaggedOption ("numPermutations", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "numPermutations",
                		"Number of permutations to use for PAC-like methods."),
                new FlaggedOption ("permutationsFile", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "permutationsFile",
                		"Specfiy file(s) containing a list of permutations to be used to analyze each chunk respectively."),

                new FlaggedOption("startPoint", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "startPoint",
                		"The starting point for the EM."),                		
        		
                new FlaggedOption("numberIterationsEM", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "numberIterationsEM",
        			"EM terminates after the given number of iterations."),
        		new FlaggedOption("numberIterationsMstep", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "numberIterationsMstep",
    				"M-step search terminates after given number of iterations"),	
    			new FlaggedOption("parallel", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "parallel",
            			"Specifies number of parallel threads"),
            			
//				new Switch ("printEmPath", JSAP.NO_SHORTFLAG, "printEmPath", "If this switch is set, we print the parameter values at each EM step."),
            	new Switch ("disableCoordinateWiseMStep", JSAP.NO_SHORTFLAG, "disableCoordinateWiseMStep", "Default mode is a coordinatewise M Step optimization. The number of iterations and the reletive error thresholds are applied coordinatewise. If disabled, the M step optimization uses a general multidimensional NM algorithm."),
            new FlaggedOption("coordinateOrder", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "coordinateOrder",
                		"Order to update parameters for coordinateWiseMStep. If not specified, use random order.")          		

        };
        
        Parameter[] additionalHiddenParams = new Parameter[] {
        		new Switch ("estimateRecomScaling", JSAP.NO_SHORTFLAG, "estimateRecomScaling",
        				"If this switch is set, the EM will estimate a scaling factor for the recombination rate, so inputRecomRate * scaling = recomRate."),
        		
        		new Switch ("estimateTheta", JSAP.NO_SHORTFLAG, "estimateTheta",
        				"If this switch is set, the EM will estimate theta for each interval"),
        				
        		new FlaggedOption("thetaIntervals", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "thetaIntervals",
        	        			"factory type by which to build the theta intervals, recommend customFixed"),
        	    new FlaggedOption("thetaIntervalParams", JSAP.STRING_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "thetaIntervalParams",
        	        			"params for the thetaIntervals factory."),
        			
                new FlaggedOption("relativeErrorM", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "relativeErrorM",
        			"M-step optimization terminates if relative improvement is less then the given value."),  			
                new FlaggedOption("relativeErrorE", JSAP.DOUBLE_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "relativeErrorE",
        			"EM optimization terminates if relative likelihood improvement is less then the given value."),  			
    			new Switch ("useParamRelErr", JSAP.NO_SHORTFLAG, "useParamRelErr", "If this switch is set, we use the relative error of the parameters (instead of log-likelihood) as stopping criteria of E-step."),
        		new Switch ("diffPermsPerChunk", JSAP.NO_SHORTFLAG, "diffPermsPerChunk", "If this switch is set, use different set of permutations for each independent chunk."),
            new FlaggedOption ("numCsdsPerPerm", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "numCsdsPerPerm", "Number of CSDs used for each permutation."),
    			new Switch ("useParamRelErrM", JSAP.NO_SHORTFLAG, "useParamRelErrM", "If this switch is set, we use the relative error of parameters (instead of log-likelihood) as stopping criteria for M-step.")
        };

        // register some meta stuff
        MetaOptimizationArguments metaArguments = new MetaOptimizationArguments();
        
        additionalParams = metaArguments.register (additionalParams);
        additionalHiddenParams =  metaArguments.registerHidden (additionalHiddenParams);
        
		DiCalParamSet dicalParams = new DiCalParamSet (args, additionalParams, additionalHiddenParams, System.out);        
		
		boolean verbose = dicalParams.jsapParams.getBoolean("verbose");
		
		// double negative make positive
		boolean coordinateWiseM = !dicalParams.jsapParams.getBoolean("disableCoordinateWiseMStep");

		Integer numCsdsPerPerm = null;
        if (dicalParams.jsapParams.contains("numCsdsPerPerm")) {
        	numCsdsPerPerm = dicalParams.jsapParams.getInt ("numCsdsPerPerm");
        	if (numCsdsPerPerm < 1) throw new IOException ("Have to use at least one csd per permutation.");
        }
		
		// intialize the randomness
        // maybe this only sets a flag, and then actual randomness gets initialized when need
        Long seeed = null;
        if (dicalParams.jsapParams.contains("seed")) {
        		// we have a seed, so things looking gucci
        		seeed= new Long (dicalParams.jsapParams.getLong("seed"));
        }
        
        // and get the delayed random guy
		DelayedRandom randomGenerator = new DelayedRandom (seeed);
				
		// print it
//		boolean printFullEmPath = dicalParams.jsapParams.getBoolean ("printEmPath");
		boolean printFullEmPath = true;

		Double relativeErrorE = null;
		Integer numberIterationsEM = null;
		if (dicalParams.jsapParams.contains("relativeErrorE") && !dicalParams.jsapParams.contains("numberIterationsEM")) {
			relativeErrorE = dicalParams.jsapParams.getDouble ("relativeErrorE");
		}
		else if (dicalParams.jsapParams.contains("numberIterationsEM") && !dicalParams.jsapParams.contains("relativeErrorE")) {
			numberIterationsEM = dicalParams.jsapParams.getInt ("numberIterationsEM");
		}
		else {
			// should not be
			System.err.println ("Need to specify exactly one of relativeErrorE or numberIterationsEM.");
			System.exit(-1);
		}			

		// default is no parallel
		Integer parallelThreads = null;
		if (dicalParams.jsapParams.contains("parallel")) {
			parallelThreads = dicalParams.jsapParams.getInt("parallel");
			if (parallelThreads < 1) {
				System.err.println ("Need to a positive number of parallel threads (Not " + parallelThreads + ").");
				System.exit(-1);
			}
		}
		
		
		Double relativeErrorM = null;
		Integer numberIterationsMstep = null;
		if (dicalParams.jsapParams.contains("relativeErrorM") && !dicalParams.jsapParams.contains("numberIterationsMstep")) {
			relativeErrorM = dicalParams.jsapParams.getDouble ("relativeErrorM");
		}
		else if (dicalParams.jsapParams.contains("numberIterationsMstep") && !dicalParams.jsapParams.contains("relativeErrorM")) {
			numberIterationsMstep = dicalParams.jsapParams.getInt ("numberIterationsMstep");
		}
		else if (numberIterationsEM != 0) {
			// should not happen
			System.err.println ("Need to specify exactly one of relativeErrorM or numberIterationsMstep or 0 EmSteps.");
			System.exit(-1);
		}			

		
		
		boolean estimateRecomScaling = dicalParams.jsapParams.getBoolean("estimateRecomScaling");
		boolean estimateTheta= dicalParams.jsapParams.getBoolean("estimateTheta");
		
		
		if (estimateRecomScaling) {
			if (dicalParams.jsapParams.contains("bounds")) throw new IOException ("Bounds on command line incompatible with estimating reco scaling.");
		}
		if (estimateTheta) {
			if (dicalParams.jsapParams.contains("bounds")) throw new IOException ("Bounds on command line incompatible with estimating theta.");
			if (! dicalParams.jsapParams.contains("thetaIntervals")) throw new IOException("You must specify thetaIntervals when estimating theta");
			if (!dicalParams.jsapParams.getString("thetaIntervals").equals("customFixed")) throw new IOException("For now thetaIntervals must be customFixed");
			if(! dicalParams.jsapParams.contains("thetaIntervalParams")) throw new IOException("thetaIntervals need thetaIntervalParams");
			double lastPoint = 0d;
			String factoryParams = dicalParams.jsapParams.getString("thetaIntervalParams");
			String[] intervalFactoryParams = factoryParams.split(",");
			double[] points = new double[intervalFactoryParams.length];
			if(intervalFactoryParams.length < 1) throw new IOException ("You messed up the thetaIntervalParams");
			for (int i = 0; i < intervalFactoryParams.length; i++) {
				double nextPoint = Double.parseDouble(intervalFactoryParams[i]);
				if (nextPoint <= lastPoint) throw new IOException("Something in your thetaIntervalParams is messed up");
				points[i] = nextPoint;
				lastPoint = nextPoint;
			}
			thetaIntervalFactory = new IntervalFactory.CustomFixedIntervalFactory (points, false);
			otherIntervalFactory = dicalParams.intervalFactory;
		}
		
		
		
		//TODO: Add an assertion that makes sure the length of startPoint is the same as
		//the number of demography variables to be inferred
		String startPointString = dicalParams.jsapParams.getString("startPoint");
		double[] startPoint = null;
		if(startPointString == null){
			if (estimateTheta) startPoint = new double[0];
			// otherwise just leave it null
		} else {
			String[] startPointCoords = startPointString.split(",");
		
			startPoint = new double[startPointCoords.length];
			for (int i = 0; i < startPoint.length; i++){
				startPoint[i] = Double.parseDouble(startPointCoords[i]);
			}
			
			double[][] themBounds = dicalParams.demoFactory.getBounds();
			if ((themBounds != null) && (themBounds.length != startPoint.length)) throw new IOException("Dimensions of startPoint and bounds don't match");
		}
		assert ((startPoint != null) || (dicalParams.demoFactory.getBounds() != null));
		// check dimensions
		if ((startPoint != null) && (startPoint.length != dicalParams.demoFactory.getDimension())) {
			throw new IOException("Dimension of starting point does not match dimension of demoFile.");
		}
		if ((dicalParams.demoFactory.getBounds() != null) && (dicalParams.demoFactory.getBounds().length != dicalParams.demoFactory.getDimension())) {
			throw new IOException("Dimension of starting point does not match dimension of demoFile.");
		}
		
		
		// get the order to optimize coordinates in MStep
		int[] coordinateOrder = null;
		if (dicalParams.jsapParams.contains("coordinateOrder")) {
			if (!coordinateWiseM) throw new IOException("coordinateOrder should only be specified when using coordinateWiseMStep flag.");
			
			String[] coordinateOrderString = dicalParams.jsapParams.getString("coordinateOrder").split(",");
			
			coordinateOrder = new int[coordinateOrderString.length];
			for (int i = 0; i < coordinateOrder.length; i++){
				coordinateOrder[i] = Integer.parseInt(coordinateOrderString[i]);
			}

			// check that all coordinates are there
			boolean[] coordinatePresent = new boolean[coordinateOrder.length];
			for (int coord : coordinateOrder) {
				coordinatePresent[coord] = true;
			}
			
			boolean allPresent = true;
			for (boolean coordPresent : coordinatePresent) allPresent = allPresent && coordPresent;
			
			if (!allPresent || coordinateOrder.length != dicalParams.demoFactory.getDimension()) {
				throw new IOException("Must specify exactly one of each coordinate in coordinateOrder");
			}
		}

		if (estimateRecomScaling && dicalParams.conditionalObjectiveFunction == ConditionalObjectiveFunctionType.MarginalKLDivergence) {
			throw new IOException("Cannot estimate recom rate using marginalKL.");
		}
		
		boolean useParamRelErr = dicalParams.jsapParams.getBoolean("useParamRelErr");
		if (useParamRelErr && dicalParams.jsapParams.contains("numberIterationsEM")) {
			throw new IOException ("Cannot use --useParamRelErr flag with --numberIterationsEM flag");
		}
		boolean useParamRelErrM = dicalParams.jsapParams.getBoolean("useParamRelErrM");
		if (useParamRelErrM && dicalParams.jsapParams.contains("numberIterationsMstep")) {
			throw new IOException ("Cannot use --useParamRelErrM flag with --numberIterationsMstep flag");
		}

		// the composite likelihood is now in dical
		String compositeLikelihoodStringLow = dicalParams.compositeLikelihoodStringLow;
		CompositeLikelihoodType compositeLikelihood = dicalParams.compositeLikelihood;
		List<List<Integer>> loadedCsdList = dicalParams.loadedCsdList;
		
		// get the number of permutations
		int numChunks = dicalParams.extendedConfigInfoList.size();
		List<List<Permutation>> metaPermutationList = null;
		String permutationsFileString = null;
		if (dicalParams.jsapParams.contains("numPermutations")) {
			List<Permutation> permutationList = null;
			if ((compositeLikelihood != CompositeLikelihoodType.PAC) && (compositeLikelihood != CompositeLikelihoodType.OLD_PAC) && (compositeLikelihood != CompositeLikelihoodType.SuperPAC)) throw new IOException("numPermutations not allowed for composite likelihood: " + compositeLikelihoodStringLow);
			if (dicalParams.jsapParams.contains("permutationsFile")) throw new IOException("Specify either numPermutations or permutations files.");
			metaPermutationList = new ArrayList<List<Permutation>>();
			int numPermutations = dicalParams.jsapParams.getInt("numPermutations");
			// fill permutation list
			for (int chunk=0; chunk < numChunks; chunk++) {
				// new one every time, unless explicitly wanted otherwise
				if (permutationList == null || dicalParams.jsapParams.getBoolean("diffPermsPerChunk")) {
					int numHap = dicalParams.extendedConfigInfoList.get(chunk).structuredConfig.getTotalGeneticTypes();
					permutationList = Permutation.generatePermutationList (numPermutations, numHap, randomGenerator);
					metaPermutationList.add(permutationList);
				}
			}

		}
		else if (dicalParams.jsapParams.contains("permutationsFile")) {
			List<Permutation> permutationList = null;
			if ((compositeLikelihood != CompositeLikelihoodType.PAC) && (compositeLikelihood != CompositeLikelihoodType.OLD_PAC) && (compositeLikelihood != CompositeLikelihoodType.SuperPAC)) throw new IOException("permutationsFile not allowed for composite likelihood: " + compositeLikelihoodStringLow);
			metaPermutationList = new ArrayList<List<Permutation>>();
			// load the list of permutations from file
			permutationsFileString = dicalParams.jsapParams.getString("permutationsFile");
			String[] permutationsFile = permutationsFileString.split(",");
			if (permutationsFile.length == 1) {
				if ((numChunks != 1) && dicalParams.jsapParams.getBoolean("diffPermsPerChunk")) {
					throw new IOException("Not enough permutation files for the chunks given.");
				}
				FileReader permReader = new FileReader (permutationsFile[0]);
				permutationList = Permutation.readPermutationList (permReader);
				permReader.close();
				metaPermutationList.add(permutationList);
			}
			else {
				if (numChunks != permutationsFile.length) throw new IOException ("You have to give as many permutations as chunks.");
				for (int perm=0; perm < permutationsFile.length; perm++) {
					FileReader permReader = new FileReader (permutationsFile[perm]);
					permutationList = Permutation.readPermutationList (permReader);
					permReader.close();
					metaPermutationList.add(permutationList);
				}
			}
		}			
		else {
			if ((compositeLikelihood != CompositeLikelihoodType.PCL)
					&& (compositeLikelihood != CompositeLikelihoodType.LOL)
					&& (compositeLikelihood != CompositeLikelihoodType.PCLOL)
					&& (compositeLikelihood != CompositeLikelihoodType.FILE)) {
				throw new IOException("permutations needed for composite likelihood: " + compositeLikelihoodStringLow);
			}
		}
		// should be it
		

		// then initialize
		List<List<List<CSDConfig<HFSAXFullHaplotypeShell>>>> csdConfigs = new ArrayList<List<List<CSDConfig<HFSAXFullHaplotypeShell>>>>();
		for (int c=0; c<numChunks; c++) {
			List<List<CSDConfig<HFSAXFullHaplotypeShell>>> thisChunkCsdConfigs = new ArrayList<List<CSDConfig<HFSAXFullHaplotypeShell>>>();
			if ((compositeLikelihood != CompositeLikelihoodType.OLD_PAC) || (c == 0)) {
				csdConfigs.add(thisChunkCsdConfigs);
			}
			// deal with this chunk
			if ((compositeLikelihood == CompositeLikelihoodType.LOL)
					|| (compositeLikelihood == CompositeLikelihoodType.SuperPAC)
					|| (compositeLikelihood == CompositeLikelihoodType.PCL)
					|| (compositeLikelihood == CompositeLikelihoodType.PCLOL)
					|| (compositeLikelihood == CompositeLikelihoodType.FILE)) {
				thisChunkCsdConfigs.add(new ArrayList<CSDConfig<HFSAXFullHaplotypeShell>>());
			} else if (compositeLikelihood == CompositeLikelihoodType.PAC) {
				int realC = (metaPermutationList.size() <= 1 ? 0 : c);
				for (int perm = 0; perm < metaPermutationList.get(realC).size(); perm++) {
					thisChunkCsdConfigs.add(new ArrayList<CSDConfig<HFSAXFullHaplotypeShell>>());
				}
			} else if (compositeLikelihood == CompositeLikelihoodType.OLD_PAC) {
				if (c == 0) {
					for (int perm = 0; perm < metaPermutationList.get(c).size(); perm++) {
						thisChunkCsdConfigs.add(new ArrayList<CSDConfig<HFSAXFullHaplotypeShell>>());
					}
				}
			} 
		}
		// should be good
		assert (csdConfigs.size() == numChunks || compositeLikelihood == CompositeLikelihoodType.OLD_PAC);
		
		
		// get the trunk process factory
		TrunkProcessFactory trunkFactory = dicalParams.trunkFactory;
		
		// go through
		for (int chunk = 0; chunk < numChunks; chunk++) {
			
			DemoConfiguration<HFSAXFullHaplotypeShell> structuredConfig = dicalParams.extendedConfigInfoList.get(chunk).structuredConfig;
			List<DemeHapPair<HFSAXFullHaplotypeShell>> demeHaplotypeList = dicalParams.extendedConfigInfoList.get(chunk).demeHaplotypeList;
			EigenParamSet pSet = dicalParams.extendedConfigInfoList.get(chunk).pSet;
			HmmStepHandler<HFSAXFullHaplotypeShell> fancyTransitionMap = dicalParams.extendedConfigInfoList.get(chunk).fancyTransitionMap;

			
			// get the right list of permutations
			List<Permutation> permutationList = null;
			if (metaPermutationList != null) {
				if (metaPermutationList.size() == 1) {
					permutationList = metaPermutationList.get(0);
				}
				else {
					permutationList = metaPermutationList.get(chunk);
				}
				assert (permutationList != null);
			
				// print out permutations
				// this is not always necessary (for LOL and PCL)
				outStream.println("# CHUNK " + chunk + " PERMUTATIONS");
				ArrayList<SimplePair<Integer,Integer>> hapDemeIdxPairList = new ArrayList<SimplePair<Integer,Integer>>();
				for (int i = 0; i < demeHaplotypeList.size(); i++) hapDemeIdxPairList.add(new SimplePair<Integer, Integer>(i, demeHaplotypeList.get(i).demeNumber));
				for (Permutation p : permutationList) {
					ArrayList<SimplePair<Integer,Integer>> permutedHapDemePairs = p.permute(hapDemeIdxPairList);
					
					String hapOrder = "";
					String demeOrder = "";
					
					for (SimplePair<Integer,Integer> pair : permutedHapDemePairs) {
						hapOrder = hapOrder + pair.first() + "\t";
						demeOrder = demeOrder + pair.second() + "\t";
					}
					
					outStream.println("# HAP PERMUTATION:\t" + hapOrder);
					outStream.println("# DEME PERMUTATION:\t" + demeOrder + "\n#");
				}
			}
			else {
				assert ((compositeLikelihood == CompositeLikelihoodType.FILE)
						|| (compositeLikelihood == CompositeLikelihoodType.LOL)
						|| (compositeLikelihood == CompositeLikelihoodType.PCLOL)
						|| (compositeLikelihood == CompositeLikelihoodType.PCL));
			}
	

			// which composite likelihood
			if (compositeLikelihood == CompositeLikelihoodType.LOL){
				System.out.println ("# [WARNING] permutations will be ignored.");
				csdConfigs.get(chunk).get(0).addAll(getLOLConfigList (structuredConfig, pSet, fancyTransitionMap, false));
			} else if (compositeLikelihood == CompositeLikelihoodType.PCL) {
				System.out.println ("# [WARNING] permutations will be ignored.");
				csdConfigs.get(chunk).get(0).addAll(getPCLConfigList (structuredConfig, demeHaplotypeList, pSet, fancyTransitionMap, false));
			} else if (compositeLikelihood == CompositeLikelihoodType.PCLOL) {
				System.out.println ("# [WARNING] permutations will be ignored.");
				// add PCL and LOL =)
				csdConfigs.get(chunk).get(0).addAll(getPCLConfigList (structuredConfig, demeHaplotypeList, pSet, fancyTransitionMap, true));
				csdConfigs.get(chunk).get(0).addAll(getLOLConfigList (structuredConfig, pSet, fancyTransitionMap, true));
			} else if (compositeLikelihood == CompositeLikelihoodType.FILE) {
				System.out.println ("# [WARNING] permutations will be ignored.");
				csdConfigs.get(chunk).get(0).addAll(getFileConfigList (loadedCsdList, demeHaplotypeList, pSet, fancyTransitionMap, structuredConfig.getNumberDemes(), structuredConfig.numLoci(), structuredConfig.numAlleles()));
			} else if (compositeLikelihood == CompositeLikelihoodType.OLD_PAC
					|| compositeLikelihood == CompositeLikelihoodType.PAC
					|| compositeLikelihood == CompositeLikelihoodType.SuperPAC){
				
				// make a new visitor
				MakePACConfigList permVisitor = new MakePACConfigList (demeHaplotypeList, pSet, fancyTransitionMap, numCsdsPerPerm);
	
				// some numbers to permute
				ArrayList<Integer> listOfNumbers = new ArrayList<Integer>();
				for (int i=0; i< demeHaplotypeList.size(); i++) {
					listOfNumbers.add(i);
				}
				// visit the permutations
				VisitPermuations.visitPermutations (listOfNumbers , permVisitor, permutationList);
	
				List<List<CSDConfig<HFSAXFullHaplotypeShell>>> thisChunkConfigs = permVisitor.getCsdConfigs(); 
					
				assert (thisChunkConfigs.size() == permutationList.size());
				
				if (compositeLikelihood == CompositeLikelihoodType.OLD_PAC) {
					for (int perm = 0; perm < permutationList.size(); perm++) {
						csdConfigs.get(0).get(perm).addAll(thisChunkConfigs.get(perm));
					}
				} else if (compositeLikelihood == CompositeLikelihoodType.PAC) {
					// i think this is right
					assert (permutationList.size() == csdConfigs.get(chunk).size());
					for (int perm = 0; perm < permutationList.size(); perm++) {
						csdConfigs.get(chunk).get(perm).addAll(thisChunkConfigs.get(perm));
					}
				} else if (compositeLikelihood == CompositeLikelihoodType.SuperPAC) {
					for (List<CSDConfig<HFSAXFullHaplotypeShell>> thisPerm : thisChunkConfigs) {
						assert (csdConfigs.get(chunk).size() == 1);
						csdConfigs.get(chunk).get(0).addAll(thisPerm);
					}
				} else {
					assert (false);
				}
				
			} else {
				System.err.println("Invalid composite likelihood");
				System.exit(-1);
			}
		}	
		
		
		
		// this can be done cleaner
		for(int i= 0; i< dicalParams.extendedConfigInfoList.get(0).structuredConfig.numLoci(); i++){
			assert(csdConfigs.get(0).get(0).get(0).pSet.getMutationRate(0) == csdConfigs.get(0).get(0).get(0).pSet.getMutationRate(i));
		}
		int thetaDim;
		double[] thetaInitial;
		if(estimateTheta){
			thetaDim = dicalParams.jsapParams.getString("thetaIntervalParams").split(",").length + 1;
			thetaInitial= new double[thetaDim];
			for(int i=0; i< thetaDim; i++){
				thetaInitial[i]= csdConfigs.get(0).get(0).get(0).pSet.getMutationRate(0);
			}
		} else{
			thetaDim= 0;
			thetaInitial= new double[1];
			thetaInitial[0]= csdConfigs.get(0).get(0).get(0).pSet.getMutationRate(0);
		}
		
		String permutationString = null;
		if (metaPermutationList != null) {
			permutationString = "" + metaPermutationList.get(0).size();
			for (int i=1; i<metaPermutationList.size(); i++) {
				permutationString += ("," + metaPermutationList.get(i).size());
			}
		}
		
		// get them meta arguments
		metaArguments.get (dicalParams.jsapParams, dicalParams.demoFactory.getBounds());
		
		if (startPoint == null) {
			if ((metaArguments.numStartPoints == null) || (metaArguments.numStartPoints <= 1)) {
				if (metaArguments.metaStartFile == null) {
					throw new IOException("Provide either startpoint, or a number of start points to choose that is greater than 1, or a start point file.");
				}
			}
		}
		else {
			if (((metaArguments.numStartPoints != null) && (metaArguments.numStartPoints > 1)) || (metaArguments.metaStartFile != null)) {
				throw new IOException("Provide either startpoint, or a number of start points to choose that is greater than 1, or a start point file.");
			}
		}
		if ((metaArguments.numStartPoints != null) && (metaArguments.numStartPoints > 1) && (metaArguments.startingRanges == null)) {
			throw new IOException ("Need starting ranges to construct the starting points");
		}

		if (seeed != null) {
			outStream.println("# seed = " + seeed.longValue());
		}
		outStream.println("# numPermutations = " + permutationString);
		outStream.println("# startPoint = " + startPointString);
		outStream.println("# compositeLikelihoods = " + compositeLikelihoodStringLow);
		outStream.println("# useParamRelErr = " + useParamRelErr);
		outStream.println("# useParamRelErrM = " + useParamRelErrM);
		outStream.println ("# EM threshold = " + relativeErrorE);
		outStream.println ("# EM number iterations = " + numberIterationsEM);
		outStream.println ("# M-Step threshold = " + relativeErrorM);
		outStream.println ("# M-Step number iterations = " + numberIterationsMstep);
		outStream.println("# Estimate recom scaling = " + estimateRecomScaling);
		outStream.println("# Estimate theta = " + estimateTheta);
		outStream.println("# Coordinatewise M Step = " + coordinateWiseM);
		outStream.println("# Coordinate order for optimization = " + coordinateOrder);
		outStream.println("# parallel threads = " + parallelThreads);
		outStream.println("# permutations file = " + permutationsFileString);
		outStream.println("# number CSDs per permutation = " + numCsdsPerPerm);
		
		
		// also print the meta arguments
		metaArguments.print (outStream);

		
		// call the meta wrapper
		metaEmWraper (csdConfigs, startPoint, estimateRecomScaling, dicalParams.demoFactory, dicalParams.demoStateFactory, dicalParams.conditionalObjectiveFunction, trunkFactory, printFullEmPath, relativeErrorE, numberIterationsEM, relativeErrorM, numberIterationsMstep, useParamRelErr, useParamRelErrM, dicalParams.nelderMeadFraction, coordinateWiseM, coordinateOrder, parallelThreads, randomGenerator, estimateTheta, thetaInitial, metaArguments, verbose, dicalParams.useEigenCore);
		
		// final timing
		long endTime = System.currentTimeMillis();
		System.out.println ("# Total time elapsed (for whole execution of the program): " + (endTime - programStartTime) + "    (" + String.format("%.3f", (endTime - programStartTime)/1000d/60d/60d) + " hrs)");
		
	}

	
	private static void metaEmWraper (List<List<List<CSDConfig<HFSAXFullHaplotypeShell>>>> csdConfigs, double[] metaStartingPoint, boolean estimateRecomScaling, DemographyFactory demoFactory, DemoStateFactory demoStateFactory, ConditionalObjectiveFunctionType objectiveType, TrunkProcessFactory trunkFactory, boolean printFullEmPath, Double relativeErrorE, Integer numberIterationsEM, Double relativeErrorM, Integer numMstepIters, boolean useParamRelError, boolean useParamRelErrorM, Double nelderMeadFraction, boolean coordinateWiseM, int[] coordinateOrder, Integer parallelThreads, DelayedRandom masterRandom, boolean estimateTheta, double[] thetaInitial, MetaOptimizationArguments metaArguments, boolean verbose, boolean useEigenCore) throws IOException {
		
		long startTime = System.currentTimeMillis();
		
		// some things that we might need for theta
		int thetaDim = (estimateTheta ? thetaInitial.length : 0);
		double defaultMutRate = thetaInitial[0];
		if (thetaDim == 0) {
			assert (!estimateTheta);
			assert (thetaInitial.length == 1);
			assert (defaultMutRate == csdConfigs.get(0).get(0).get(0).pSet.getMutationRate(0));
		}
		
		
		// get (a bunch) of starting points
		List<double[]> startingPoints = new ArrayList<double[]>();
		if (((metaArguments.numStartPoints == null) || (metaArguments.numStartPoints == 1)) && (metaArguments.metaStartFile == null)) {
			// only a single point
			assert (metaStartingPoint != null);
			assert (metaArguments.numMetaIterations == 1);
			double[] currPoint = MetaOptimization.addRecoMutEstimation (metaStartingPoint, estimateRecomScaling, thetaDim, thetaInitial);
			startingPoints.add (currPoint);			
		}
		else {
			assert (metaStartingPoint == null);
			if (metaArguments.gridStart) {
				// build a grid
				startingPoints = MetaOptimization.buildGrid (metaArguments.startingRanges, metaArguments.numStartPoints);
			}
			else {
				if (metaArguments.metaStartFile != null) {
					// read start points from file
					startingPoints = MetaOptimization.readPointsFromFile (metaArguments.metaStartFile, demoFactory.getDimension());
				}
				else {
					// draw some random points
					startingPoints = MetaOptimization.randomPoints (metaArguments.startingRanges, metaArguments.numStartPoints, demoFactory, estimateRecomScaling, estimateTheta, thetaDim, masterRandom);
				}
			}
		}
		
		
		List<double[]> currStartingPoints = startingPoints;
		List<EmResult> finalResults = null;
		// additional steps?
		for (int currMetaStep = 0; currMetaStep<metaArguments.numMetaIterations; currMetaStep++) {

			System.out.println ("# [META_STEP] " + currMetaStep);
		
			
			// parallel, if we want to
			ExecutorService taskExecutor = null;
			if (parallelThreads != null) {
				taskExecutor = new ForkJoinPool (parallelThreads);
			}
			
			
			// run EMs for all the start points
			List<EmResult> results = new ArrayList<EmResult>();
			int currStartIdx = 0;
			
			while (currStartIdx < currStartingPoints.size()) {
								
				List<Future<EmResult>> currBatch = new ArrayList<Future<EmResult>>();
				// start them in batches, so we can control the memory consumption somehow
				for ( ; (currStartIdx < currStartingPoints.size()) && (currBatch.size() < metaArguments.metaParallelEmSteps); currStartIdx++) {
					
					// some randomness
					DelayedRandom slaveRandom = masterRandom.spawnOffspring();

					// the coordinate order and also the direction of first step is chosen randomly each time,
					// thus something can change after being stuck
				
					// get a thread for this one
					EmThread thread = new EmThread (csdConfigs, currStartingPoints.get(currStartIdx), demoFactory, demoStateFactory, objectiveType, trunkFactory, printFullEmPath, relativeErrorE, numberIterationsEM, relativeErrorM, numMstepIters, useParamRelError, useParamRelErrorM, nelderMeadFraction, coordinateWiseM, coordinateOrder, estimateRecomScaling, estimateTheta, thetaDim, defaultMutRate, currMetaStep, currStartIdx, slaveRandom, verbose, useEigenCore);
				
					// and remember some future task for this
					Future<EmResult> theFuture = null;
					if (taskExecutor != null) {
						// this is parallel
						theFuture = taskExecutor.submit (thread);
					}
					else {
						FutureTask<EmResult> futureTask = new FutureTask<EmResult> (thread);
						futureTask.run();
						theFuture = futureTask;
					}
					
					// add it
					assert (theFuture != null);
					currBatch.add (theFuture);
				}
				
				
				// then collect for this batch
				for (Future<EmResult> f : currBatch) {
					// try to get the result
					EmResult result = null;
					try {
						// get it
						result = f.get();
					} catch (InterruptedException e) {
						System.err.println("Interrupted excpetion in em thread:");
						e.getCause().printStackTrace(System.err); // this is more informative than printing the stack trace from e
						System.exit(-1);
					} catch (ExecutionException e) {
						System.err.println("execution exception in em thread:");
						e.getCause().printStackTrace(System.err); // this is more informative than printing the stack trace from e
						System.exit(-1);
					}
					assert (result != null);
					// and add it to total list
					results.add (result);
				}
				
				// and do it again, if there is more
			}
			
			
			// if we are at last meta step, we stop here
			if (currMetaStep >= metaArguments.numMetaIterations - 1) {
				finalResults = results;
				break;
			}
		
			// in case we want more meta steps, resample the next points from you current set
			currStartingPoints = MetaOptimization.nextGeneration (results, metaArguments.metaKeepBest, metaArguments.metaNumPoints, metaArguments.stretchProportion, metaArguments.disperseFactor, metaArguments.sdPercentage, demoFactory, estimateRecomScaling, estimateTheta, thetaDim, masterRandom);
			// that's it, go to next round
		}
		
		
		// get the time at the end
		long endTime = System.currentTimeMillis();
		long thisTime = endTime - startTime;
		
		// print the results
		System.out.println ("# ++++++++++ final results ++++++++++");
		
		if (!printFullEmPath) {
			for (int i=0; i<finalResults.size(); i++) {
				printPoint ("final_" + i, finalResults.get(i), thisTime);
			}
		}
		else {
			// only time
			System.out.println ("# Total time elapsed (for the EM): " + thisTime);
		}
	}

	private static void printPoint (String idString, EmResult result, long time) {
		// header
		String debugString = "";
		debugString += "# LogLikelihood\tTime";
		for (int i = 0; i < result.mle.length; i++) debugString += "\tcoord" + i;
		debugString += "\t[idString]";		
		debugString += "\n";
		
		// print out result
		debugString += "" + result.logLikelihood + "\t" + time;
		for (int i = 0; i < result.mle.length; i++) debugString += "\t" + result.mle[i];
		debugString += "\t["+ idString + "]";
		
		MetaOptimization.synchronizedPrintln (debugString);
	}

	public static class EmResult {
		public final double[] mle;
		public final double logLikelihood;

		public EmResult (double[] mle, double logLikelihood) {
			this.mle = mle;
			this.logLikelihood = logLikelihood;
		}
	}

	public static class EmThread implements Callable<EmResult>{
	
		// some parameters
		final private List<List<List<CSDConfig<HFSAXFullHaplotypeShell>>>> csdConfigs;
		final private double[] startingPoint;
		final private DemographyFactory demoFactory;
		final private DemoStateFactory demoStateFactory;
		final private ConditionalObjectiveFunctionType objectiveType;
		final private TrunkProcessFactory trunkFactory;
		final private boolean printFullEmPath;
		final private Double relativeErrorE;
		final private Integer numberIterationsEM;
		final private Double relativeErrorM;
		final private Integer numMstepIters;
		final private boolean useParamRelError;
		final private boolean useParamRelErrorM;
		final private Double nelderMeadFraction;
		final private boolean coordinateWiseM;
		final private int[] coordinateOrder;
		final private boolean estimateRecomScaling;
		final private boolean estimateTheta;
		final private int thetaDim;
		final private double defaultMutRate;
		final private int currMetaStep;
		final private int currIdx;
		final private DelayedRandom rGen;
		final private boolean verbose;
		final private boolean useEigenCore;
		

		public EmThread (List<List<List<CSDConfig<HFSAXFullHaplotypeShell>>>> csdConfigs, double[] startingPoint, DemographyFactory demoFactory, DemoStateFactory demoStateFactory, ConditionalObjectiveFunctionType objectiveType, TrunkProcessFactory trunkFactory, boolean printFullEmPath, Double relativeErrorE, Integer numberIterationsEM, Double relativeErrorM, Integer numMstepIters, boolean useParamRelError, boolean useParamRelErrorM, Double nelderMeadFraction, boolean coordinateWiseM, int[] coordinateOrder, boolean estimateRecomScaling, boolean estimateTheta, int thetaDim, double defaultMutRate, int currMetaStep, int currIdx, DelayedRandom rGen, boolean verbose, boolean useEigenCore) {
			super();
			this.csdConfigs = csdConfigs;
			this.startingPoint = startingPoint;
			this.demoFactory = demoFactory;
			this.demoStateFactory = demoStateFactory;
			this.objectiveType = objectiveType;
			this.trunkFactory = trunkFactory;
			this.printFullEmPath = printFullEmPath;
			this.relativeErrorE = relativeErrorE;
			this.numberIterationsEM = numberIterationsEM;
			this.relativeErrorM = relativeErrorM;
			this.numMstepIters = numMstepIters;
			this.useParamRelError = useParamRelError;
			this.useParamRelErrorM = useParamRelErrorM;
			this.nelderMeadFraction = nelderMeadFraction;
			this.coordinateWiseM = coordinateWiseM;
			this.coordinateOrder = coordinateOrder;
			this.estimateRecomScaling = estimateRecomScaling;
			this.estimateTheta = estimateTheta;
			this.thetaDim = thetaDim;
			this.defaultMutRate = defaultMutRate;
			this.currMetaStep = currMetaStep;
			this.currIdx = currIdx;
			this.rGen = rGen;
			this.verbose = verbose;
			
			this.useEigenCore = useEigenCore;
		}
		
		@Override
		public EmResult call() throws Exception {
			
			// some locals
			double[] currPoint = this.startingPoint;
			double[] prevPoint = null;
			double currLogLike = Double.NEGATIVE_INFINITY;
			double prevLogLike = Double.NEGATIVE_INFINITY;
			int numIterations = 0;
			
			// do loop
			while (true) {
				
				// last iteration?
				boolean lastIteration = false;
				if ((this.numberIterationsEM != null) && (this.numberIterationsEM == numIterations)) {
					// we are at last
					lastIteration = true;
				}
				
				// remember previous likelihood (deep)
				prevLogLike = currLogLike;	

				// E step
				String idString = this.currMetaStep + "_" + numIterations + "_" + this.currIdx;
				MetaOptimization.synchronizedPrintln ("# [EXPECTATION_" + idString + "]\t" + Arrays.toString(currPoint));

				// get time
				long eStartTime = System.currentTimeMillis();

				// get the new likelihoods
				// that is, add them to the executor
				CompositeLikelihoodObjectiveFunction objectiveFunction = null;
				// valid start point?
				if ((currPoint != null) && (this.demoFactory.getDemography (DiCalObjectiveFunction.getDemoPointStatic (currPoint, this.estimateRecomScaling, this.estimateTheta, this.thetaDim)) != null)) {
					objectiveFunction = new DiCalObjectiveFunction (this.csdConfigs, this.demoFactory, this.demoStateFactory, this.objectiveType, currPoint, this.estimateRecomScaling, this.trunkFactory, this.estimateTheta, this.defaultMutRate, this.thetaDim, lastIteration, this.verbose, this.useEigenCore);
				}

				// get the likelihood (set to -infty if no proper objective function was calculated)
				currLogLike = Double.NEGATIVE_INFINITY;
				if (objectiveFunction != null) {
					currLogLike = objectiveFunction.getLogLikelihood();
				}

				// timing
				long eEndTime = System.currentTimeMillis();
				long eTime = eEndTime - eStartTime;

				EmResult currEmResult = new EmResult (currPoint, currLogLike);

				assert (!Double.isNaN(currEmResult.logLikelihood));
				// we have to be somewhat lenient here since we switched to storing only floats
				// assert (prevLogLike <= currEmResult.logLikelihood);
				if (prevLogLike > currEmResult.logLikelihood) {
					assert (prevLogLike/currEmResult.logLikelihood < 1+5e-8);
				}

				// what to print
				if (!this.printFullEmPath) {
					// just some commented things
					String debugString = "# [LOGLIKE_" + idString + "] at " + Arrays.toString(currEmResult.mle) + ": " + currLogLike + "\n";
					// print time
					debugString += "# [EXPECT_TIME_" + idString + "] " + eTime;
					MetaOptimization.synchronizedPrintln (debugString);
				}
				else {
					// print it nicely for R
					printPoint (idString, currEmResult, eTime);
				}
			
				// do we want to stop?
				// stopping criterion
				boolean stop = stopEm (this.numberIterationsEM, this.relativeErrorE, this.useParamRelError, numIterations, currLogLike, currPoint, prevLogLike, prevPoint);
								
				// only go on if there is more to do
				if (stop) break;
				
				MetaOptimization.synchronizedPrintln ("# [MAXIMIZATION_" + idString + "]");
				
				// get a new point
				prevPoint = Arrays.copyOf(currPoint, currPoint.length);
				if (!this.coordinateWiseM) {
					currPoint = updatePoint (currPoint, objectiveFunction, this.numMstepIters, this.relativeErrorM, this.nelderMeadFraction, this.useParamRelErrorM, this.rGen, verbose);
				}
				else {
					currPoint = updatePointCoordinatewise  (currPoint, objectiveFunction, this.numMstepIters, this.relativeErrorM, this.nelderMeadFraction, this.useParamRelErrorM, this.coordinateOrder, this.rGen, verbose);
				}

				// and round it goes
				numIterations++;

				// print time
				double MEndTime = System.currentTimeMillis();
				String debugString = "# [MAX_TIME_" + idString + "] " + (MEndTime - eEndTime) + "\n";
				// and iteration
				debugString += "# [ITERATION_" + idString + "] " + numIterations + " iterations finished.";
				MetaOptimization.synchronizedPrintln (debugString);
			}
			
			return new EmResult (currPoint, currLogLike);
		}
	}
	
	private static Boolean stopEm (Integer numberIterationsEM, Double relativeErrorE, boolean useParamRelError, int numIterations, double currLogLike, double[] currPoint, double prevLogLike, double[] prevPoint) {
		
		// this is the stopping criterion
		if (currLogLike == Double.NEGATIVE_INFINITY) {
			return true;
		}
				
		// if number iteration is stopping criterion, and the number is used up, then stop
		if ((numberIterationsEM != null) && (numIterations >= numberIterationsEM)) {
			return true;
		}
	
		// if relative improvement is stopping criterion, and if there is no more improvement, then stop
		if ((relativeErrorE != null) && (prevPoint != null)) {
			
			double relErr;
			
			if (!useParamRelError) {
				// use the relative error of the likelihood
				relErr = relativeError (currLogLike, prevLogLike, UberDemographyCore.ONELOCUS_EPSILON);
			} else {
				// use the relative error of the parameters
				relErr = maxRelativeError (currPoint, prevPoint, UberDemographyCore.ONELOCUS_EPSILON);
			}
			
			System.out.println ("# [RELATIVE EM ERROR] current: " + relErr + "\tgoal: " + relativeErrorE);

			if (relErr < relativeErrorE) {
				// stop
				return true;
			}
		}
		return false;
	}

	private static double relativeError (double valueOne, double valueTwo, double EPSILON) {
		// get some nice relative error with the right sign
		double relErr = Math.abs (valueOne - valueTwo);
		if (relErr < EPSILON && Math.abs(valueTwo) < EPSILON) {
			return 0d;
		}
		else {
			relErr /= Math.abs (valueTwo);
		}
		// return it
		return relErr;
	}

	private static double maxRelativeError (double[] vectorOne, double[] vectorTwo, double EPSILON) {
		assert (vectorOne.length == vectorTwo.length);
		
		double maxRelativeError = 0d;
		for (int i = 0; i < vectorOne.length; i++) {
			maxRelativeError = Math.max(maxRelativeError, relativeError(vectorOne[i], vectorTwo[i], EPSILON));
		}
		
		return maxRelativeError;
	}
	
	private static double[] getInitialStepSizes (double[] currPoint, Double fraction, DelayedRandom rand) {
		// first get the dimension
		int dim = currPoint.length;
		
		// then something to return
		double[] toReturn = new double[dim];
		// now go through
		for (int i=0; i<dim; i++) {
			// store some percentage as the stepsize or the minimum
			toReturn[i] = Math.max (NM_MIMUM_STEPSIZE, (fraction * currPoint[i]) * (rand.nextBoolean() ? 1d : -1d));
		}
		// should be good to go		
		return toReturn;
	}

	private static class ClampedObjectiveFunction implements MultivariateFunction {
		final CompositeLikelihoodObjectiveFunction unclamped;
		final double[] startPoint;
		final int unclampedIdx;
		
		public ClampedObjectiveFunction(
				CompositeLikelihoodObjectiveFunction unclamped,
				double[] startPoint, int unclampedIdx) {
			super();
			this.unclamped = unclamped;
			this.startPoint = startPoint;
			this.unclampedIdx = unclampedIdx;
		}

		@Override
		public double value(double[] arg0) {
			assert arg0.length == 1;
			
			double[] pointToEval = Arrays.copyOf(startPoint, startPoint.length);
			pointToEval[unclampedIdx] = arg0[0];
			
			return unclamped.value(pointToEval);
		}
		
	}
	
	
	private static double[] updatePointCoordinatewise (double[] currPoint, CompositeLikelihoodObjectiveFunction objective, Integer numMstepIters, Double relativeErrorM, Double nelderMeadFraction, boolean useParamRelErrorM, int[] coordinateOrder, DelayedRandom rand, boolean verbose) {

		double[] localCurrPoint = Arrays.copyOf(currPoint, currPoint.length);
		
		List<Integer> perm = new ArrayList<Integer>();
		if (coordinateOrder == null) {
			for (int i=0; i<localCurrPoint.length; i++) {
				perm.add(i);
			}
			Collections.shuffle (perm, rand.getInternalRandom());
		} else {
			for (int coord : coordinateOrder) {
				perm.add(coord);
			}
		}
		
		for (int i : perm) {
			
			// make a clamped one
			ClampedObjectiveFunction clam = new ClampedObjectiveFunction (objective, localCurrPoint, i);
			localCurrPoint[i] = updatePoint(new double[] {localCurrPoint[i]}, clam, numMstepIters, relativeErrorM, nelderMeadFraction, useParamRelErrorM, rand, verbose)[0];
			
		}
	
		// should be fine
		return localCurrPoint;	
	}


	
	private static double[] updatePoint (double[] currPoint, MultivariateFunction objectiveFunction, Integer numMstepIters, Double relativeMError, Double nelderMeadFraction, boolean useParamRelErr, DelayedRandom rand, boolean verbose) {
		// maximize the objective function
		
		// get a Nelder Mead object (for a simplex with dimension that of currPoint)
		NelderMeadSimplex nm = null;
		if (nelderMeadFraction == null) {
			// no special initial stepsize
			nm = new NelderMeadSimplex (currPoint.length);
		}
		else {
			// initial stepsize should be percentage of current point
			nm = new NelderMeadSimplex (getInitialStepSizes (currPoint, nelderMeadFraction, rand));
		}
		
		// initialize properly
		// for now start at currPoint
		nm.build (currPoint);
		
		// and evaluate the first simplex
		nm.evaluate (objectiveFunction, new PointValueCompare());
		
		int iterations = 0;
		while (true) {
			
			// see how good we are
			// get best and worst Q value on the simplex [assume them ordered]
			double bestQ = nm.getPoint(0).getValue();
			double worstQ = nm.getPoint(nm.getSize()-1).getValue();
			
			// find maxima and minima in the coordinates of the simplex
			double[] maximums = Arrays.copyOf (nm.getPoint(0).getPoint(), nm.getPoint(0).getPoint().length);
			double[] minimums = Arrays.copyOf (maximums, maximums.length);
			for (PointValuePair pv : nm.getPoints()) {
				for (int i=0; i<maximums.length; i++) {
					maximums[i] = Math.max (maximums[i], pv.getPoint()[i]);
					minimums[i] = Math.min (minimums[i], pv.getPoint()[i]);
				}
			}
			// the maximal error
			double errorCoord = maxRelativeError (maximums, minimums, UberDemographyCore.ONELOCUS_EPSILON);

			double errorQ = relativeError (bestQ, worstQ, UberDemographyCore.ONELOCUS_EPSILON);
			
			// the relative error between those
			if (relativeMError != null) {
				assert (numMstepIters == null);
				if (iterations > 0) {
					if (!useParamRelErr) {
						if (errorQ < relativeMError) {
							// stop it
							break;
						}
					} else {
						if (errorCoord < relativeMError) {
							break;
						}
					}
				}
			} else {
				assert numMstepIters != null;
				if (iterations >= numMstepIters) break;
			}
			
			// otherwise, go on
			if (verbose) {
				MetaOptimization.synchronizedPrintln ("# [NELDER_MEAD_ITERATION_" + iterations + "]");
			}
			
			// make step
			nm.iterate (objectiveFunction, new PointValueCompare());
			
			// increase iterations
			iterations++;

		}

		if (verbose) {
			// put some debug string together
			String debugString = "# [LAST_SIMPLEX] " + iterations + " iterations\n";
			for (PointValuePair pv : nm.getPoints()) {
				String tmpString = "# ";
				for (double x :pv.getPoint()) {
					tmpString += x + "\t";
				}
				//the (-1) gets the sign right on the likelihood.
				tmpString += (-1) * pv.getValue() + "\n";
				debugString += tmpString;
			}
			debugString += "# [MAX_AT] " + Arrays.toString(nm.getPoint(0).getPoint());
			// and print it threadsafe
			MetaOptimization.synchronizedPrintln(debugString);
		}
		
		// return the maximal point
		return nm.getPoint(0).getPoint();
	}


	private static SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> getObjectiveFunction (ConditionalObjectiveFunctionType objectiveType, SMCSDemo<HFSAXFullHaplotypeShell> csd, HFSAXFullHaplotypeShell additionalHap, boolean likelihoodOnly){
		// if the conditional config is empty, we have to get a special objective function
		if (csd.core.isConditionalConfigEmpty()) {
			return new EmptyConditionalConfigOF<HFSAXFullHaplotypeShell> (csd, additionalHap);
		}
		
		// maybe likelihood only
		if (likelihoodOnly) {
			return new SMCSDemoEMObjectiveFunction.LikelihoodOnlyObjectiveFunction<HFSAXFullHaplotypeShell> (csd, additionalHap);
		}
		
		// otherwise we have to switch it
		switch (objectiveType) {
			case ConditionLineageTransitionType: return new SMCSDemoEMObjectiveFunction.ConditionLineageTransitionTypeOF<HFSAXFullHaplotypeShell> (csd, additionalHap);
			case ConditionLineage: return new SMCSDemoEMObjectiveFunction.ConditionLineageOF<HFSAXFullHaplotypeShell> (csd, additionalHap);
			case MarginalKLDivergence: return new SMCSDemoEMObjectiveFunction.MarginalKLDivergence<HFSAXFullHaplotypeShell>(csd, additionalHap);
			default: assert(false); return null;
		}
	}
	
	public static class CSDConfig<H extends FSAHaplotype> {
		public final DemoConfiguration<H> trunkConfig;
		public final DemeHapPair<H> additionalHapDeme;
		public final HmmStepHandler<H> fancyTransitionMap;
		public final EigenParamSet pSet;
		private Double logFactor;
		
		public CSDConfig(DemoConfiguration<H> trunkConfig, DemeHapPair<H> additionalHapDeme, EigenParamSet pSet, HmmStepHandler<H> fancyTransitionMap) {
			this.trunkConfig = trunkConfig;
			this.additionalHapDeme = additionalHapDeme;
			this.pSet = pSet;
			this.fancyTransitionMap = fancyTransitionMap;
		}

		public Double getLogFactor() {
			// TODO Auto-generated method stub
			return this.logFactor;
		}
		
		public void setLogFactor (Double f) {
			this.logFactor = f;
		}
	}
	
	public static List<CSDConfig<HFSAXFullHaplotypeShell>> getLOLConfigList(DemoConfiguration<HFSAXFullHaplotypeShell> structuredConfig, EigenParamSet pSet, HmmStepHandler<HFSAXFullHaplotypeShell> fancyTransitionMap, boolean pcLol) {		
		List<CSDConfig<HFSAXFullHaplotypeShell>> toReturn = new ArrayList<CSDConfig<HFSAXFullHaplotypeShell>>();
		
		// also copy the configuration
		DemoConfiguration<HFSAXFullHaplotypeShell> workingConfig = DemoConfiguration.copyConfiguration (structuredConfig);		
		
		// get the expectations
		// loop over demes
		for (int d = 0; d < structuredConfig.getNumberDemes(); d++) {
			// then loop over that deme
			for (GeneticTypeMultiplicity<HFSAXFullHaplotypeShell> popEntry : structuredConfig.getPopulation(d)) {				
				// now take it out
				workingConfig.getPopulation(d).adjustType(popEntry.geneticType, -1);
				
				// and make a clone
				DemoConfiguration<HFSAXFullHaplotypeShell> thisConfig = DemoConfiguration.copyConfiguration (workingConfig);
				
				for (int i = 0; i < popEntry.multiplicity; i++) {
					toReturn.add(new CSDConfig<HFSAXFullHaplotypeShell>(thisConfig, new DemeHapPair<HFSAXFullHaplotypeShell>(d, popEntry.geneticType), pSet, fancyTransitionMap));
				}
				
				// put the guy back in
				workingConfig.getPopulation(d).adjustType(popEntry.geneticType, 1);
			}
		}

		if (pcLol) {
			for (int i=0; i<toReturn.size(); i++) {
				toReturn.get(i).setLogFactor(1d/toReturn.size());
			}
		}
		
		return toReturn;
	}
	

	public static List<CSDConfig<HFSAXFullHaplotypeShell>> getPCLConfigList (DemoConfiguration<HFSAXFullHaplotypeShell> structuredConfig, List<DemeHapPair<HFSAXFullHaplotypeShell>> demeHaplotypeList, EigenParamSet pSet, HmmStepHandler<HFSAXFullHaplotypeShell> fancyTransitionMap, boolean pcLol) {
		// return thing
		List<CSDConfig<HFSAXFullHaplotypeShell>> toReturn = new ArrayList<CSDConfig<HFSAXFullHaplotypeShell>>();
		
		// now all pairs
		// unordered, because for assymetric migration it matters, i think
		// empirically it does matter (a bit)
		for (int i=0; i<demeHaplotypeList.size(); i++) {
			for (int j=0; j<demeHaplotypeList.size(); j++) {
				// that wouldn't make sense
				if (i != j) {
					
					// now make a trunk configuration with one i in there
					DemoConfiguration<HFSAXFullHaplotypeShell> singleTrunkConfig = new  DemoConfiguration<HFSAXFullHaplotypeShell>(structuredConfig.getNumberDemes(), structuredConfig.numLoci(), structuredConfig.numAlleles());
					singleTrunkConfig.adjustHaplotypeMultiplicity (demeHaplotypeList.get(i).demeNumber, demeHaplotypeList.get(i).haplotype, 1);
					
					// then j is the additional one
					DemeHapPair<HFSAXFullHaplotypeShell> addPair = new DemeHapPair<HFSAXFullHaplotypeShell>(demeHaplotypeList.get(j).demeNumber, demeHaplotypeList.get(j).haplotype);
					toReturn.add(new CSDConfig<HFSAXFullHaplotypeShell>(singleTrunkConfig, addPair, pSet, fancyTransitionMap));
				}
			}
		}

		if (pcLol) {
			for (int i=0; i<toReturn.size(); i++) {
				toReturn.get(i).setLogFactor(1d/toReturn.size());
			}
		}

		return toReturn;
	}

	private static List<CSDConfig<HFSAXFullHaplotypeShell>> getFileConfigList (List<List<Integer>> loadedCsdList, List<DemeHapPair<HFSAXFullHaplotypeShell>> demeHaplotypeList, EigenParamSet pSet, HmmStepHandler<HFSAXFullHaplotypeShell> fancyTransitionMap, int numDemes, int numLoci, int numAlleles) {
		// make a list
		List<CSDConfig<HFSAXFullHaplotypeShell>> toReturn = new ArrayList<CSDConfig<HFSAXFullHaplotypeShell>>();
		// go through list
		for (List<Integer> list : loadedCsdList) {
			// first is the additional
			DemeHapPair<HFSAXFullHaplotypeShell> additionalHapDeme = demeHaplotypeList.get(list.get(0));
			// then we need the trunk config
			List<DemeHapPair<HFSAXFullHaplotypeShell>> subList = new ArrayList<DiCalParamSet.DemeHapPair<HFSAXFullHaplotypeShell>>();
			for (int i=1; i<list.size(); i++) {
				subList.add (demeHaplotypeList.get(list.get(i)));
			}
			DemoConfiguration<HFSAXFullHaplotypeShell> trunkConfig = new DemoConfiguration<HFSAXFullHaplotypeShell> (numDemes, numLoci, numAlleles, subList);
			// make config
			CSDConfig<HFSAXFullHaplotypeShell> config = new CSDConfig<HFSAXFullHaplotypeShell>(trunkConfig, additionalHapDeme, pSet, fancyTransitionMap);
			// and it it to be returned
			toReturn.add(config);
		}
		// return it
		return toReturn;
	}
	
	public static class MakePACConfigList implements PermutationVisitor<Integer> {
		private final List<List<CSDConfig<HFSAXFullHaplotypeShell>>> csdConfigs;
		private final List<DemeHapPair<HFSAXFullHaplotypeShell>> demeHaplotpyeList;
		private final EigenParamSet pSet;
		private final HmmStepHandler<HFSAXFullHaplotypeShell> fancyTransitionMap;
		private final int numDemes;
		private final int numLoci;
		private final int numAlleles;
		private final Integer numCsdsPerPerm;
		
		public MakePACConfigList (List<DemeHapPair<HFSAXFullHaplotypeShell>> demeHaplotpyeList, EigenParamSet pSet, HmmStepHandler<HFSAXFullHaplotypeShell> fancyTransitionMap, Integer numCsdsPerPerm) {
			super();
			this.demeHaplotpyeList = demeHaplotpyeList;
			this.csdConfigs = new ArrayList<List<CSDConfig<HFSAXFullHaplotypeShell>>>();
			this.pSet = pSet;
			this.fancyTransitionMap = fancyTransitionMap;
			// get the number of demes
			TIntSet demeSet = new TIntHashSet();
			for (DemeHapPair<HFSAXFullHaplotypeShell> pair : demeHaplotpyeList) {
				demeSet.add(pair.demeNumber);
			}
			this.numDemes = demeSet.size();
			// get the number of loci
			this.numLoci = pSet.numLoci();
			// get the number of alleles
			this.numAlleles = pSet.getMutationMatrix(0).getColumnDimension();
			this.numCsdsPerPerm = numCsdsPerPerm;
		}

		@Override
		public boolean visitPermutation(List<Integer> permutation) {
			
			// before messing around with the configuration, clone it
			DemoConfiguration<HFSAXFullHaplotypeShell> workingConfiguration = new DemoConfiguration<HFSAXFullHaplotypeShell>(this.numDemes, this.numLoci, this.numAlleles);
			
			// for this permutation
			List<CSDConfig<HFSAXFullHaplotypeShell>> thisPerm = new ArrayList<CSDConfig<HFSAXFullHaplotypeShell>>();
			
			
			// make a set of all the trunk sizes that should be used
			int totalNumHaps = this.demeHaplotpyeList.size();
			Set<Integer> trunkSizes = new TreeSet<Integer> ();
			if ((this.numCsdsPerPerm == null) || (this.numCsdsPerPerm >= totalNumHaps)) {
				// just all of them
				for (int size=0; size < totalNumHaps; size++) {
					trunkSizes.add (size);
				}
			}
			else if (this.numCsdsPerPerm == 1) {
				// only use the biggest one
				trunkSizes.add (totalNumHaps-1);
			}
			else if (this.numCsdsPerPerm == 2) {
				// biggest and smallest that is not empty trunk
				trunkSizes.add (totalNumHaps-1);
				trunkSizes.add (1);
			}
			else {
				// what ever
				assert ((this.numCsdsPerPerm > 2) && (this.numCsdsPerPerm < totalNumHaps));
				
				// this is the difficult one
				int highestIdx = this.numCsdsPerPerm-1;
				for (int idx=0; idx <= highestIdx; idx++) {
					int sizeToAdd = (int)(1 + (idx * ((totalNumHaps-2)/((float)highestIdx))));
					trunkSizes.add(sizeToAdd);
				}
			}
			
			System.out.println ("# trunk sizes: " + trunkSizes);
			
			for (int i = permutation.size() - 1; i >= 0; i--) {
				int hapIdx = permutation.get (i);
				
				// get the deme and haplotype for this number
				DemeHapPair<HFSAXFullHaplotypeShell> thisPair = this.demeHaplotpyeList.get(hapIdx);
				
				int currSize = workingConfiguration.getTotalGeneticTypes();

				// put it there only if we want to
				if (trunkSizes.contains(currSize)) {
					DemoConfiguration<HFSAXFullHaplotypeShell> thisConfig = DemoConfiguration.copyConfiguration (workingConfiguration);
					thisPerm.add (new CSDConfig<HFSAXFullHaplotypeShell>(thisConfig, thisPair, this.pSet, this.fancyTransitionMap));
				}				
				// update the working configuration
				workingConfiguration.getPopulation(thisPair.demeNumber).adjustType(thisPair.haplotype, 1);				
			}
			
			// remember things from this permutatioan
			this.csdConfigs.add(thisPerm);
			
			// success
			return true;
		}

		public List<List<CSDConfig<HFSAXFullHaplotypeShell>>> getCsdConfigs() {
			return this.csdConfigs;
		}
	}

	public static class DiCalThread implements Callable<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>> {

		private final CSDConfig<HFSAXFullHaplotypeShell> config;
		private final double recoScaling;
		private final double[] theta;
		private final ConditionalObjectiveFunctionType objectiveType;
		private final CoreCache<HFSAXFullHaplotypeShell> demoCoreCache;
		private final boolean lastIteration;
		private final DemoStateFactory demoStateFactory;
		private final TrunkProcessFactory trunkFactory;
		private final Demography demo;
		private final boolean useEigenCore;
		
		
		public DiCalThread (CSDConfig<HFSAXFullHaplotypeShell> config, CoreCache<HFSAXFullHaplotypeShell> demoCoreCache, ConditionalObjectiveFunctionType objectiveType, double recoScaling, double[] theta, boolean lastIteration, boolean useEigenCore) {
			this(config, demoCoreCache, objectiveType, recoScaling, theta, lastIteration, null, null, null, useEigenCore);
		}
		public DiCalThread (CSDConfig<HFSAXFullHaplotypeShell> config, CoreCache<HFSAXFullHaplotypeShell> demoCoreCache, ConditionalObjectiveFunctionType objectiveType, double recoScaling, double[] theta, boolean lastIteration, DemoStateFactory demoStateFactory, TrunkProcessFactory trunkFactory, Demography demo, boolean useEigenCore) {
			this.config = config;
			this.recoScaling = recoScaling;
			this.theta= theta;
			this.objectiveType = objectiveType;
			this.demoCoreCache = demoCoreCache;
			this.lastIteration = lastIteration;
			
			//either you are caching cores, or you are not and then you need all of these factories and a demo
			assert(demoCoreCache == null || (demoStateFactory == null && trunkFactory == null && demo == null));
			this.demoStateFactory = demoStateFactory;
			this.trunkFactory = trunkFactory;
			this.demo = demo;
			
			this.useEigenCore = useEigenCore;
		}
		
		@Override
		public SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> call() throws Exception {

			UberDemographyCore core = null;
			
			assert ((demoCoreCache == null) || (demoCoreCache.isEigenCore() == this.useEigenCore));
			
			if(demoCoreCache == null){
				assert (this.demo != null);
				DemoStateCollection demoStates = demoStateFactory.getDemoStates (this.demo, this.config.trunkConfig, this.config.additionalHapDeme.demeNumber);
				TrunkProcess trunk = trunkFactory.getTrunk(this.config.trunkConfig.getSampleSizes(), demoStates);
				core = UberDemographyCore.getDemographyCore (trunk, demoStates,  this.config.additionalHapDeme.demeNumber, false, this.useEigenCore);
			} else {
				core  = this.demoCoreCache.getCore (this.config.trunkConfig, this.config.additionalHapDeme.demeNumber);
			}
			

			// now get some conditional sampling distribution
			SMCSDemo<HFSAXFullHaplotypeShell> csd = new SMCSDemo<HFSAXFullHaplotypeShell> (this.config.trunkConfig, core, this.config.pSet, this.recoScaling, this.config.fancyTransitionMap, this.theta);

			// remember the expected values
			// this should not really need the demoStateFactory, does it? maybe it does
			// anyways, if we are the last iteration, then don't need full table stuff
			SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> objectiveFunction = getObjectiveFunction (this.objectiveType, csd, this.config.additionalHapDeme.haplotype, lastIteration);
			
			if (this.config.getLogFactor() != null) {
				objectiveFunction.setLogFactor (this.config.getLogFactor());
			}
			
			// and return the objective function together with the index
			return objectiveFunction;
		}
	}

	public static class DiCalObjectiveFunction extends DiCalObjectiveFunctionBlueprint {

		public DiCalObjectiveFunction (List<List<List<CSDConfig<HFSAXFullHaplotypeShell>>>> csdConfigList, DemographyFactory demoFactory, DemoStateFactory demoStateFactory, ConditionalObjectiveFunctionType objectiveType, double[] startPoint, boolean estimateRecomScaling, TrunkProcessFactory trunkFactory, boolean estimateTheta, double defaultTheta, int thetaDim, boolean lastIteration, boolean verbose, boolean useEigenCore) {
			super(demoFactory, demoStateFactory, trunkFactory, estimateRecomScaling, estimateTheta, defaultTheta, thetaDim, verbose, useEigenCore);
			
			assert (trunkFactory.trunkStyle != TrunkProcess.TrunkProcessFactory.TrunkStyle.RecursiveTrunk);
			
			// also get some demography
			Demography demo = this.demoFactory.getDemography (this.getDemoPoint(startPoint));
			// should not happen
			assert (demo != null);

			// prepare the cache
			CoreCache<HFSAXFullHaplotypeShell> demoCoreCache;
			if(!estimateTheta){
				demoCoreCache = new CoreCache<HFSAXFullHaplotypeShell> (demo, trunkFactory, this.demoStateFactory, UberDemographyCore.DEFAULT_RENORMALIZATION_EPSILON, super.useEigenCore);
			} else {
				demoCoreCache = null;
			}
			
			// make space for the future
			this.futureExpectations = new ArrayList<List<List<Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>>>();

			// go through the configs and make the corres ponding threads
			for (List<List<CSDConfig<HFSAXFullHaplotypeShell>>> currChunkConfigs : csdConfigList) {
				List<List<Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>> chunkFutures = new ArrayList<List<Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>>();
				this.futureExpectations.add (chunkFutures);
				for (List<CSDConfig<HFSAXFullHaplotypeShell>> currPermConfigs : currChunkConfigs) {
					List<Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>> permsFutures = new ArrayList<Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>> ();
					chunkFutures.add(permsFutures);
					for (CSDConfig<HFSAXFullHaplotypeShell> config : currPermConfigs) {
						
						assert ((demoCoreCache == null) | (demoCoreCache.isEigenCore() == super.useEigenCore));

						// make this thread
						DiCalThread thread;
						if(!estimateTheta){
							thread = new DiCalThread (config, demoCoreCache, objectiveType, this.getRecomScaling(startPoint), this.getTheta(startPoint), lastIteration, demoCoreCache.isEigenCore());
						} else {
							thread = new DiCalThread (config, null, objectiveType, this.getRecomScaling(startPoint), this.getTheta(startPoint, demo, config.trunkConfig, config.additionalHapDeme.demeNumber), lastIteration, demoStateFactory, trunkFactory, demo, super.useEigenCore);
						}

						// submit it to the executor and remmeber the future
						Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>> thisFuture = null;
						// are we in a fork-join pool?
						if (ForkJoinTask.inForkJoinPool()) {
							// Parallel!
							// get a fork
							ForkJoinTask<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>> myFork = ForkJoinTask.adapt(thread);
							// fork it
							myFork.fork();
							// and put it in the drawer
							thisFuture = myFork;
						}
						else {
							// no parallel
							FutureTask<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>> futureTask = new FutureTask<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>(thread);
							futureTask.run();
							thisFuture = futureTask;
						}
						
						// remeber the future anyway
						assert (thisFuture != null);
						permsFutures.add(thisFuture);
					}
				}
			}
		}
		
		
		private Map<List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>, Double> chunkLogLikelihoodTimesK;
		
		private List<List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>> expectations;

		private final List<List<List<Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>>> futureExpectations;
		
		private Double logLikelihood;
		
		private boolean tablesFilled() {
			if ((this.chunkLogLikelihoodTimesK == null) && (this.expectations == null) && (this.logLikelihood == null)) {
				return false;
			}
			else if ((this.chunkLogLikelihoodTimesK != null) && (this.expectations != null) && (this.logLikelihood != null)) {
				return true;
			}
			else {
				// not possible
				assert false;
				return false;
			}
		}
		
		private void fillTables() {
			// make some map
			this.chunkLogLikelihoodTimesK = new HashMap<List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>, Double>();

			// collect them expectations
			// make the expectations
			double fullLogLike = 0d;
			this.expectations = new ArrayList<List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>>();
			// loop over chunks
			for (List<List<Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>> chunksFutures : this.futureExpectations) {
				List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>> currChunkObjectiveFunctions = new ArrayList<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>();
				this.expectations.add(currChunkObjectiveFunctions);
				
				// loop over perms
				for (List<Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>> permFutures : chunksFutures) {
					List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>> currPermObjectiveFunctions = new ArrayList<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>();
					currChunkObjectiveFunctions.add(currPermObjectiveFunctions);
					
					// loop over internal CSD configs
					for (Future<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>> daFuture : permFutures) {

						// get the current guy
						SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> objectiveFunction = null;
						try {
							// get it
							objectiveFunction = daFuture.get();
						} catch (InterruptedException e) {
							System.err.println ("Interrupted excpetion in dical thread:");
							e.getCause().printStackTrace(System.err); // this is more informative than printing the stack trace from e
							System.exit(-1);
						} catch (ExecutionException e) {
							System.err.println ("execution exception in dical thread:");
							e.getCause().printStackTrace(System.err); // this is more informative than printing the stack trace from e
							System.exit(-1);
						}
						assert (objectiveFunction != null);
						
						// and add it to the list in hopefully the right place
						currPermObjectiveFunctions.add (objectiveFunction);
					}
				}
		
				// remember the log likelihood for this chunk
				LogSum logLikSum = new LogSum(currChunkObjectiveFunctions.size());
				for (List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>> currExpectList : currChunkObjectiveFunctions) {
					double currLogLik = 0d;
					for (SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> currExpect : currExpectList) {
						if (currExpect.getLogFactor() == null) {
							currLogLik += currExpect.getLogLikelihood();
						}
						else {
							currLogLik += currExpect.getLogFactor() * currExpect.getLogLikelihood();
						}
					}
					
					logLikSum.addLogSummand(currLogLik);
				}
			
				// we might want it with the normalizing factor
				double d = logLikSum.retrieveLogSum();
				
				assert (!Double.isNaN(d));
				this.chunkLogLikelihoodTimesK.put(currChunkObjectiveFunctions, d);
				// and remember it for full
				// WARNING this factor of K should be here to get a true log-likelihood
				// however, we don't want to store it for the M-step (cause there it cancels)
				// d -= Math.log(currChunkObjectiveFunctions.size());
				fullLogLike += d;
			}
			// should be done
			// store the real logLikelihood
			this.logLikelihood = fullLogLike;
			
		}
		
		@Override
		public double getChunkLikelihood(List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>> currChunkList) {
			if (!this.tablesFilled()) {
				this.fillTables();
			}
			return this.chunkLogLikelihoodTimesK.get(currChunkList);
		}

		@Override
		public List<List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>> getExpectations() {
			if (!this.tablesFilled()) {
				this.fillTables();
			}
			assert (this.expectations != null);
			return this.expectations;
		}

		@Override
		public double getLogLikelihood() {
			if (!this.tablesFilled()) {
				this.fillTables();
			}
			assert (this.logLikelihood != null);
			return this.logLikelihood;
		}

	}
		
	public static abstract class DiCalObjectiveFunctionBlueprint extends CompositeLikelihoodObjectiveFunction {

		abstract public double getChunkLikelihood(List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>> currPermList);
		
		abstract public List<List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>>> getExpectations();
		
		public DiCalObjectiveFunctionBlueprint(DemographyFactory demoFactory, DemoStateFactory demoStateFactory, TrunkProcessFactory trunkFactory, boolean estimateRecomScaling, boolean estimateTheta, double defaultTheta, int thetaDim, boolean verbose, boolean useEigenCore) {
			// just redirect to super
			super(demoFactory, demoStateFactory, trunkFactory, estimateRecomScaling, estimateTheta, defaultTheta, thetaDim, verbose, useEigenCore);
		}

		public static class KlThread implements Callable<Pair<Double,Double>> {

			private final SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> csdOj;
			private final CoreCache<HFSAXFullHaplotypeShell> demoCoreCache;
			private final double recoScaling;
			private final double[] thetas;
			private final DemoStateFactory demoStateFactory;
			private final TrunkProcessFactory trunkFactory;
			private final Demography demo;
			private final boolean useEigenCore;

			public KlThread (SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> csdOj, CoreCache<HFSAXFullHaplotypeShell> demoCoreCache, double recoScaling, double[] thetas, boolean useEigenCore) {
				this(csdOj, demoCoreCache, recoScaling, thetas, null, null, null, useEigenCore);
			}
			public KlThread (SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> csdOj, CoreCache<HFSAXFullHaplotypeShell> demoCoreCache, double recoScaling, double[] thetas, DemoStateFactory demoStateFactory, TrunkProcessFactory trunkFactory, Demography demo, boolean useEigenCore) {
				this.csdOj = csdOj;
				this.demoCoreCache= demoCoreCache;
				this.recoScaling = recoScaling;
				this.thetas = thetas;
				assert(demoCoreCache == null || (demoStateFactory == null && trunkFactory == null && demo == null));
				this.demoStateFactory = demoStateFactory;
				this.trunkFactory = trunkFactory;
				this.demo = demo;
				
				this.useEigenCore = useEigenCore;
				assert ((demoCoreCache == null) || (demoCoreCache.isEigenCore() == useEigenCore));
			}
			
			
			@Override
			public Pair<Double,Double> call() throws Exception {

				double mValue;
				
				assert ((demoCoreCache == null) || (demoCoreCache.isEigenCore() == this.useEigenCore));
				
				if(demoCoreCache != null){
					mValue = this.csdOj.value (this.demoCoreCache, this.recoScaling, this.thetas);
				} else {
					DemoStateCollection demoStates = demoStateFactory.getDemoStates(demo, csdOj.getPresentConfig(), csdOj.oldCore.getObservedPresentDeme());
					TrunkProcess trunk = trunkFactory.getTrunk(csdOj.getPresentConfig().getSampleSizes(), demoStates);
					UberDemographyCore core = UberDemographyCore.getDemographyCore(trunk, demoStates, csdOj.oldCore.getObservedPresentDeme(), false, this.useEigenCore);
					mValue = this.csdOj.value(core, recoScaling, thetas);
				}
				double likelihood = this.csdOj.getLogLikelihood();
				
				if (this.csdOj.getLogFactor() == null) {
					return new Pair<Double,Double> (mValue, likelihood);
				}
				else {
					return new Pair<Double,Double> (this.csdOj.getLogFactor() * mValue, this.csdOj.getLogFactor() * likelihood);
				}
			}
		}
		
		
		@Override
		protected double getValue(double[] point) {

			// get new demography
			Demography newDemo = this.demoFactory.getDemography(this.getDemoPoint(point));
			
			// did we get a valid demography?
			if (newDemo == null) {
				return Double.POSITIVE_INFINITY;
			}

			// do a cache
			CoreCache<HFSAXFullHaplotypeShell> demoCoreCache;
			if(!estimateTheta){
				demoCoreCache= new CoreCache<HFSAXFullHaplotypeShell> (newDemo, this.trunkFactory, this.demoStateFactory, UberDemographyCore.DEFAULT_RENORMALIZATION_EPSILON, super.useEigenCore);
			} else {
				demoCoreCache = null;
			}

			// some map
			Map<List<List<Future<Pair<Double,Double>>>>, Double> futureChunkLogLikelihood = new HashMap<List<List<Future<Pair<Double,Double>>>>, Double>();

			
			// and loop over all the inner objective functions
			List<List<List<Future<Pair<Double,Double>>>>> futureExpectations = new ArrayList<List<List<Future<Pair<Double,Double>>>>>();
			for (List<List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>>> currChunkList : this.getExpectations()) {
				List<List<Future<Pair<Double,Double>>>> futureChunk = new ArrayList<List<Future<Pair<Double,Double>>>>();
				futureExpectations.add (futureChunk);
				for (List<SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell>> currPermList : currChunkList) {
					List<Future<Pair<Double,Double>>> futurePerm = new ArrayList<Future<Pair<Double,Double>>>();
					futureChunk.add (futurePerm);
					for (SMCSDemoEMObjectiveFunction<HFSAXFullHaplotypeShell> currCSD : currPermList) {
						// make a thread
						KlThread thread;
						
						assert ((demoCoreCache == null) || (demoCoreCache.isEigenCore() == super.useEigenCore));

						if(!estimateTheta){
							thread = new KlThread (currCSD, demoCoreCache, this.getRecomScaling(point), this.getTheta(point), super.useEigenCore);
						} else {
							thread = new KlThread (currCSD, demoCoreCache, this.getRecomScaling(point), this.getTheta(point, demoFactory.getDemography(this.getDemoPoint(point)), currCSD.getPresentConfig(), currCSD.oldCore.getObservedPresentDeme()), demoStateFactory, trunkFactory, demoFactory.getDemography(this.getDemoPoint(point)), super.useEigenCore);
						}
						
						// you have a bright future behind you
						Future<Pair<Double,Double>> theFuture = null;
						// and submit it or run it
						if (ForkJoinTask.inForkJoinPool()) {
							// Parallel!
							// get a fork
							ForkJoinTask<Pair<Double,Double>> myFork = ForkJoinTask.adapt(thread);
							// fork it
							myFork.fork();
							// and put it in the drawer
							theFuture = myFork;
						}
						else {
							FutureTask<Pair<Double,Double>> futureTask = new FutureTask<Pair<Double,Double>>(thread);
							futureTask.run();
							theFuture = futureTask;							
						}
						assert (theFuture != null);
						// into the future
						futurePerm.add (theFuture);
					}
				}			
				// IMPORTANT after it is filled, cause of definition of hashCode and equals 
				futureChunkLogLikelihood.put(futureChunk, this.getChunkLikelihood (currChunkList));
			}
			
			// now get the actual values
			double prodValue = 0d;
			// and loop over all the inner objective functions
			for (List<List<Future<Pair<Double,Double>>>> currChunkList : futureExpectations) {
				double currChunkValue = 0d;
				for (List<Future<Pair<Double,Double>>> currPermList : currChunkList) {
					double currPermValue = 0d;
					double currPermLogLikelihood = 0d;
					for (Future<Pair<Double,Double>> currCSD : currPermList) { 
						
						// get it from the future
						Pair<Double, Double> valuePair = null;
						try {
							valuePair = currCSD.get();
						} catch (InterruptedException e) {
							System.err.println("Interrupted excpetion in kl thread:");
							e.getCause().printStackTrace(System.err); // this is more informative than printing the stack trace from e
							System.exit(-1);
						} catch (ExecutionException e) {
							System.err.println("execution exception in kl thread:");
							e.getCause().printStackTrace(System.err); // this is more informative than printing the stack trace from e
							System.exit(-1);
						}
								
								
						// get intervals
						currPermValue += valuePair.first();
						assert (!Double.isNaN(currPermValue));
						
						
						// this might happen if some step was allowed under the old parameters,
						// but it is prohibited under the new ones
						if (Double.isInfinite(currPermValue)) {
							// we can never recover from this, so might as well stop
							// also, we have prevent this Infinity from getting multiplied with 0
							return Double.POSITIVE_INFINITY;
						}
						
						
						currPermLogLikelihood += valuePair.second();
						assert (!Double.isNaN(currPermLogLikelihood));
					}
					// we have to normalize for numerical stability [by the log-likelihood]
					// [WARNING] for large differences this produces NaNs
					double weightOfCurrPerm = Math.exp(currPermLogLikelihood - futureChunkLogLikelihood.get(currChunkList));
					assert (!Double.isNaN(weightOfCurrPerm));
					currPermValue *= weightOfCurrPerm;
					currChunkValue += currPermValue;
					assert (!Double.isNaN(currChunkValue));
				}
				
				
				//return it (remember, we wan to use a minimizer)
				assert (!Double.isNaN(currChunkValue));
				// divide by size in non-log space since it is an average of sorts
				prodValue += (Math.log(currChunkList.size()) - currChunkValue);
			}
			
			// give it away
			return prodValue;
		}

	}

	// composite likelihood
	public static abstract class CompositeLikelihoodObjectiveFunction implements MultivariateFunction{

		CompositeLikelihoodObjectiveFunction (DemographyFactory demoFactory, DemoStateFactory demoStateFactory, TrunkProcessFactory trunkFactory, boolean estimateRecomScaling, boolean estimateTheta, double defaultTheta, int thetaDim, boolean verbose, boolean useEigenCore) {
			// store the factory
			this.demoFactory = demoFactory;
			this.demoStateFactory = demoStateFactory;
			this.trunkFactory = trunkFactory;
			this.estimateRecomScaling = estimateRecomScaling;
			assert (estimateTheta == (thetaDim > 0));
			this.estimateTheta = estimateTheta;
			this.defaultTheta = defaultTheta;
			this.thetaDim = thetaDim;
			this.verbose = verbose;
			
			this.useEigenCore = useEigenCore;
			
			this.memoizedSimplex = new HashMap<SimplexPointKey, Double>();
		}

		// function to get the likelihood
		public abstract double getLogLikelihood ();
		
		// wrapper for the value function
		@Override
		public double value (double[] point) {
			// get the time
			long startTime = System.currentTimeMillis();
						
			double result = Double.NaN;
			
			// have we already computed it?
//			MetaOptimization.synchronizedPrintln (String.valueOf (this.memoizedSimplex.size()));
			SimplexPointKey key = new SimplexPointKey (point);
			Double daResult = this.memoizedSimplex.get (key);
			
			if (daResult == null) {
				// no, so compute it

				// call the real function
				result = this.getValue (point);
				
				// and also store it
				this.memoizedSimplex.put(key, result);
			}
			else {
				// we have already computed it
				result = daResult.doubleValue();
			}
			
			if (this.verbose) {
				// print it (but with a minus in front of the result, since we Maximize the logLike)
				// also, print with time
				long elapsedTime = System.currentTimeMillis() - startTime;
				String debugString = "# " + Arrays.toString(point) + ": " + (-1)*result + "\n";
				debugString += "# [Q_TIME] " + elapsedTime;
				// memoized or not
				if (daResult == null) {
					debugString += " (new)";
				}
				else {
					debugString += " (memoized)";
				}
				MetaOptimization.synchronizedPrintln(debugString);
			}

			assert (!Double.isNaN(result));
			return result;
		}
		
		protected double[] getDemoPoint(double[] point) {
			return getDemoPointStatic (point, this.estimateRecomScaling, this.estimateTheta, this.thetaDim);
		}
		
		public static double[] getDemoPointStatic (double[] point, boolean estimateRecomScaling, boolean estimateTheta, int thetaDim) {
			double[] demoPoint = new double[point.length - (estimateRecomScaling ? 1 : 0) - (estimateTheta ? thetaDim : 0)];
			for (int i = 0; i < demoPoint.length; i++) {
				demoPoint[i] = point[i];
			}
			return demoPoint;
		}
		
		protected double getRecomScaling(double[] point) {
			return getRecomScalingStatic (point, this.estimateRecomScaling, this.estimateTheta, this.thetaDim);
		}
		
		public static double getRecomScalingStatic (double[] point, boolean estimateRecomScaling, boolean estimateTheta, int thetaDim) {
			if (estimateRecomScaling) {
				return point[point.length - 1 - (estimateTheta ? thetaDim : 0)];
			} else {
				return 1d;
			}
		}
		
		protected double[] getTheta(double[] point, Demography demo, DemoConfiguration<HFSAXFullHaplotypeShell> presentConfig, int additionalDeme){
			return getThetaStatic(point, this.estimateTheta, this.thetaDim, this.defaultTheta, demo, presentConfig, additionalDeme);
		}
		
		protected double[] getTheta(double[] point){
			return getTheta(point, null, null, -1);
		}
		
		public static double[] getThetaStatic (double[] point, boolean estimateTheta, int thetaDim, double defaultTheta, Demography demo
				,DemoConfiguration<HFSAXFullHaplotypeShell> presentConfig, int additionalDeme) {
			double[] theta;
			if(estimateTheta){
				assert demo != null;
				assert presentConfig != null;
				double[] unrefinedTheta= new double[thetaDim];
				for(int i=0; i< unrefinedTheta.length; i++){
					unrefinedTheta[i]= point[point.length - thetaDim + i];
				}
				
				List<IntervalFactory> facts = new ArrayList<IntervalFactory>(2);
				facts.add(StructureEstimationEM.otherIntervalFactory);
				facts.add(StructureEstimationEM.thetaIntervalFactory);
				IntervalFactory refiner = new IntervalFactory.CombinedIntervalFactory(facts, 1e-13, false);
				Interval[] refinedIntervals = refiner.getNontrivialIntervals(demo, presentConfig, additionalDeme);
				Interval[] thetaIntervals = StructureEstimationEM.thetaIntervalFactory.getNontrivialIntervals(demo, presentConfig, additionalDeme);
				if(thetaIntervals[thetaIntervals.length - 1].endPoint != Double.POSITIVE_INFINITY){
					throw new RuntimeException("bad theta intervals");
				}
				int currUnrefinedIntervalIdx = 0;
				theta = new double[refinedIntervals.length];
				for(int i=0; i < refinedIntervals.length; i++){
					while(refinedIntervals[i].endPoint > thetaIntervals[currUnrefinedIntervalIdx].endPoint){
						currUnrefinedIntervalIdx++;
					}
					theta[i] = unrefinedTheta[currUnrefinedIntervalIdx];
				}
				//
			} else {
				assert demo == null;
				assert presentConfig == null;
				assert additionalDeme == -1;
				// no estimation wanted
				theta = new double[] {defaultTheta};
			}
			return theta;
		}
		
		// the real value function
		protected abstract double getValue (double[] point);
		
		// get the dimension
		public int getDimension() {
			return this.demoFactory.getDimension() + (this.estimateRecomScaling ? 1 : 0) + (this.estimateTheta ? this.thetaDim : 0);
		}

		// here is the factory
		protected final DemographyFactory demoFactory;
		protected final DemoStateFactory demoStateFactory;
		protected final TrunkProcessFactory trunkFactory;
		protected boolean estimateRecomScaling;
		protected boolean estimateTheta;
		private double defaultTheta;
		private int thetaDim;
		private boolean verbose;
		protected final boolean useEigenCore;
		
		protected final Map<SimplexPointKey, Double> memoizedSimplex;
		
		// lookup class for points on simplex 
		private static class SimplexPointKey {
			// variables for look up
			private double[] key;
			
			public double[] getKey () {
				return key;
			}
					
			public SimplexPointKey (double[] point) {
				// make it
				this.key = new double[point.length];
				// fill it
				for (int i=0; i<point.length; i++) {
					this.key[i] = point[i];
				}
			}

		    @Override
		    public int hashCode() {
		        final int prime = 31;
		        int result = 1;
				for (int i=0; i<this.key.length; i++) {
					int mantissa = (int)(Double.doubleToLongBits (this.key[i]) & 0x000fffffffffffffL);
					result = prime * result + mantissa;
				}
//				MetaOptimization.synchronizedPrintln ("mantissa " + result);
		        return result;
		    }

			@Override
			public boolean equals (Object o) {
//				MetaOptimization.synchronizedPrintln ("equals");
				assert (o instanceof SimplexPointKey);
				double[] newKey =  ((SimplexPointKey) o).getKey();
				assert (newKey.length == this.key.length);
				// have to literally be the same
				boolean same = true;
				for (int i=0; i<newKey.length; i++) {
					same &= (this.key[i] == newKey[i]);
//					MetaOptimization.synchronizedPrintln (String.valueOf (this.key[i] - newKey[i]));
				}
				// give it away now
				return same;
			}
		}

	}

	
	// comparing things
	public static class PointValueCompare implements Comparator<PointValuePair> {
		
		@Override
		public int compare(PointValuePair o1, PointValuePair o2) {
			if (o1.getValue() < o2.getValue()) {
				return -1;
			}
			else if (o1.getValue() == o2.getValue()) {
				return 0;
			}
			else {
				return 1;
			}
		}
	}	
}
