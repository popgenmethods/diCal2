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

package edu.berkeley.diCal2.csd;

import java.util.ArrayList;
import java.util.Arrays;

import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.haplotype.FSAHaplotype;
import edu.berkeley.diCal2.haplotype.FSAParamSet;
import edu.berkeley.diCal2.haplotype.GeneticType;
import edu.berkeley.diCal2.haplotype.GeneticType.GeneticTypeMultiplicity;
import edu.berkeley.diCal2.utility.LogSum;
import edu.berkeley.diCal2.utility.MemSaveArrays.TwoDimDoubleArray;
import edu.berkeley.diCal2.utility.RealPartition.Interval;

public class SMCSDemo<H extends FSAHaplotype> {

	public static <H extends FSAHaplotype> Interval[] getDefaultIntervals (int numAdditionalIntervals, Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme, boolean weddingCake, boolean printIntervals) {
		return (new IntervalFactory.OldIntervalFactory(numAdditionalIntervals, printIntervals)).getIntervals(demo, presentConfig, additionalDeme);
	}

	// baby baby constructor
	public SMCSDemo (DemoConfiguration<H> presentConfig, int observedDeme, EigenParamSet pSet, Demography demo, int numAdditionalIntervals, boolean useEigenCore) {
		// get some core
		this (presentConfig, observedDeme, pSet, 1d, demo, numAdditionalIntervals, useEigenCore);
	}
	
	// baby constructor
	public SMCSDemo (DemoConfiguration<H> presentConfig, int observedDeme, EigenParamSet pSet, double recoScaling, Demography demo, int numAdditionalIntervals, boolean useEigenCore) {
		// get some core
		this (presentConfig, observedDeme, pSet, recoScaling, demo, getDefaultIntervals (numAdditionalIntervals, demo, presentConfig, observedDeme, false, false), useEigenCore);
	}
	
	//bigger baby constructor
	public SMCSDemo (DemoConfiguration<H> presentConfig, int observedDeme, EigenParamSet pSet, Demography demo, int numAdditionalIntervals, boolean useEigenCore, HmmStepHandler<H> decodeLocusTransitionMap) {
		// get some core
		this (presentConfig, observedDeme, pSet, 1d, demo, getDefaultIntervals (numAdditionalIntervals, demo, presentConfig, observedDeme, false, false), useEigenCore, decodeLocusTransitionMap);
	}

	
	
	public SMCSDemo (DemoConfiguration<H> presentConfig, int observedDeme, EigenParamSet pSet, double recoScaling, Demography demo, Interval[] stateIntervals, boolean useEigenCore) {
		this (presentConfig, UberDemographyCore.getDemographyCore(presentConfig.getSampleSizes(), observedDeme, demo, stateIntervals, false, useEigenCore), pSet, recoScaling, new HmmStepHandler.SingleStepLocusTransitionMap<H>(pSet));
	}
	
	public SMCSDemo (DemoConfiguration<H> presentConfig, int observedDeme, EigenParamSet pSet, double recoScaling, Demography demo, Interval[] stateIntervals, boolean useEigenCore, HmmStepHandler<H> decodeLocusTransitionMap) {
		this(presentConfig, UberDemographyCore.getDemographyCore(presentConfig.getSampleSizes(), observedDeme, demo, stateIntervals, false, useEigenCore), pSet, recoScaling, decodeLocusTransitionMap);
	}
	
	public SMCSDemo (DemoConfiguration<H> presentConfig, UberDemographyCore core, EigenParamSet pSet, double recoScaling, HmmStepHandler<H> decodeLocusTransitionMap){
		this(presentConfig, core, pSet, recoScaling, decodeLocusTransitionMap, null);  //
	}
	
	
	public SMCSDemo (DemoConfiguration<H> presentConfig, UberDemographyCore core, EigenParamSet pSet, double recoScaling, HmmStepHandler<H> decodeLocusTransitionMap, double[] theta) {
		
		this.presentConfig = presentConfig;
		this.trunkAbsorptionLogFractions = computeLogTrunkAbsorptionFractions (this.presentConfig);
		assert (core.getTrunk() == null || Arrays.equals(presentConfig.getSampleSizes(), core.getTrunk().getSampleSizes()));
		
		assert (decodeLocusTransitionMap.isValidConfig (this.getPresentConfig().getTotalPopulation()));
		
		// save the decode loci
		this.decodeLocusTransitionMap = decodeLocusTransitionMap;
		
		// save the pSet
		this.pSet = pSet;
		// be sure
		// with the new type for all parameters same everywhere, this is given
//		assert (checkMutationMatricesEqual (pSet, UberDemographyCore.ONELOCUS_EPSILON, UberDemographyCore.RELATIVE_ERROR));
		
		// save the core
		this.core = core;
		
		// and use it to get proberbilities
		// again, if the conditional config is empty, we don't have to do anything
		if (!this.core.isConditionalConfigEmpty()) {
			
			// initialize the initial probabilities
			this.initialLogProbabilities = this.core.getDemoStateLogMarginal();
			
			// initialize the transition probabilities
			decodeLocusTransitionMap.updateTheta(theta);
			this.hmmProbsList = decodeLocusTransitionMap.computeTransitions (this.core, recoScaling);
			
			// and the emission probabilities
			this.logEmissionArray = computeEmissionLogArray (this.decodeLocusTransitionMap, pSet, this.core, theta);
			
			// set up the hidden states
			this.firstHapStateIdxMap = createFirstHapStateIdxMap (this.core, this.presentConfig);
		}
	}
	
	public static <H extends GeneticType> int[] createFirstHapStateIdxMap (UberDemographyCore core, DemoConfiguration<H> presentConfig) {
		// make a new one
		int[] toReturn = new int[core.numDemoStates()+1];

		// and fill it
		int firstHapIdx = 0;
		for (int demoState=0; demoState < core.numDemoStates(); demoState++) {
			// store it
			toReturn[demoState] = firstHapIdx;
			
			// and go further
			int presentDemeIdx = core.getDemoStateCollection().getPresentDeme(demoState);
			firstHapIdx += presentConfig.getPopulation(presentDemeIdx).distinctGeneticTypes();
		}
		
		// and put the total number in at the end
		toReturn[core.numDemoStates()] = firstHapIdx;
		
//		MetaOptimization.synchronizedPrintln (Arrays.toString (toReturn));
		
		// give it away now
		return toReturn;
	}

	public static boolean checkMutationMatricesEqual (FSAParamSet pSet, double EPSILON, double RELATIVE_ERROR) {
		// everything's fine before we see anything
		boolean equal = true;
		
		// get first Matrix
		double[][] firstMatrix = pSet.getMutationMatrix(0).getArray();
		// check for square of first matrix
		int tmp = firstMatrix.length;
		for (double[] row : firstMatrix) {
			equal &= (tmp == row.length);
		}		

		// check whether the mutation matrices are the same on every loci
		// if that is not true, we don't do anything
		// now go through everything and check
		for (int l=1; l<pSet.numLoci(); l++) {
			double[][] thisMatrix = pSet.getMutationMatrix(l).getArray();
			// check yourself before you wreck yourself
			equal &= matrixEqual (firstMatrix, thisMatrix, EPSILON, RELATIVE_ERROR);
		}

		// whatcha say?
		return equal;
	}
	
	public static boolean matrixEqual (double[][] firstMatrix, double[][] thisMatrix, double EPSILON, double RELATIVE_ERROR) {
		//  I guess we can take squaredness of matrices for granted 
		boolean equal = (firstMatrix.length == thisMatrix.length);
		// now check the entries
		for (int i=0; i<firstMatrix.length; i++) {
			for (int j=0; j<firstMatrix.length; j++) {
				equal &= (((Math.abs (thisMatrix[i][j]) < EPSILON) && (Math.abs (firstMatrix[i][j]) < EPSILON))
						|| (Math.abs ((thisMatrix[i][j] - firstMatrix[i][j])/firstMatrix[i][j]) < RELATIVE_ERROR));
			}
		}
		// show it
		return equal;
	}

	public double getLogTrunkAbsorptionFraction(int presentDeme, int haplotypeIdx){
		return this.trunkAbsorptionLogFractions.get(presentDeme)[haplotypeIdx];
	}
	
	public static <H extends GeneticType> ArrayList<double[]> computeLogTrunkAbsorptionFractions (DemoConfiguration<H> config) {
		ArrayList<double[]> toReturn = new ArrayList<double[]>();
		
		for (int demeIdx = 0; demeIdx < config.getNumberDemes(); demeIdx++){
			HaplotypeConfiguration<H> currDeme = config.getPopulation(demeIdx);
			int currDemeSize = currDeme.totalGeneticTypes();
			
			double[] currDemeFactors = new double[currDemeSize];
			toReturn.add(currDemeFactors);
			
			double currDemeLogSize =  Math.log(currDemeSize);
			
			int currHapIdx = 0;
			for ( GeneticTypeMultiplicity<H> currHap : currDeme) {
				
				currDemeFactors[currHapIdx] = Math.log(currHap.multiplicity) - currDemeLogSize;
				
				currHapIdx++;
			}
		}
		
		return toReturn;
	}
	
	public static <H extends FSAHaplotype> double[][][][] computeEmissionLogArray(HmmStepHandler<H> decodeLocusTransitionMap, EigenParamSet pSet, UberDemographyCore core, double[] theta) {
	
		// how many decode loci?
		int numDecodeLoci = decodeLocusTransitionMap.numHmmSteps();

		// keep a list of mutation rates associated with the first position they occur
		ArrayList<RateIndexPair> firstRateOccurence = new ArrayList<RateIndexPair>();
		
		// initialize the emission array to nothing
		double[][][][] logEmissionArray = new double[numDecodeLoci][][][];
		
		// now go through and set everything
		// use old ones if you know them, calculate new one if you don't
		for (int decodeLocusIdx=0; decodeLocusIdx<numDecodeLoci; decodeLocusIdx++) {
			
			// Use the provided theta if we are estimating it, otherwise get one
			double[] mutationRate = theta != null ? theta : decodeLocusTransitionMap.getMutationRate(decodeLocusIdx);
			
			// try to find it
			RateIndexPair pair = findRate (mutationRate, firstRateOccurence, UberDemographyCore.RELATIVE_ERROR);
			
			// locals
			double[][][] emission;
			
			// did you find anything
			if (pair == null) {
				// apparently not
				
				// so add it to the list (remember the decode locus index where it is)
				firstRateOccurence.add (new RateIndexPair (mutationRate, decodeLocusIdx));
				
				// and actually compute the new ones
				// compute some emission probabilities for the current locus
				// use the eigenstuff from first matrix, cause should be the same everywhere
				emission = core.getLogEmission(mutationRate, pSet.getMutMatrixEigenDecomp());
				
			}
			else {
				// yes, we did find something, so use it
				emission = logEmissionArray[pair.index];
			}

			// now that we got them somehow, actually associate it with this locus
			logEmissionArray[decodeLocusIdx] = emission;
		}
		// should be all fine now
		return logEmissionArray;
	}

	public double computeConditionalProb(H addHap) {
		// compute via the logarithm
		return Math.exp(this.computeConditionalLogProb(addHap));
	}

	public double computeConditionalLogProb (H additionalHaplotype) {
		
		// if the conditional configuration is empty, then we just return the stationary sampling probability of one guy
		if (this.core.isConditionalConfigEmpty()) {

			// return the stationary distribution
			return this.decodeLocusTransitionMap.computeLogStationaryProbability(additionalHaplotype);
		}

		// otherwise
		// just pass the on the haplotype to more fancy stuff
		return this.forwardAlgorithm (additionalHaplotype, false, null, false);
	}
	
	public void fillForwardBackward (H additionalHaplotype) {
		
		// remember haplotype
		this.mostRecentAdditionalHaplotype = additionalHaplotype;
		
		// if the conditional configuration is empty, the forward backward thing should never be called
		assert (!this.core.isConditionalConfigEmpty());

		// allocate space for forward variables
		this.forwardLogProbabilities = new TwoDimDoubleArray (this.decodeLocusTransitionMap.numHmmSteps(), this.getNumHiddenStates());
		// and compute it (also store the likelihood)
		this.logLikelihood =  this.forwardAlgorithm (additionalHaplotype, true, this.forwardLogProbabilities, false);
		
		
		// allocate space for pre backward variables
		TwoDimDoubleArray preBackwardVariables = new TwoDimDoubleArray (this.decodeLocusTransitionMap.numHmmSteps(), this.getNumHiddenStates());
		// compute them
		this.forwardAlgorithm (additionalHaplotype, true, preBackwardVariables, true);
		
		// and modify them
		for (int decodeIdx=0; decodeIdx<this.decodeLocusTransitionMap.numHmmSteps(); decodeIdx++) {
			for (int demoState=0; demoState<this.core.numDemoStates(); demoState++) {
				int currDemeIdx = this.core.getDemoStateCollection().getPresentDeme(demoState);
				HaplotypeConfiguration<H> currPop = this.getPresentConfig().getPopulation(currDemeIdx);
				
				int hapIndex = 0;
				for (GeneticTypeMultiplicity<H> trunkHap : currPop) {
					// modify the pre variables accordingly (normalizing)
					double logDenominator = this.getLogTrunkAbsorptionFraction(currDemeIdx,hapIndex);
					logDenominator +=  this.decodeLocusTransitionMap.getLogEmissionProb (this.logEmissionArray[decodeIdx][demoState], decodeIdx, trunkHap.geneticType, currDemeIdx, additionalHaplotype);
					logDenominator +=  this.initialLogProbabilities[demoState];
					
					// if denominator is 0, assert preBackward is 0, and continue
					if (logDenominator == Double.NEGATIVE_INFINITY) {
						assert (preBackwardVariables.get (decodeIdx, this.getFirstHapStateIndex (demoState) + hapIndex) == Double.NEGATIVE_INFINITY);
						continue;
					}
					
					preBackwardVariables.adjust (decodeIdx, this.getFirstHapStateIndex (demoState) + hapIndex, -logDenominator);
					
					// increase index
					hapIndex++;
				}
			}
		}

		// and store 'em
		this.backwardLogProbabilities = preBackwardVariables;
		
		
		// remember that we did this
		this.FwBwComputed = true;
	}

	private int getNumHiddenStates() {
		assert (this.firstHapStateIdxMap.length == this.core.numDemoStates()+1);
		// ther last entry (one after regular) has the total number
		return this.firstHapStateIdxMap[this.core.numDemoStates()];
	}

	public static int getNextDecodeIdx (int currDecodeIdx, boolean reverseOrder) {
		return 	currDecodeIdx + (reverseOrder ? -1 : 1);
	}
	
	public static boolean continueForwardAlgorithm (int currDecodeIdx, boolean reverseOrder, int numDecodeLoci) {
		if (!reverseOrder) {
			return (currDecodeIdx < numDecodeLoci);
		}
		else {
			return (currDecodeIdx >= 0);
		}
	}
	
	private double forwardAlgorithm (H additionalHaplotype, boolean storeForwardVariables, TwoDimDoubleArray forwardVariables, boolean reverseOrder) {
		
		assert (storeForwardVariables || (forwardVariables == null));
		
		// just to be sure
		assert (!storeForwardVariables || (forwardVariables.getFirstSize() == this.decodeLocusTransitionMap.numHmmSteps()));
		assert (!storeForwardVariables || (forwardVariables.getSecondSize() == this.getNumHiddenStates()));
		
		// just some checking
		assert (this.decodeLocusTransitionMap.isValidHaplotype (additionalHaplotype));
		
		// get number of loci
		int numDecodeLoci = this.decodeLocusTransitionMap.numHmmSteps();
		
		// if the conditional configuration is empty, the forward backward thing should never be called
		assert (!this.core.isConditionalConfigEmpty());
		
		// now that we have (more or less) the probabilities for our hidden Markov model
		// let's calculate the forward of observing the given haplotype
	
		// start at first decode locus
		int currDecodeLocusIdx = reverseOrder ? numDecodeLoci - 1 : 0;
		
		// get the first guy
		double[] currLogF = new double[this.getNumHiddenStates()];

		// for now and testing and fun just the first locus, calculate probability while filling
		// which demoState can we be absorbed in?
		for (int absorbDemoState=0; absorbDemoState < this.core.numDemoStates(); absorbDemoState++) {
			
			int presentDemeIdx = this.core.getDemoStateCollection().getPresentDeme(absorbDemoState);
			HaplotypeConfiguration<H> currPop = this.getPresentConfig().getPopulation(presentDemeIdx);
			
			// get some space for this deme
			// which haplotype can we be absorbed into?
			int hapIndex = 0;
			for (GeneticTypeMultiplicity<H> trunkHap : currPop) {
				// do it in log scale, so addition
				double value = initialLogProbabilities[absorbDemoState];
				value += this.getLogTrunkAbsorptionFraction (presentDemeIdx, hapIndex);
				value += this.decodeLocusTransitionMap.getLogEmissionProb(this.logEmissionArray[currDecodeLocusIdx][absorbDemoState], currDecodeLocusIdx, trunkHap.geneticType, presentDemeIdx, additionalHaplotype);

				currLogF[this.getFirstHapStateIndex (absorbDemoState) + hapIndex] = value;
				hapIndex++;
			}
		}

		// save it, if you want it
		if (storeForwardVariables) {
			forwardVariables.setRow (currDecodeLocusIdx, currLogF);
		}
		
		// go to next locus (or is it the previous?)
		int prevDecodeLocusIdx = currDecodeLocusIdx;
		currDecodeLocusIdx = getNextDecodeIdx (currDecodeLocusIdx, reverseOrder) ;

		
		// the big loop
		while (continueForwardAlgorithm (currDecodeLocusIdx, reverseOrder, numDecodeLoci)) {
			
			// get the right Fs
			double[] prevLogF = currLogF;
			currLogF = new double[this.getNumHiddenStates()];
			
			
			// and then start with the actual calculations
			// first calculate the probabilities of being in a certain time in a certain deme for the current locus
			// cause when we have recombination, the transition is independent of the haplotype
			double[] logQ = new double[this.core.numDemoStates()];

			for (int absorbDemoState=0; absorbDemoState < this.core.numDemoStates(); absorbDemoState++) {
					// initialize
					// do computation in log scale
					logQ[absorbDemoState] = Double.NEGATIVE_INFINITY;
					// and fill
					for (int hapIdx=this.getFirstHapStateIndex(absorbDemoState); hapIdx<this.getFirstHapStateIndex(absorbDemoState+1); hapIdx++) {
						logQ[absorbDemoState] = LogSum.computePairLogSum (logQ[absorbDemoState], prevLogF[hapIdx]);
					}
			}
			// Q should be filled
			
			// now proceed to next locus
			
			// again we calculate the haplotype independent part here of the recombination-transition here
			double[] logR = new double[this.core.numDemoStates()];
			
			// get the right reco transitions (remember, it's at the previous locus to current)
			double[][] recoTransition = this.hmmProbsList.getLogReco (this.decodeLocusTransitionMap.getTransitionIndex(prevDecodeLocusIdx, reverseOrder));
			
			
			// where do we go
			for (int dstDemoState=0; dstDemoState < this.core.numDemoStates(); dstDemoState++) {
			
				// do all calculations in log scale
				// first there is nothing here
				logR[dstDemoState] = Double.NEGATIVE_INFINITY;
				// where do we come from possibly
				for (int srcDemoState=0; srcDemoState < this.core.numDemoStates(); srcDemoState++) {
					// log scale
					// [currDemoState][prevDemoState]
					logR[dstDemoState] = LogSum.computePairLogSum (logR[dstDemoState], recoTransition[srcDemoState][dstDemoState] + logQ[srcDemoState]);
				}
			}

			
			// get some fancy no reco pants (remember, it's at the previous locus to current)
			double[] noRecoTransition = this.hmmProbsList.getLogNoReco (this.decodeLocusTransitionMap.getTransitionIndex(prevDecodeLocusIdx, reverseOrder));

			// now the real hidden states times emission probability
			// what is the deme at next locus?
			
			for (int dstDemoState=0; dstDemoState < this.core.numDemoStates(); dstDemoState++) {
				
				int dstDemeIdx = this.core.getDemoStateCollection().getPresentDeme(dstDemoState);
				HaplotypeConfiguration<H> dstPop = this.getPresentConfig().getPopulation(dstDemeIdx);
				
				// and what is the haplotype at the next locus
				int hapIndex = 0;
				for (GeneticTypeMultiplicity<H> trunkHap : dstPop) {
					// do all calculations in log scale
					// what happened?
					// either there is no recombination and we stay
					// currLocus - 1, cause that's where we store it
					double tmp = noRecoTransition[dstDemoState] + prevLogF[this.getFirstHapStateIndex(dstDemoState) + hapIndex];
					// or there is recombination, we can come from everywhere and have to be absorbed in this haplotype
					// NOTE the following line adds to the value
					tmp = LogSum.computePairLogSum(tmp, logR[dstDemoState] + this.getLogTrunkAbsorptionFraction(this.core.getDemoStateCollection().getPresentDeme(dstDemoState), hapIndex));
					// either way, we have to emit the allele at the given position [multiply with that]
					tmp += this.decodeLocusTransitionMap.getLogEmissionProb(this.logEmissionArray[currDecodeLocusIdx][dstDemoState], currDecodeLocusIdx, trunkHap.geneticType, dstDemeIdx, additionalHaplotype);

					currLogF[this.getFirstHapStateIndex(dstDemoState) + hapIndex] = tmp;
					
					// increase index
					hapIndex++;
				}
			}
			
			// maybe remember it
			if (storeForwardVariables) {
				forwardVariables.setRow (currDecodeLocusIdx, currLogF);
			}
			
			// step to the next locus (in whatever direction you want)
			prevDecodeLocusIdx = currDecodeLocusIdx;
			currDecodeLocusIdx = getNextDecodeIdx (currDecodeLocusIdx, reverseOrder) ;
		}
		
		
		// compute the likelihood, necessary stuff should be in currLogF, just sum it
		LogSum thisLog = new LogSum (this.getNumHiddenStates());
		for (int i=0; i<currLogF.length; i++) {
			thisLog.addLogSummand (currLogF[i]);
		}
		double logLikelihood = thisLog.retrieveLogSum();

		// I actually don't think that this is needed anymore
//		// we have to reverse the order of the storage in some cases
//		if (storeForwardVariables && reverseOrder) {
//			Collections.reverse (forwardVariables);
//		}
		
		// give it away now
		return logLikelihood;
	}

	
	public int getFirstHapStateIndex (int demoState) {
		return this.firstHapStateIdxMap[demoState];
	}

	// method to find some index pair
	private static RateIndexPair findRate (double[] mutationRate, ArrayList<RateIndexPair> firstRateOccurence, double relativeError) {
		// go through linearly [for now] and find it
		// i guess in the future we can implement some binary search
		for (RateIndexPair pair : firstRateOccurence) {
			// is this the one?
			boolean match= true;
			for(int i= 0; i < mutationRate.length; i++){
				if (Math.abs ((mutationRate[i] - pair.rate[i])/pair.rate[i]) < relativeError) {
					
				}	
			}
			// yes it is, so return it
			if(match){
				return pair;
			}
		}
		// not found, so return null
		return null;
	}

	public boolean expectationsComputed() {
		return expectationsComputed;
	}

	/// sub class
	private static class RateIndexPair {
		// constructor
		public RateIndexPair (double[] rate, int index) {
			this.rate = rate;
			this.index = index;
		}
		// some member function stuff
		public double[] rate;
		public int index;
	}

	public double getLogLikelihood () {
		assert (this.FwBwComputed);

		return this.logLikelihood;
	}
	
	public double getPosteriorProb(int currDecodeIdx, int absDemoState, int absHapIdx) {
		assert (this.FwBwComputed);
		
		double jointLogProb = this.forwardLogProbabilities.get (currDecodeIdx, this.getFirstHapStateIndex(absDemoState) + absHapIdx);
		jointLogProb += this.backwardLogProbabilities.get (currDecodeIdx, this.getFirstHapStateIndex(absDemoState) + absHapIdx);
		double condProb = Math.exp (jointLogProb - this.logLikelihood);

		return condProb;
	}
	
	public ArrayList<double[]> getPosteriorProbs(int decodeLocusIdx) {
		ArrayList<double[]> currLocusPosterior = new ArrayList<double[]>();
		
		// this should also sum to one
		double posteriorSum = 0d;
		
		// iterate absorption interval deme
		for (int absDemoState=0; absDemoState < this.core.numDemoStates(); absDemoState++) {
			// iterate over absorption haplotype
			int currDemeIdx = this.core.getDemoStateCollection().getPresentDeme(absDemoState);
			HaplotypeConfiguration<H> currConfig = this.getPresentConfig().getPopulation(currDemeIdx);
			
			double[] currDemoStatePosterior = new double[currConfig.distinctGeneticTypes()];
			currLocusPosterior.add(currDemoStatePosterior);
			
			for (int absHapIdx=0; absHapIdx<currConfig.distinctGeneticTypes(); absHapIdx++) {
				
				double currProb = this.getPosteriorProb(decodeLocusIdx, absDemoState, absHapIdx);
				posteriorSum += currProb;
				
				currDemoStatePosterior[absHapIdx] = currProb;
			}
		}
		
		assert (Math.abs(posteriorSum - 1d) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
		
		return currLocusPosterior;
	}
	
	// returns [locus][demoState][haplotype]
	public ArrayList<ArrayList<double[]>> getPosteriorProbs() {
		assert (this.FwBwComputed);
		
		ArrayList<ArrayList<double[]>> posterior = new ArrayList<ArrayList<double[]>>();
		
		// now iterate over the loci
		for (int currDecodeIdx=0; currDecodeIdx<this.decodeLocusTransitionMap.numHmmSteps(); currDecodeIdx++) {
			posterior.add(getPosteriorProbs(currDecodeIdx));
		}
		
		return posterior;
	}
	
	//returns present config
	public DemoConfiguration<H> getPresentConfig(){
		return this.presentConfig;
	}
	
	// some private stuff
	public UberDemographyCore core;
	H mostRecentAdditionalHaplotype = null;

	// decode loci
	final HmmStepHandler<H> decodeLocusTransitionMap;

	// flags
	boolean FwBwComputed = false;
	private boolean expectationsComputed = false;
	
	// some global forward/backward things
	public TwoDimDoubleArray forwardLogProbabilities;
	public TwoDimDoubleArray backwardLogProbabilities;

	private double logLikelihood;
	
	// the parameter set
	public final EigenParamSet pSet;

	// locus specific stuff HMM stuff
	// initial [demoState] probabilities
	public double[] initialLogProbabilities;
	// some transitions
	public SMCSDTransitionList hmmProbsList;
	// and finally the emission [decodeLocusIdx,demoState,trunkType,emissionType]
	public double[][][][] logEmissionArray;
	
	//configuration
	private DemoConfiguration<H> presentConfig;
	
	private final ArrayList<double[]> trunkAbsorptionLogFractions;

	// this is to get the hidden states organized
	private int[] firstHapStateIdxMap;


}
