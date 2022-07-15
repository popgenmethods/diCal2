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

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import edu.berkeley.diCal2.csd.DemoState.DemoStateCollection;
import edu.berkeley.diCal2.csd.TrunkProcess.SimpleTrunk;
import edu.berkeley.diCal2.csd.auxiliary.MyEigenDecomposition;
import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.utility.LogSum;
import edu.berkeley.diCal2.utility.RealPartition.Interval;
import edu.berkeley.diCal2.utility.SimplePair;
import edu.berkeley.diCal2.utility.SumArray;

public abstract class UberDemographyCore {
	
	public UberDemographyCore(DemoStateCollection demoStateCollection, boolean cacheTransitionEmissions, double renormalizeEpsilon) {
		super();
		this.demoStateCollection = demoStateCollection;
		this.recoCache = cacheTransitionEmissions ? new TreeMap<Double, RecoLogProbs>() : null;
		this.logEmissionCache = cacheTransitionEmissions ? new HashMap<SimplePair<Double, MyEigenDecomposition>, double[][][]>() : null;
		this.renormalizeEpsilon = renormalizeEpsilon;
	}

	// some of them wrappers
	public static UberDemographyCore getDemographyCore  (int[] sampleSizes, int observedPresentDeme, Demography demo, Interval[] hiddenStateIntervals, boolean cacheTransitionEmissions, boolean useEigenCore) {
		return getDemographyCore (sampleSizes, new DemoStateCollection(demo, hiddenStateIntervals, null, UberDemographyCore.ONELOCUS_EPSILON), observedPresentDeme, cacheTransitionEmissions, useEigenCore);
	}
	
	public static UberDemographyCore getDemographyCore (int[] sampleSizes, DemoStateCollection demoStates, int observedPresentDeme, boolean cacheTransitionEmissions, boolean useEigenCore) {
		return getDemographyCore (new SimpleTrunk(sampleSizes, demoStates.refinedDemography), demoStates, observedPresentDeme, cacheTransitionEmissions, UberDemographyCore.DEFAULT_RENORMALIZATION_EPSILON, useEigenCore);
	}
	
	public static UberDemographyCore getDemographyCore (TrunkProcess trunk, DemoStateCollection demoStates, int observedPresentDeme, boolean cacheTransitionEmissions, boolean useEigenCore) {
		return getDemographyCore(trunk, demoStates, observedPresentDeme, cacheTransitionEmissions, UberDemographyCore.DEFAULT_RENORMALIZATION_EPSILON, useEigenCore);
	}
	
	// the real thing
	public static UberDemographyCore getDemographyCore (TrunkProcess trunk, DemoStateCollection demoStates, int observedPresentDeme, boolean cacheTransitionEmissions, double renormalizeEpsilon, boolean useEigenCore) {
		// for the moment we decide between two possible cores
		// good old eigenstuff core, if we do not have some exponential growth rates given
		// WARNING base case (empty cond. sample) has trunk == null, in this case the type of core does not even matter 

		
//		if ((trunk == null) || (trunk.refinedDemography.expGrowthRates == null)) {
//			return new EigenCore (trunk, demoStates, observedPresentDeme, cacheTransitionEmissions, renormalizeEpsilon);
//		}
//		// new core with ode solver, if we have exponential rates
//		else {
//			//TODO make this an option
//			return new ODECore (trunk, demoStates, observedPresentDeme, cacheTransitionEmissions, renormalizeEpsilon, false);
//		}
		if (useEigenCore) {
			// this one can't do exponential rates
			// assertion here, proper input checking when reading params
			assert (trunk.refinedDemography.expGrowthRates == null);

			// create core
			return new EigenCore (trunk, demoStates, observedPresentDeme, cacheTransitionEmissions, renormalizeEpsilon);
		}
		// new core with ode solver, if we have exponential rates
		else {
			
			if (trunk.refinedDemography.hasPulseMigration()) {
				// we have to implement this properly. Not done yet.
				throw new RuntimeException ("ODE core cannot be used if the demographic model has pulse migration/admixture. Use: --useEigenCore instead.");
			}
			// create core
			return new ODECore (trunk, demoStates, observedPresentDeme, cacheTransitionEmissions, renormalizeEpsilon, false);
		}

	}
	
	
	protected double[] getRefinedTheta(double[] mutationRate){
		List<Double> refinedMutationRate= new ArrayList<Double>();
		if(mutationRate.length == 1){
			double[] newMut= new double[this.getDemoStateCollection().getNumHiddenTimes()];
			for(int i= 0; i < newMut.length; i++){
				newMut[i]= mutationRate[0];
			}
			mutationRate= newMut;
		}
		
		for(int i= 0; i < mutationRate.length; i++){
			for(int j= 0; j < this.getDemoStateCollection().getNumIntervalsInHiddenTime(i); j++){
				refinedMutationRate.add(mutationRate[i]);
			}
		}
		
		double[] toReturn = new double[refinedMutationRate.size()];
		for(int i = 0; i < toReturn.length; i++){
			toReturn[i] = refinedMutationRate.get(i);
		}
		return toReturn;
	}
	
	// some sub-classes 
	public static class RecoLogProbs {
		final double[][] logReco;
		final double[] logNoReco;
		
		public RecoLogProbs(double[][] logReco, double[] logNoReco) {
			super();
			this.logReco = logReco;
			this.logNoReco = logNoReco;
		}
	}

	public DemoStateCollection getDemoStateCollection() {
		return demoStateCollection;
	}

	private final DemoStateCollection demoStateCollection;
	
	
	// some epsilon stuff
	/// the precision
	public final static double ONELOCUS_EPSILON = 1e-12;
	
	// handles multilocus, and multi-multilocus, epsilons
	public static double MULTILOCUS_EPSILON(int numDecodeLoci) {
		// we use 100 here, because empirically we saw that we have too
		if (numDecodeLoci < 100) {
			numDecodeLoci = 100;
		}
		return Math.sqrt(numDecodeLoci*Math.log(Math.log(numDecodeLoci))) * 1e-7;
	}
	
	/// we allow this relative error margin
	/// the reference is always the previously computed value
	public final static double RELATIVE_ERROR = 0.0001;
	public final static double DEMO_FACTORY_EPSILON = 1e-4;
	
	public final static double DEFAULT_RENORMALIZATION_EPSILON = 1e-5;
	protected final double renormalizeEpsilon;
	
	// this is apparently used
	// it being static is not too nice, but ok for now
	// this being public is actually the worst thing ever
	public static boolean RENORMALIZE_PROBS = true;
	
	
	// here comes all the stuff that the derived classes have to implement
	public abstract boolean isConditionalConfigEmpty();
	public abstract int getObservedPresentDeme();
	public abstract TrunkProcess getTrunk();
	
	public abstract int numDemoStates();
	public abstract Interval getTimeInterval(int intervalDeme);
	
	public abstract double[] getDemoStateLogMarginal();
	
	private final TreeMap<Double, RecoLogProbs> recoCache;
	
	synchronized public RecoLogProbs getRecoLogProbs (double recoRate) {
		if (recoCache == null) return this.computeRecoLogProbs(recoRate);
		
		RecoLogProbs toReturn = recoCache.get(recoRate);
		if (toReturn == null) {
			toReturn = this.computeRecoLogProbs(recoRate);
			recoCache.put(recoRate, toReturn);
		}
		return toReturn;
	}
	protected abstract RecoLogProbs computeRecoLogProbs (double recoRate);
	// [demoState][trunkType][emissionType]
	
	private final Map<SimplePair<Double, MyEigenDecomposition>, double[][][]> logEmissionCache;
	
	synchronized public double[][][] getLogEmission (double[] mutRate, MyEigenDecomposition firstEigenMutMatrix) {
		return this.computeLogEmission (mutRate, firstEigenMutMatrix);
	}

	protected abstract double[][][] computeLogEmission (double[] mutRate, MyEigenDecomposition firstEigenMutMatrix);
	
	// returns prior probability of [allele]
	// given counts for each [presentDeme][allele]
	// log emission probabilities for [demoState][fromAllele][toAllele]
	public double[] imputeLocus (int[][] observedAlleleCounts, double[][][] logEmissionProb) {
		
		int numAlleles = observedAlleleCounts[0].length;
		int numDemes = this.getTrunk().getSampleSizes().length;
		
		double[] imputedProbs = new double[numAlleles];
		
		// get fraction of alleles in each deme
		double[][] observedAlleleFractions = new double[numDemes][numAlleles];
		for (int deme = 0; deme < numDemes; deme++) {
			int totalCount = SumArray.getSum(observedAlleleCounts[deme]);
			assert totalCount == this.getTrunk().getSampleSizes()[deme] : "Missing trunk alleles not handled for this imputation.";
			
			for (int allele = 0; allele < numAlleles; allele++) {
				if (totalCount != 0d) observedAlleleFractions[deme][allele] = (double) observedAlleleCounts[deme][allele] / (double) totalCount;
			}

		}
		
		// integrate over all hidden states and alleles
		for (int toAllele = 0; toAllele < numAlleles; toAllele++) {
			for (int demoState = 0; demoState < this.numDemoStates(); demoState++ ) {
				for (int fromAllele = 0; fromAllele < numAlleles; fromAllele++) {
					imputedProbs[toAllele] += observedAlleleFractions[this.getDemoStateCollection().getPresentDeme(demoState)][fromAllele] * Math.exp(this.getDemoStateLogMarginal()[demoState] + logEmissionProb[demoState][fromAllele][toAllele]);
				}		
			}
		}
		
		return imputedProbs;
	}
	
	
	public abstract void dump(double recoRate, double[] mutRate, double[][] mutMatrix, PrintStream out);
	public abstract void checkAncientEmissions (double[] mutationRate, MyEigenDecomposition eigenMutation);
	public abstract void checkAncientDemeProbs(double recoRate);
	public abstract void checkAncientMarginal();
	
	// makes sure rows of matrix sum to 1
	public static void renormalizeStochasticMatrix(double[][] stochasticMatrix, double renormalizeEpsilon) {
		int numRows = stochasticMatrix.length;
		assert (numRows > 0);
		
		for (int row=0; row < numRows; row++) {
			int numCols = stochasticMatrix[row].length;
			assert (numRows == numCols);
			
			double sum = 0d;
			for (int column = 0; column < stochasticMatrix[row].length; column++) {
				
				// assert that the value is not too negative, so ideally it should be zero
				assert (stochasticMatrix[row][column] > -renormalizeEpsilon);
				// set it to zero if it is morally zero
				if (stochasticMatrix[row][column] < 0) stochasticMatrix[row][column] = 0;
				
				sum += stochasticMatrix[row][column];
			}
			
			assert (Math.abs(1d - sum) < renormalizeEpsilon);
			
			for (int column = 0; column < numCols; column++) {
				stochasticMatrix[row][column] /= sum;
			}
		}
	}

	public static void renormalizeLogStochasticMatrix(double[][] logStochasticMatrix, double renormalizeEpsilon) {
		int numRows = logStochasticMatrix.length;
		assert (numRows > 0);
		
		for (int row=0; row < numRows; row++) {
			int numCols = logStochasticMatrix[row].length;
			assert (numRows == numCols);
			
			renormalizeLogStochasticVector(logStochasticMatrix[row], renormalizeEpsilon);
		}
	}
	
	public static void renormalizeLogTransitions(RecoLogProbs recomLogProbs, double[] logMarginal, double renormalizeEpsilon) {
		double[][] logReco = recomLogProbs.logReco;
		double[] logNoReco = recomLogProbs.logNoReco;
		
		int numRows = logReco.length;
		assert(logNoReco.length == numRows);
		
		for (int row = 0; row < numRows; row++) {
			assert (logReco[row].length == numRows);
			
			
			LogSum logSum = new LogSum(numRows + 1);
			for (double value : logReco[row]) logSum.addLogSummand(value);
			logSum.addLogSummand(logNoReco[row]);
			
			double logSumValue = logSum.retrieveLogSum();
			if (logMarginal[row] == Double.NEGATIVE_INFINITY){
				assert (logSumValue == Double.NEGATIVE_INFINITY);
				continue;
			}
			
			for (int col = 0; col < logReco[row].length; col++) {
				logReco[row][col] -= logSumValue;
			}
			logNoReco[row] -= logSumValue;
		}
	}
	
	public static void renormalizeLogStochasticVector(double[] logStochasticVector, double renormalizeEpsilon) {
		
		LogSum logSum = new LogSum(logStochasticVector.length);
		for (double value : logStochasticVector) logSum.addLogSummand(value);
		
		double logSumValue = logSum.retrieveLogSum();
		if (Math.abs(1d - Math.exp(logSumValue)) > renormalizeEpsilon){
			// TODO: cleaner solution
//			System.out.println("# [ASSERTION_ERROR] The probabilities don't sum to one, but: " + Math.exp(logSumValue) + ". Trying to proceed anyways and renormalize it to one.");
			assert false;
		}
		
		for (int i = 0; i < logStochasticVector.length; i++) {
			logStochasticVector[i] -= logSumValue;
		}
	}

}
