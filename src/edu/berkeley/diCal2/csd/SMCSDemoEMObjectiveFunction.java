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
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import edu.berkeley.diCal2.csd.HmmStepHandler.HasPartiallyMissing;
import edu.berkeley.diCal2.csd.HmmStepHandler.SingleStepMissingAllelesTransitionMap.CoreCache;
import edu.berkeley.diCal2.csd.auxiliary.MyEigenDecomposition;
import edu.berkeley.diCal2.haplotype.FSAHaplotype;
import edu.berkeley.diCal2.haplotype.GeneticType.GeneticTypeMultiplicity;
import edu.berkeley.diCal2.utility.LogSum;

public abstract class SMCSDemoEMObjectiveFunction<H extends FSAHaplotype> {

	public SMCSDemoEMObjectiveFunction (SMCSDemo<H> smcsd) {
		// remember core
		this.oldCore = smcsd.core;
		this.presentConfig= smcsd.getPresentConfig();
		// remember the transition things
		this.decodeLocusTransitionMap = smcsd.decodeLocusTransitionMap;
	}
	public DemoConfiguration<H> getPresentConfig(){
		return this.presentConfig;
	}
		
	public double value (CoreCache<H> demoCoreCache, double recoScaling, double[] theta){
		// build a new core
		UberDemographyCore newCore = demoCoreCache.getCore (this.presentConfig, oldCore.getObservedPresentDeme());
		
		// and get the value
		double toReturn = this.value(newCore, recoScaling, theta);
		
		return toReturn;
	}
	
	public abstract double value (UberDemographyCore newCore, double recoScaling, double[] theta);	

	public abstract double getLogLikelihood();

	private static boolean allZero(double[][] ds, double EPSILON) {
		for (double[] row : ds) {
			for (double value : row) {
				if (Math.abs(value) > EPSILON) return false;
			}
		}
		return true;
	}
	
	// pair class
	public static class EmissionExpectation {
		
		public final double[] mutRate;
		public final MyEigenDecomposition eigenMutation;
		// [demoState][trunkAllele][observedAllele]
		public final double[][][] expect;
	
		public EmissionExpectation (double[] mutRate, MyEigenDecomposition eigenMutation, double[][][] expect) {
			this.mutRate = mutRate;
			this.eigenMutation = eigenMutation;
			this.expect = expect;
		}
	}


	// the old core
	public final UberDemographyCore oldCore;
	protected HmmStepHandler<H> decodeLocusTransitionMap;
	protected DemoConfiguration<H> presentConfig;
	private Double logFactor = null;
	
	// if true, doesn't draw new demoStates for M-step
	// unfortunately fixing demo states doesn't work yet, because demo states are currently defined by integer index of refined interval/ancient deme,
	// which don't necessarly match up between different demographies
	public static class EmptyConditionalConfigOF<H extends FSAHaplotype> extends SMCSDemoEMObjectiveFunction<H> {
		
		public EmptyConditionalConfigOF (SMCSDemo<H> smcsd, H additionalHaplotype) {
			super (smcsd);
			
			//  calculate likelihood
			assert (!smcsd.FwBwComputed);
			// get it straight
			this.logLikelihood = smcsd.computeConditionalLogProb (additionalHaplotype);
		}
		
		@Override
		public double value (UberDemographyCore newCore, double recoScaling, double[] theta) {
			// it is just the log likelihood
			return logLikelihood;
		}
		
		@Override
		public double getLogLikelihood() {
			return logLikelihood;
		}

		public final double logLikelihood;
	}
	
	public static class LikelihoodOnlyObjectiveFunction<H extends FSAHaplotype> extends SMCSDemoEMObjectiveFunction<H> {

		public final double logLikelihood;
		
		public LikelihoodOnlyObjectiveFunction(SMCSDemo<H> smcsd, H additionalHaplotype) {
			super(smcsd);

			//  calculate likelihood
			assert (!smcsd.FwBwComputed);
			// get it straight
			this.logLikelihood = smcsd.computeConditionalLogProb (additionalHaplotype);
		}

		@Override
		public double value (UberDemographyCore newCore, double recoScaling, double[] theta) {
			// this should never be used
			assert (false);
			return Double.NEGATIVE_INFINITY;
		}

		@Override
		public double getLogLikelihood() {
			return this.logLikelihood;
		}
	}
	
	public static abstract class StructMigEmHMMObjectiveFunction<H extends FSAHaplotype> extends SMCSDemoEMObjectiveFunction<H> {
		
		public final ArrayList<EmissionExpectation> emissionExpectation;
		protected final int numDemoStates;
		public final double logLikelihood;
		
		@Override
		public double getLogLikelihood() {
			return this.logLikelihood;
		}


		public StructMigEmHMMObjectiveFunction(SMCSDemo<H> smcsd, H additionalHaplotype){
			super (smcsd);
						
			this.numDemoStates = this.oldCore.numDemoStates();

			// check whether forward/backward has been computed correctly, do so if not
			if (!smcsd.FwBwComputed || (additionalHaplotype != smcsd.mostRecentAdditionalHaplotype)) {
				// compute it
				smcsd.fillForwardBackward(additionalHaplotype);
			}
			
			// make sure (assume the haplotype is correct)
			assert (smcsd.forwardLogProbabilities.getFirstSize() == this.decodeLocusTransitionMap.numHmmSteps());
			assert (smcsd.backwardLogProbabilities.getFirstSize() == this.decodeLocusTransitionMap.numHmmSteps());
			
			// before we iterate over loci, let's get the likelihood of the sequence
			// by just summing over last entry in forward probabilities
			
			// and see how big this error gets
			this.logLikelihood = smcsd.getLogLikelihood();

			// expected emissions here, other expectation stuff somewhere else
			// figure out how many different mutation rates
			// and make a container for everybody
			Map<double[][][], EmissionExpectation> emissionMap = new LinkedHashMap<double[][][], EmissionExpectation>();
			
			// now iterate over the loci
			for (int currDecodeIdx=0; currDecodeIdx<this.decodeLocusTransitionMap.numHmmSteps(); currDecodeIdx++) {
				
				// get the right container to add the expected transition
				EmissionExpectation currPair = emissionMap.get(smcsd.logEmissionArray[currDecodeIdx]);

				// if not there, we make a new one
				if (currPair == null) {
					// make a new container
					currPair = new EmissionExpectation (this.decodeLocusTransitionMap.getMutationRate(currDecodeIdx), smcsd.pSet.getMutMatrixEigenDecomp(), new double[numDemoStates][smcsd.pSet.numAlleles()][smcsd.pSet.numAlleles()]);
					// and add it
					emissionMap.put (smcsd.logEmissionArray[currDecodeIdx], currPair);
				}
				
				double[][][] currEmissionExpectation = currPair.expect;
	
				
				// this should also sum to one
				double emissionSum = 0d;
				
				// now we have a valid container, so actually record the expected transition number added by transition from this locus to the next one
				// iterate absorption interval deme
				for (int absDemoState=0; absDemoState<this.oldCore.numDemoStates(); absDemoState++) {
					// iterate over absorption haplotype
					int absHapIdx = 0;
					int currDemeIdx = this.oldCore.getDemoStateCollection().getPresentDeme(absDemoState);
					HaplotypeConfiguration<H> currConfig = super.presentConfig.getPopulation(currDemeIdx);
					
					for (GeneticTypeMultiplicity<H> absHap : currConfig) {
								
						// now remember this emission's share of the expected emission
						double condProb = smcsd.getPosteriorProb (currDecodeIdx, absDemoState, absHapIdx);
						
						this.decodeLocusTransitionMap.updateEmissionTable(condProb, currEmissionExpectation[absDemoState], currDecodeIdx, absHap.geneticType, additionalHaplotype, smcsd.logEmissionArray[currDecodeIdx][absDemoState]);
						
						// and add to sum
						emissionSum += condProb;
						
						// increase idx
						absHapIdx++;
					}
				}
				
				// see what we got (might also not work with too many loci)
				assert (Math.abs(emissionSum-1d) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
				// go to next locus
			}
			
			// in the end, copy it
			this.emissionExpectation = new ArrayList<EmissionExpectation>();
			for (EmissionExpectation emissionExpect : emissionMap.values()) {
				emissionExpectation.add (emissionExpect);
			}
			
			// check whether everything is fine (doesn't work with a lot of loci) [seems to be ok actually]
			assert(checkValidEmissions());
		}
	
		private boolean checkValidEmissions(){
			// sum the emissions, should sum to numLoci
			double emissionSum = 0d;
			for (EmissionExpectation emissionExpect : this.emissionExpectation) {
				for (double[][] thisDemoState : emissionExpect.expect) {
					for (double[] row : thisDemoState) {
						for (double value : row) {
							emissionSum += value;
						}
					}
				}
			}
			// check it
			int numLoci = this.presentConfig.numLoci();
			int numNonEmisssionLoci = this.decodeLocusTransitionMap.numNonEmissionLoci();
			// this condition might be better implemented with an interface
			if (this.decodeLocusTransitionMap instanceof HasPartiallyMissing) {
				// we has to deal with partial missing loci
				// not very sophisticated, but for now ok
				int numPartiallyMissing = ((HasPartiallyMissing) this.decodeLocusTransitionMap).numPartiallyMissing();
				boolean ok = emissionSum - UberDemographyCore.MULTILOCUS_EPSILON(numLoci - numNonEmisssionLoci) <= numLoci - numNonEmisssionLoci;
				ok &= emissionSum + UberDemographyCore.MULTILOCUS_EPSILON(numLoci - numNonEmisssionLoci - numPartiallyMissing) >= numLoci - numNonEmisssionLoci - numPartiallyMissing;
				if (!ok) return false;
			}
			else {
				// just the normal thing
				int numEmissionLoci = numLoci- numNonEmisssionLoci;
				if (Math.abs (numEmissionLoci - emissionSum) > numEmissionLoci * UberDemographyCore.MULTILOCUS_EPSILON(numEmissionLoci)) return false;
			}
			return true;
		}
	
		protected abstract boolean checkValidMarginalsAndTransitions(int numDecodeLoci);
	
		protected double expectedLogEmissionLikelihood (UberDemographyCore newCore, double[] theta){
			int numDemoStates = newCore.numDemoStates();
			double value = 0d;
			for (EmissionExpectation mutStuff : this.emissionExpectation) {
				// get the emissison properbilities
				// [demoState,trunkType,emissionType]
				double[][][] logEmissionProb = newCore.getLogEmission (theta, mutStuff.eigenMutation);
				
				// now loop
				for (int demoState=0; demoState<numDemoStates; demoState++) {
					int presentDeme = newCore.getDemoStateCollection().getPresentDeme(demoState);
					
					// is it a ghost deme in the present?
					if (super.presentConfig.isGhost(presentDeme)) {
						// we should not observe it
						assert (allZero (mutStuff.expect[demoState], UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps())));
						continue;
					}
					
					for (int trunkAllele=0; trunkAllele< logEmissionProb[demoState].length; trunkAllele++) {
						for (int emissionAllele=0; emissionAllele < logEmissionProb[demoState][trunkAllele].length; emissionAllele++) {
							// add the right term
							double logProb = logEmissionProb[demoState][trunkAllele][emissionAllele];
							
							assert (mutStuff.expect[demoState][trunkAllele][emissionAllele] >= 0);
							// -Infinity * zero is zero
							if (mutStuff.expect[demoState][trunkAllele][emissionAllele] > UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps())) {
								value += logProb * mutStuff.expect[demoState][trunkAllele][emissionAllele];
							}
						}
					}
				}
			}
			
			return value;
		}
	}	
	
	// This objective function conditions on the absorbing lineage, but does not condition on whether transitions were recombinations or not
	public static class ConditionLineageOF<H extends FSAHaplotype> extends ConditionLineageTransitionTypeOF<H> {

		public ConditionLineageOF(SMCSDemo<H> smcsd, H additionalHaplotype) {
			super(smcsd, additionalHaplotype, true);
		}

		@Override
		protected double expectedLogNoRecombinationLikelihood (UberDemographyCore newCore, double[] logNoReco, double[][] logReco, TransitionExpectation transitionExpect){
			assert(countRecoToSelfAsNoReco);
			
			double value = 0d;
			
			// do the no reco part
			for (int demoState=0; demoState<numDemoStates; demoState++) {
				int presentDeme = newCore.getDemoStateCollection().getPresentDeme(demoState);
				
				// if it is a ghost interval deme
				if (this.presentConfig.isGhost(presentDeme)) {
					// we should not observe it
					assert (Math.abs(transitionExpect.noRecoExpect[demoState]) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
					continue;
				}
				
				// 0 * Inf = 0
				if (Math.abs(transitionExpect.noRecoExpect[demoState]) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps())) {
					continue;
				}
				
				// add the right term
				value += LogSum.computePairLogSum(logNoReco[demoState], logReco[demoState][demoState] - Math.log (newCore.getDemoStateCollection().getPresentHapConfig(this.presentConfig, demoState).totalGeneticTypes())) * transitionExpect.noRecoExpect[demoState];
			}
			
			return value;
		}		
	}
	
	// This objective function conditions on the absorbing lineage, and whether transitions were recombinations or not
	public static class ConditionLineageTransitionTypeOF<H extends FSAHaplotype> extends StructMigEmHMMObjectiveFunction<H> {

		public ConditionLineageTransitionTypeOF (SMCSDemo<H> smcsd, H additionalHaplotype) {
			this(smcsd, additionalHaplotype, false);
		}
		
		protected ConditionLineageTransitionTypeOF (SMCSDemo<H> smcsd, H additionalHaplotype, boolean countRecoToSelfAsNoReco) {
			super(smcsd, additionalHaplotype);
			
			this.countRecoToSelfAsNoReco = countRecoToSelfAsNoReco;
						
			// should sum to one
			double initSum = 0d;
			
			// first the expectation for the initial state
			initialExpectation = new double[numDemoStates];
			// iterate over interval demes
			for (int initDemoState=0; initDemoState<initialExpectation.length; initDemoState++) {
				// iterate over absorbing haplotype
				double demoStateProb = 0d;
				
				for (int initHap=0; initHap<this.oldCore.getDemoStateCollection().getPresentHapConfig(this.presentConfig, initDemoState).distinctGeneticTypes(); initHap++) {
					// get the marginal probability
					demoStateProb += smcsd.getPosteriorProb(0, initDemoState, initHap);
				}

				// and store it
				initialExpectation[initDemoState] += demoStateProb;
				
				// also update sum
				initSum += demoStateProb;
			}
			
			// look at sum
			// hacky
			assert (Math.abs(initSum-1d) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
			
			// now iterate over the loci
			for (int currDecodeIdx=0; currDecodeIdx<this.decodeLocusTransitionMap.numHmmSteps()-1; currDecodeIdx++) {
				
				// get some indices
				int nextDecodeIdx = currDecodeIdx+1;
				// and some transitions
				double[] noReco = smcsd.hmmProbsList.getLogNoReco (this.decodeLocusTransitionMap.getTransitionIndex(currDecodeIdx, false));
				double[][] reco = smcsd.hmmProbsList.getLogReco (this.decodeLocusTransitionMap.getTransitionIndex(currDecodeIdx, false));
				
				// what is the current transition index (no reverse direction)
				int transitionIndex = this.decodeLocusTransitionMap.getTransitionIndex (currDecodeIdx, false);
				
				// get the right container to add the expected transition
				TransitionExpectation thisTransitionExpectation = transitionExpectationMap.get (transitionIndex);
				// if not there, we make a new one
				if (thisTransitionExpectation == null) {
					// make a new container
					thisTransitionExpectation = new TransitionExpectation (new double[numDemoStates],new double[numDemoStates][numDemoStates]);
					// and add it
					transitionExpectationMap.put (transitionIndex, thisTransitionExpectation);
				}
				

				// here some things should also sum to one
				double transitionSum = 0d; 

				
				// reco
				
				// get a reference
				double[][] currRecoExpectation = thisTransitionExpectation.recoExpect;
			
				// now we have a valid container, so actually record the expected reco transition number added by transition from this locus to the next one
				// go through demoStates to precompute some variables
				double[] demoStateForwardLogProb = new double[numDemoStates];
				double[] demoStateBackwardLogStuff = new double[numDemoStates];
				// fill them
				for (int demoState=0; demoState<demoStateForwardLogProb.length; demoState++) {
					
					int currPresentDeme = this.oldCore.getDemoStateCollection().getPresentDeme(demoState);
					HaplotypeConfiguration<H> currConfig = this.oldCore.getDemoStateCollection().getPresentHapConfig(this.presentConfig, demoState);
					int numDistinctHaplotypes = currConfig.distinctGeneticTypes();
					
					// some logsum objects
					LogSum forwardLogProb = new LogSum (numDistinctHaplotypes);
					LogSum backwardLogStuff = new LogSum (numDistinctHaplotypes);
					// iterate over trunk haplotype
					int trunkHapIdx = 0;
					for (GeneticTypeMultiplicity<H> trunkHap : currConfig) {
						
						// add some forward stuff
						forwardLogProb.addLogSummand(smcsd.forwardLogProbabilities.get (currDecodeIdx, smcsd.getFirstHapStateIndex(demoState) +trunkHapIdx));
						
						// add some backward stuff
						double tmp = smcsd.backwardLogProbabilities.get (nextDecodeIdx, smcsd.getFirstHapStateIndex(demoState) + trunkHapIdx);
						tmp += this.decodeLocusTransitionMap.getLogEmissionProb(smcsd.logEmissionArray[nextDecodeIdx][demoState], nextDecodeIdx, trunkHap.geneticType, currPresentDeme, additionalHaplotype);
						tmp += smcsd.getLogTrunkAbsorptionFraction(currPresentDeme, trunkHapIdx);

						backwardLogStuff.addLogSummand(tmp);
							
						// increase idx
						trunkHapIdx++;
					}
					
					// get out the log sums
					demoStateForwardLogProb[demoState] = forwardLogProb.retrieveLogSum();
					demoStateBackwardLogStuff[demoState] = backwardLogStuff.retrieveLogSum();
				}
				
				
				// now fill the actual expectations
				// iterate over source deme
				for (int srcDemoState=0; srcDemoState<this.oldCore.numDemoStates(); srcDemoState++) {
					// iterate over destination deme
					for (int dstDemoState=0; dstDemoState<this.oldCore.numDemoStates(); dstDemoState++) {
											
						// calculate the actual joint probability of a recombination between these two demoStates
						double jointLogProb = demoStateForwardLogProb[srcDemoState];
						jointLogProb += demoStateBackwardLogStuff[dstDemoState];
						jointLogProb += reco[srcDemoState][dstDemoState];
						// normalize
						double condProb = Math.exp (jointLogProb - logLikelihood);
						
						// store it
						currRecoExpectation[srcDemoState][dstDemoState] += condProb;
												
						// and update the sum
						transitionSum += condProb;
					}
				}

				
				// no reco
				
				// get a reference
				double[] currNoRecoExpectation = thisTransitionExpectation.noRecoExpect;
							
				// now we have a valid container, so actually record the expected no reco transition number added by transition from this locus to the next one
				// iterate over interval deme (it's source and destination)
				for (int demoState=0; demoState<this.oldCore.numDemoStates(); demoState++) {
					
					int currPresentDeme = this.oldCore.getDemoStateCollection().getPresentDeme(demoState);
					HaplotypeConfiguration<H> currConfig = this.oldCore.getDemoStateCollection().getPresentHapConfig(this.presentConfig, demoState);
					
					// iterate over trunk haplotype
					int trunkHapIdx = 0;
					for (GeneticTypeMultiplicity<H> trunkHap : currConfig) {
								
						// now remember this transition's share of the expected transition
						// first the part that does not depend on the transition prob
						double tmp = smcsd.forwardLogProbabilities.get (currDecodeIdx, smcsd.getFirstHapStateIndex(demoState) + trunkHapIdx);
						tmp += smcsd.backwardLogProbabilities.get (nextDecodeIdx, smcsd.getFirstHapStateIndex(demoState) + trunkHapIdx);
						tmp += this.decodeLocusTransitionMap.getLogEmissionProb(smcsd.logEmissionArray[nextDecodeIdx][demoState], nextDecodeIdx, trunkHap.geneticType, currPresentDeme, additionalHaplotype);

						// now the no reco transition
						double transitionLogProb = noReco[demoState];
						
						// and multiply by it
						double jointLogProb = tmp + transitionLogProb;
						
						// divide by the log likelihood and store it
						double condProb = Math.exp(jointLogProb - logLikelihood);
						currNoRecoExpectation[demoState] += condProb;
						
						// and add to sum
						transitionSum +=  condProb;
						
						
						
						// check if we "count" a recombination to yourself as the same as a non-recombination
						if (this.countRecoToSelfAsNoReco) {
																				
							double selfRecoJointLogProb = tmp + reco[demoState][demoState] - Math.log(currConfig.totalGeneticTypes());
							
							double selfRecoCondProb = Math.exp(selfRecoJointLogProb - logLikelihood);
							currNoRecoExpectation[demoState] += selfRecoCondProb; 
							currRecoExpectation[demoState][demoState] -= selfRecoCondProb;
						}
						
						
						// increase idx
						trunkHapIdx++;
					}
				}
				
				// see what we got
				// hacky
				assert (Math.abs(transitionSum-1d) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
				// go to next locus
			}
			
			
			// hacky
			assert(checkValidMarginalsAndTransitions(this.decodeLocusTransitionMap.numHmmSteps()));
			
		}
		
		
		@Override
		public boolean checkValidMarginalsAndTransitions(int numDecodeLoci) {
			
			// sum the initial guys, should just sum to one
			double initSum = 0d;
			for (double value : initialExpectation) {
				initSum += value;
			}
			// check it
			if (Math.abs (1d - initSum) > UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps())) {
				assert (false);
				return false;
			}
			
			// sum the transition guys, should sum to numLoci-1
			double transitionSum = 0d;
			// no reco
			for (TransitionExpectation transitionExpect : transitionExpectationMap.values()) {
				for (double value : transitionExpect.noRecoExpect) {
					transitionSum += value;
				}
				
				for (double[] row : transitionExpect.recoExpect) {
					for (double value : row){
						transitionSum += value;
					}
				}
			}
			// check it
			if (Math.abs (numDecodeLoci - 1d - transitionSum) > (numDecodeLoci - 1d) * UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps())) {
				assert (false);
				return false;
			}
			
	
			// everything fine
			return true;
		}
		
		
		protected double expectedLogInitialLikelihood (UberDemographyCore newCore){
			int numDemoStates = newCore.numDemoStates();

			double value = 0d;

			assert (this.initialExpectation.length == numDemoStates);
			assert (numDemoStates == newCore.getDemoStateLogMarginal().length);

			// iterate over initial demoState
			for (int demoState=0; demoState<numDemoStates; demoState++) {
				
				int currDemeIdx = newCore.getDemoStateCollection().getPresentDeme(demoState);
				HaplotypeConfiguration<H> currDeme = newCore.getDemoStateCollection().getPresentHapConfig(this.presentConfig, demoState);
				
				// if it is a ghost interval deme
				if (this.presentConfig.isGhost (currDemeIdx)) {
					// we should not observe it
					assert (Math.abs(this.initialExpectation[demoState]) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
					continue;
				}

				// compute the term using the new core (initial probabilities)
				double logAbsProb = newCore.getDemoStateLogMarginal()[demoState] - Math.log (currDeme.totalGeneticTypes());
				// and add it
				// 0 * Inf = 0
				if (this.initialExpectation[demoState] >= UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps())){
					value += logAbsProb * this.initialExpectation[demoState];
				}
			}
			
			return value;
		}
		
		protected double expectedLogRecombinationLikelihood (UberDemographyCore newCore, double[][] logRecoTransition, TransitionExpectation transitionExpect){
			
			double value = 0d;
			
			// do the reco part
			for (int srcDemoState=0; srcDemoState<numDemoStates; srcDemoState++) {
				for (int dstDemoState=0; dstDemoState<numDemoStates; dstDemoState++) {
	
					int srcPresentDeme = newCore.getDemoStateCollection().getPresentDeme(srcDemoState);
					int dstPresentDeme = newCore.getDemoStateCollection().getPresentDeme(dstDemoState);
					
					// if it is a ghost interval deme
					if (this.presentConfig.isGhost(srcPresentDeme) || this.presentConfig.isGhost (dstPresentDeme)) {
						// we should not observe it
						assert (Math.abs(transitionExpect.recoExpect[srcDemoState][dstDemoState]) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
						continue;
					}
					
					// should be nonnegative
					assert (transitionExpect.recoExpect[srcDemoState][dstDemoState] >= -UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
					
					
					// 0 * Inf = 0				
					if (transitionExpect.recoExpect[srcDemoState][dstDemoState] < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps())) {
						continue;
					}
					
					// get the reco transition log prob (with destination normalizer)
					double logTransitionProb = logRecoTransition[srcDemoState][dstDemoState] - Math.log (newCore.getDemoStateCollection().getPresentHapConfig(this.presentConfig, dstDemoState).totalGeneticTypes());
					value += logTransitionProb * transitionExpect.recoExpect[srcDemoState][dstDemoState];					
				}
			}
			
			// return it
			return value;
		}
		
		protected double expectedLogNoRecombinationLikelihood (UberDemographyCore newCore, double[] logNoRecoTransition, double[][] logRecoTransition, TransitionExpectation transitionExpect){
			assert(!countRecoToSelfAsNoReco);
			
			double value = 0d;
			
			// do the no reco part
			for (int demoState=0; demoState<numDemoStates; demoState++) {
				// if it is a ghost interval deme
				if (this.presentConfig.isGhost (newCore.getDemoStateCollection().getPresentDeme(demoState))) {
					// we should not observe it
					assert (Math.abs(transitionExpect.noRecoExpect[demoState]) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps()));
					continue;
				}
				
				if (Math.abs(transitionExpect.noRecoExpect[demoState]) < UberDemographyCore.MULTILOCUS_EPSILON(this.decodeLocusTransitionMap.numHmmSteps())) {
					continue;
				}
				
				// add the right term
				value += logNoRecoTransition[demoState] * transitionExpect.noRecoExpect[demoState];
			}
			
			return value;
		}
		
		@Override
		public double value (UberDemographyCore newCore, double recoScaling, double[] theta) {
			
			// just some assertion
			assert (newCore.numDemoStates() == this.oldCore.numDemoStates());
			
			// the return value
			double value = 0d;
			
			// the initial distribution
			value += expectedLogInitialLikelihood (newCore);

			// get a new list
			this.decodeLocusTransitionMap.updateTheta(theta);
			
			SMCSDTransitionList newTransitionList = this.decodeLocusTransitionMap.computeTransitions(newCore, recoScaling);

			// and go through stuff
			for (Entry<Integer, TransitionExpectation> transitionExpect : this.transitionExpectationMap.entrySet()) {

				// get the new transitions probabilities
				double[] logNoRecoTransition = newTransitionList.getLogNoReco (transitionExpect.getKey());
				double[][] logRecoTransition = newTransitionList.getLogReco (transitionExpect.getKey());

				value += expectedLogNoRecombinationLikelihood (newCore, logNoRecoTransition, logRecoTransition, transitionExpect.getValue());
			
				// reco part
				value += expectedLogRecombinationLikelihood (newCore, logRecoTransition, transitionExpect.getValue());
			}
			
			// the emissions
			value += expectedLogEmissionLikelihood(newCore, theta);
			
			// return it
			return value;
		}
		
		public static class TransitionExpectation{
			// [demoState]
			public final double[] noRecoExpect;
			// [srcDemoState][dstDemoState]
			public final double[][] recoExpect;
			
			public TransitionExpectation (double[] noRecoExpect, double[][] recoExpect) {
				this.noRecoExpect = noRecoExpect;
				this.recoExpect = recoExpect;
			}			
		}
		
		// here are the expectations
		public final double[] initialExpectation;
		public final TreeMap<Integer,TransitionExpectation> transitionExpectationMap = new TreeMap<Integer, TransitionExpectation>();
		
		protected final boolean countRecoToSelfAsNoReco;
		
	}
	
	public static class MarginalKLDivergence<H extends FSAHaplotype> extends ConditionLineageTransitionTypeOF<H> {
		
		public MarginalKLDivergence(SMCSDemo<H> smcsd, H additionalHaplotype) {
			super(smcsd, additionalHaplotype);
			double[] averageMarginal = new double[numDemoStates];
			for (int i = 0; i < numDemoStates; i++) {
				averageMarginal[i] = initialExpectation[i];
			}
			
			for (TransitionExpectation transitionExp : transitionExpectationMap.values()) {
				for (int dstDemoState = 0; dstDemoState < numDemoStates; dstDemoState++) {
					
					averageMarginal[dstDemoState] += transitionExp.noRecoExpect[dstDemoState];
					
					for (int srcDemoState = 0; srcDemoState < numDemoStates; srcDemoState++) {
						averageMarginal[dstDemoState] += transitionExp.recoExpect[srcDemoState][dstDemoState];
					}
				}
			}
			
			// replace initial expectation
			for (int i = 0; i < numDemoStates; i++) {
				initialExpectation[i] = averageMarginal[i];
			}
		}

		@Override
		public double value(UberDemographyCore newCore, double recoScaling, double[] theta) {
			return expectedLogInitialLikelihood (newCore) + expectedLogEmissionLikelihood(newCore, theta);
		}

		@Override
		public double getLogLikelihood() {
			//TODO: Fix this to work with the theta estimation
			assert(false);
			return this.value(this.oldCore, 1d, null);
		}
	}

	public Double getLogFactor() {
		return this.logFactor;
	}
	
	public void setLogFactor (Double f) {
		this.logFactor = f;
	}
}
