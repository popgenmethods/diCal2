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
import java.util.Arrays;

import Jama.Matrix;
import edu.berkeley.diCal2.csd.DemoState.DemoStateCollection;
import edu.berkeley.diCal2.csd.TrunkProcess.SimpleTrunk;
import edu.berkeley.diCal2.csd.auxiliary.MyEigenDecomposition;
import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.haplotype.FSAHaplotype;
import edu.berkeley.diCal2.utility.LogSum;
import edu.berkeley.diCal2.utility.RateMatrixTools;
import edu.berkeley.diCal2.utility.RealPartition.Interval;
import edu.berkeley.diCal2.utility.SimplePair;
import edu.berkeley.diCal2.utility.SumArray;

public class ODECore extends UberDemographyCore{

	public ODECore (int[] sampleSizes, int observedPresentDeme, Demography demo, Interval[] hiddenStateIntervals) {
		this (sampleSizes, new DemoStateCollection(demo, hiddenStateIntervals, null, UberDemographyCore.ONELOCUS_EPSILON), observedPresentDeme);
	}
	
	public ODECore (int[] sampleSizes, DemoStateCollection demoStates, int observedPresentDeme) {
		this (new SimpleTrunk(sampleSizes, demoStates.refinedDemography), demoStates, observedPresentDeme, false);
	}
	
	public ODECore (TrunkProcess trunk, DemoStateCollection demoStates, int observedPresentDeme, boolean cacheRecoTransitions) {
		this(trunk, demoStates, observedPresentDeme, cacheRecoTransitions, UberDemographyCore.DEFAULT_RENORMALIZATION_EPSILON, false);
	}
	
	/// the constructor, gets configuration, parameters and number of intervals
	public ODECore (TrunkProcess trunk, DemoStateCollection demoStates, int observedPresentDeme, boolean cacheRecoTransitions, double renormalizeEpsilon, boolean smcPrime) {
		
		super(demoStates, cacheRecoTransitions, renormalizeEpsilon);
		
		// remember some things
		this.smcPrime = smcPrime;
		this.trunk = trunk;
		this.observedPresentDeme = observedPresentDeme;

		// is it empty?
		this.conditionalConfigEmpty = trunk == null ? true : (SumArray.getSum(trunk.getSampleSizes()) == 0);
		
		// again, if the conditional config is empty, we don't have to do anything
		if (!this.conditionalConfigEmpty) {

			this.refinedDemography = this.trunk.getRefinedDemography();
			assert (this.refinedDemography == this.getDemoStateCollection().refinedDemography);
			
			// locus independent part of the HMM
			// precompute the non-absoprtion and absorption probabilities (need also in ODE mode, should work)
			this.demoTransitions = preComputeDemoTransitions (this.observedPresentDeme, this.trunk, this.renormalizeEpsilon);

			// initialize the marginal intervalDeme probabilities (ancient ones)
			this.epochAncientDemeMarginalLogProbabilities = computeEpochAncientDemeLogMarginal (this.observedPresentDeme, this.trunk, this.demoTransitions.Q);

			for (double[] thisEpochMarginalProbs : epochAncientDemeMarginalLogProbabilities){
				for (double value : thisEpochMarginalProbs){
					assert (Math.exp(value) > 0d - ONELOCUS_EPSILON);
				}
			}
			
			
			// initialize the marginal intervalDeme probabilities (real ones)
			this.demoStateLogMarginal = computeDemoStateLogMarginal ();

			for (double value : this.demoStateLogMarginal){
				assert (Math.exp(value) > 0d - ONELOCUS_EPSILON);
			}
		}
		else {
			
			assert(this.trunk == null);
			assert (this.getDemoStateCollection() == null);
			
			// blank
			this.demoTransitions = null;
			this.epochAncientDemeMarginalLogProbabilities = null;
			this.demoStateLogMarginal = null;
			this.refinedDemography = null;
		}	
	}
		
	private static <H extends FSAHaplotype> DemographyTransitions preComputeDemoTransitions (int observedPresentDeme, TrunkProcess trunk, double renormalizeEpsilon) {
 
		// for each interval we have a matrix (includes the absorbing states)
		// [interval][starting state][ending state]
		ArrayList<double[][]> f = new ArrayList<double[][]>();
		
		Demography refinedDemography = trunk.getRefinedDemography();
		
		// go through intervals
		for (int e=0; e<refinedDemography.numIntervals(); e++) {
			int numAncientDemes = refinedDemography.numAncientDemes(e);
			double[][] currF;
			// is pulse?
			if (refinedDemography.isPulse(e)) {
				// its a pulse!
				// initialize all to zero
				currF = new double[2*numAncientDemes][2*numAncientDemes];
				// just copy the pulse matrix for the non-absorbing demes
				double[][] pulseMatrix = refinedDemography.pulseMigrationMatrixList.get(e);
				for (int i=0; i<numAncientDemes; i++) {
					for (int j=0; j<numAncientDemes; j++) {
						currF[i][j] = pulseMatrix[i][j];
					}
					// and fill in some 1s on the diagonal to make it a stochastic matrix
					currF[numAncientDemes+i][numAncientDemes+i] = 1d;
				}
			}
			// not pulse, do more
			else {
				// compute the transition probabilities using the ODE
				// make sure it really is not a pulse
				assert (refinedDemography.migrationMatrixList.get(e) != null);
				assert (refinedDemography.popSizes.get(e) != null);
				assert (refinedDemography.expGrowthRates.get(e) != null);
				
				// the use of the absorption rates seems to be fine for now
				currF = CoalescentODEs.getMarginalTransitionMatrix (trunk.refinedDemography.epochList[e], trunk.refinedDemography.migrationMatrixList.get(e), trunk.getAbsorbRates(e), trunk.refinedDemography.expGrowthRates.get(e));
				
				// and check that the matrix has the right dimensions
				// now fill the matrix of non-absorption probabilities
				assert (currF.length == 2*numAncientDemes);
				assert (RateMatrixTools.isSquare(currF));
			}
			
			// renormalize the f's so every row sums to one (they should do so anyways, but numerical stability yada yada)
			if (UberDemographyCore.RENORMALIZE_PROBS) renormalizeStochasticMatrix(currF, renormalizeEpsilon);
			assert (RateMatrixTools.isProperStochaticMatrix (currF, UberDemographyCore.ONELOCUS_EPSILON));
			
			f.add (currF);
		}
		// so these matrices should be fine now
		
		
		int numIntervals = refinedDemography.numIntervals();
		// now we calculate some probabilities for not being absorbed for a long time
		// [first interval][second interval][first start deme][second start deme]
		ArrayList<ArrayList<double[][]>> p = new ArrayList<ArrayList<double[][]>>();
		
		// loop over first interval
		for (int firstInterval=0; firstInterval<numIntervals; firstInterval++) {
			// inner container
			ArrayList<double[][]> currP = new ArrayList<double[][]>();
			// loop over second interval
			for (int secondInterval=0; secondInterval<numIntervals; secondInterval++) {
				if (secondInterval < firstInterval) {
					// nada
					currP.add(null);
				}
				else if (firstInterval == secondInterval) {
					// initialize the dp
					currP.add(Matrix.identity(refinedDemography.numAncientDemes(firstInterval), refinedDemography.numAncientDemes(firstInterval)).getArray());
				}
				else {
					// step forwward to next matrix
					// last in list should be the previous matrix
					double[][] prevMatrix = currP.get(currP.size()-1);
					assert (prevMatrix.length == refinedDemography.numAncientDemes(firstInterval));
					assert (prevMatrix[0].length == refinedDemography.numAncientDemes(secondInterval-1));
					// now step forward
					double[][] newMatrix = new double[refinedDemography.numAncientDemes(firstInterval)][refinedDemography.numAncientDemes(secondInterval)];
					// fill it
					// loop over first starting deme
					for (int firstAncientDeme=0; firstAncientDeme<newMatrix.length; firstAncientDeme++) {
						// loop over second starting deme
						for (int secondAncientDeme=0; secondAncientDeme<newMatrix[0].length; secondAncientDeme++) {
							// initialize
							newMatrix[firstAncientDeme][secondAncientDeme] = 0d;
							for (int memberAncientDeme : refinedDemography.getMemberDemesIndices (secondAncientDeme, secondInterval)) {
								assert (memberAncientDeme <= refinedDemography.numAncientDemes(secondInterval-1));
								// loop over intermediate deme
								for (int intermediateAncientDeme=0; intermediateAncientDeme<prevMatrix[0].length; intermediateAncientDeme++) {
									newMatrix[firstAncientDeme][secondAncientDeme] += prevMatrix[firstAncientDeme][intermediateAncientDeme] * f.get(secondInterval-1)[intermediateAncientDeme][memberAncientDeme];
								}
							}
						}
					}
					// remember the matrix
					currP.add (newMatrix);
				}
			}
			// remember this P
			p.add (currP);
		}
		
		
		// now get some Q's
		// [first interval][absorbing interval][first start deme][absorbing deme]
		ArrayList<ArrayList<double[][]>> Q = new ArrayList<ArrayList<double[][]>>();
		// loop over first interval
		for (int firstInterval=0; firstInterval<numIntervals; firstInterval++) {
			// inner container
			ArrayList<double[][]> currQ = new ArrayList<double[][]>();
			// loop over second interval
			for (int absorbingInterval=0; absorbingInterval<numIntervals; absorbingInterval++) {
				if (absorbingInterval < firstInterval) {
					// nada
					currQ.add(null);
				}
				else {
					// what's the probablility of absorption
					// we don't need to treat i == j specially, cause we have the identity matrix
					double[][] newMatrix = new double[refinedDemography.numAncientDemes(firstInterval)][refinedDemography.numAncientDemes(absorbingInterval)];
					// fill it
					// loop over first starting deme
					for (int startAncientDeme=0; startAncientDeme<newMatrix.length; startAncientDeme++) {
						// loop over second starting deme
						for (int absorbingAncientDeme=0; absorbingAncientDeme<newMatrix[0].length; absorbingAncientDeme++) {
							// initialize
							newMatrix[startAncientDeme][absorbingAncientDeme] = 0d;
							// starting deme in the absorbing interval
							for (int intermediateAncientDeme=0; intermediateAncientDeme<refinedDemography.numAncientDemes(absorbingInterval); intermediateAncientDeme++) {
								// we want to go from startDeme to intermediate deme to absorbing deme
								// remember the indexing for the absorbing states in the extended matrices
								newMatrix[startAncientDeme][absorbingAncientDeme] += p.get(firstInterval).get(absorbingInterval)[startAncientDeme][intermediateAncientDeme] * f.get(absorbingInterval)[intermediateAncientDeme][refinedDemography.numAncientDemes(absorbingInterval) + absorbingAncientDeme];
							}
						}
					}
					
					// remember the matrix
					currQ.add (newMatrix);
				}
			}
			
			// remember thisQ
			Q.add (currQ);
		}

		// then return p and Q
		return new DemographyTransitions(p, Q, f);
	}


	
	// [epoch][ancient deme]
	// returns the probability of being absorbed into given (epoch,ancient deme)
	private static <H extends FSAHaplotype> ArrayList<double[]> computeEpochAncientDemeLogMarginal (int observedPresentDeme, TrunkProcess trunk, ArrayList<ArrayList<double[][]>> q ) {

		// fill the initial probabilities from q
		// it's just the aborption stuff when starting in the first interval
		ArrayList<double[]> initialProb = new ArrayList<double[]>();
		
		// go through epochs
		for (int epoch = 0; epoch < trunk.refinedDemography.epochList.length; epoch++){
			
			//make an array big enough
			double[] thisEpochLogMarginal = new double[trunk.refinedDemography.treePartitionList.get(epoch).size()];
			// add it to array list
			initialProb.add(thisEpochLogMarginal);
		
			// go through them ancient demes
			for (int ancientDeme=0; ancientDeme<thisEpochLogMarginal.length; ancientDeme++ ){
				// scary
				if (!trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
					double preU = q.get(0).get(epoch)[observedPresentDeme][ancientDeme];
					// should be positive
					assert (preU > 0d - ONELOCUS_EPSILON);
					// be nice
					if (preU < 0) preU = 0d;
	
					thisEpochLogMarginal[ancientDeme] = Math.log (preU);
				}
				else {
					thisEpochLogMarginal[ancientDeme] = Double.NEGATIVE_INFINITY;
				}
				
				assert(!Double.isNaN(thisEpochLogMarginal[ancientDeme]));
			}
		}
		// return it
		return initialProb;
	}

	// [demoState]
	// returns probability of being absorbed into the demoState (summed over all states living in the demoState)
	private double[] computeDemoStateLogMarginal() {
		
		double[] demoStateLogMarginal = new double[this.numDemoStates()];
		// fill with log 0s
		Arrays.fill(demoStateLogMarginal, Double.NEGATIVE_INFINITY);
		
		// add up the joint probabilities
		for (int epoch = 0; epoch < this.refinedDemography.numIntervals(); epoch++){
			for (int ancientDeme = 0; ancientDeme < this.refinedDemography.numAncientDemes(epoch); ancientDeme++){
					
				if (this.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
					assert (this.epochAncientDemeMarginalLogProbabilities.get(epoch)[ancientDeme] == Double.NEGATIVE_INFINITY);
					continue;
				}
				
				for (int presentDeme = 0; presentDeme < this.trunk.getSampleSizes().length; presentDeme++){
					
					int demoState = this.getDemoStateCollection().getDemoState(epoch, ancientDeme, presentDeme);
					
					// probability of whole deme
					double value = this.epochAncientDemeMarginalLogProbabilities.get(epoch)[ancientDeme];
					
					// multiply by number of lineages to get probability for refined-interval, present-deme
					value += Math.log(this.trunk.fractionAncientDemeToPresentDeme(epoch, presentDeme, ancientDeme));
					
					// add to probability of state-interval, present-deme
					demoStateLogMarginal[demoState] = LogSum.computePairLogSum(demoStateLogMarginal[demoState], value);
					
					assert(!Double.isNaN(demoStateLogMarginal[demoState]));

				}
			}
		}
		
		if (UberDemographyCore.RENORMALIZE_PROBS) renormalizeLogStochasticVector(demoStateLogMarginal, this.renormalizeEpsilon);
		
		return demoStateLogMarginal;

	}

	// the index to rCache is first <epoch> and then [startAncientDeme, firstEndAncientDeme, secondEndAncientDeme]
	private static <H extends FSAHaplotype> double W (int firstAbsorbInterval, int firstAbsorbAncientDeme, int secondAbsorbInterval, int secondAbsorbAncientDeme, ArrayList<ArrayList<double[][]>> Q, ArrayList<double[][]> rCache, Demography refinedDemography) {
		// some assertions
		// they can be everything
		assert ((0 <= firstAbsorbAncientDeme) && (firstAbsorbAncientDeme < refinedDemography.numAncientDemes(firstAbsorbInterval)));
		assert ((0 <= secondAbsorbAncientDeme) && (secondAbsorbAncientDeme < refinedDemography.numAncientDemes(secondAbsorbInterval)));

		// if they are in the wrong order, switch them
		if (secondAbsorbInterval < firstAbsorbInterval) {
			// call doubleU with changed roles
			return W (secondAbsorbInterval, secondAbsorbAncientDeme, firstAbsorbInterval, firstAbsorbAncientDeme, Q, rCache, refinedDemography);
		}
		
		// normal things
		int numAncientDemes = refinedDemography.numAncientDemes(firstAbsorbInterval);
		if (firstAbsorbInterval == secondAbsorbInterval) {
			// in this case we just have R
			return rCache.get(firstAbsorbInterval)[numAncientDemes + firstAbsorbAncientDeme][numAncientDemes + secondAbsorbAncientDeme];
		}
		else {
			assert (firstAbsorbInterval < refinedDemography.numIntervals()); 
			// we have firstAbsorbInterval < secondAbsorbInterval, so we need to multiply with Q
			double result = 0d;
			
			// loop over next start deme
			for (int nextStartAncientDeme=0; nextStartAncientDeme<refinedDemography.numAncientDemes(firstAbsorbInterval+1); nextStartAncientDeme++) {
				double tmp = 0d;
				
				// loop over content of next start deme
				for (int memberAncientDeme : refinedDemography.getMemberDemesIndices (nextStartAncientDeme, firstAbsorbInterval+1)) {
					// is the index proper?
					assert (memberAncientDeme <= refinedDemography.numAncientDemes(firstAbsorbInterval));
					// add some stuff
					tmp += rCache.get(firstAbsorbInterval)[numAncientDemes + firstAbsorbAncientDeme][memberAncientDeme];
				}
				
				// multiply it by Q
				tmp *= Q.get(firstAbsorbInterval+1).get(secondAbsorbInterval)[nextStartAncientDeme][secondAbsorbAncientDeme];
				// and add it
				result += tmp;
			}
			
			// give it away now
			assert (result > 0 - ONELOCUS_EPSILON);
			return result;
		}
	}
	
	// [epoch][ancient deme]
	// returns the joint probability of being absorbed into the given (epoch,ancient deme) and not recombining between the two loci
	protected ArrayList<double[]> computeAncientDemeNoRecoJointLog (ArrayList<double[]> jointP) {

		// WARING don't use the regular p, cause it's from the marginal, we need a no reco p

		// now calculate the actual y probabilities
		ArrayList<double[]> logJointY = new ArrayList<double[]>();

		// loop over all epochs
		for (int epoch = 0; epoch < this.refinedDemography.numIntervals(); epoch++){
			
			//make an array big enough
			double[] thisEpochY = new double[this.refinedDemography.numAncientDemes(epoch)];
			// add it to array list
			logJointY.add(thisEpochY);
		
			double[] thisJointP = jointP.get(epoch);
			// go through them ancient demes
			for (int ancientDeme=0; ancientDeme<thisEpochY.length; ancientDeme++ ){
		
				// scary ghost demes
				// this non-absorbing business should also take care of the pulse intervals
				if (this.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
					// log of zero
					thisEpochY[ancientDeme] = Double.NEGATIVE_INFINITY;
					// end next one
					continue;
				}
	
				// remember, it has to be absorbing, so + num demes
				int thisAbsorbDeme = this.refinedDemography.numAncientDemes(epoch) + ancientDeme;

				// logarithmify
				thisEpochY[ancientDeme] = Math.log (thisJointP[thisAbsorbDeme]);
				
				assert(!Double.isNaN(thisEpochY[ancientDeme]));
			}
		}
		
		// give it away now
		return logJointY;
	}

	// [prev epoch][next epoch][prev ancient deme][next ancient deme]
	// returns the joint probability of being absorbed into the given (prev epoch,prev ancient deme) at the first locus and recombining and eing absorbed into the given (next epoch,next ancient deme) at the second locus
	// the index to rCache is first <epoch> and then [startAncientDeme, firstEndAncientDeme, secondEndAncientDeme]
	protected ArrayList<ArrayList<double[][]>> computeAncientDemeRecoJoint (ArrayList<double[][]> rCache, ArrayList<double[]> jointP) {

		// WARING don't use the regular p. it's from the marginal, so not right for the joint
		// jointP should be the one you should be using
		// but the Q should be fine
		ArrayList<ArrayList<double[][]>> Q = demoTransitions.Q;

		
		// now iterate over all ancientDemes somehow
		// we have some precomputed p, Q, and R
		// so now just fill the stuff
		
		// storage
		ArrayList<ArrayList<double[][]>> joint = new ArrayList<ArrayList<double[][]>>();
		
		// go through previous epochs
		for (int prevAbsorbEpoch = 0; prevAbsorbEpoch < this.refinedDemography.numIntervals(); prevAbsorbEpoch++){
			
			// make some storage
			ArrayList<double[][]> prevEpochJoint = new ArrayList<double[][]>();
			joint.add(prevEpochJoint);
			
			// go through next epochs
			for (int nextAbsorbEpoch = 0; nextAbsorbEpoch < this.refinedDemography.numIntervals(); nextAbsorbEpoch++){
			
				// make some storage
				double[][] currEpochsJoint = new double[this.refinedDemography.numAncientDemes(prevAbsorbEpoch)][this.refinedDemography.numAncientDemes(nextAbsorbEpoch)];
				prevEpochJoint.add(currEpochsJoint);
				
				// go through them ancient demes
				for (int prevAbsorbDeme=0; prevAbsorbDeme<currEpochsJoint.length; prevAbsorbDeme++ ){				
					for (int nextAbsorbDeme=0; nextAbsorbDeme<currEpochsJoint[prevAbsorbDeme].length; nextAbsorbDeme++ ){				
		
						// scary ghost demes (all ghosts zero)
						// this should also take care of pulse intervals
						if (this.trunk.isNonAbsorbing(prevAbsorbEpoch, prevAbsorbDeme, UberDemographyCore.ONELOCUS_EPSILON) || this.trunk.isNonAbsorbing(nextAbsorbEpoch, nextAbsorbDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
							// log of zero
							currEpochsJoint[prevAbsorbDeme][nextAbsorbDeme] = Double.NEGATIVE_INFINITY;
							// end next one
							continue;
						}
		
						// which is the lowest epoch?
						int minAbsorbEpoch = Math.min(prevAbsorbEpoch, nextAbsorbEpoch);
						assert (minAbsorbEpoch <= this.refinedDemography.numIntervals());
						
						// init
						double preJoint = 0d;
						
						// and add the term for minAbsorbInterval
						double tmp = W (prevAbsorbEpoch, prevAbsorbDeme, nextAbsorbEpoch,nextAbsorbDeme, Q, rCache, this.refinedDemography);
						// remember tmp
						preJoint += tmp;
						
						// store the normalized log
						assert (preJoint > 0 - ONELOCUS_EPSILON);
						// have to be nice
						if (preJoint < 0) preJoint = 0.0;
						
						// take the log
						// and divide it by the absorbtion probs
						currEpochsJoint[prevAbsorbDeme][nextAbsorbDeme] = Math.log(preJoint);
						
						
						assert(!Double.isNaN(currEpochsJoint[prevAbsorbDeme][nextAbsorbDeme]));
						
					}
				}
			}
		}
		
		//give back the calculated transitions
		return joint;
	}

	
	// [epoch][ancient deme][trunk type][emission type]
	// This gives marginal probability of absorbing inside (epoch, ancient deme), times conditional probability of emission type given trunk type
	// this stuff will work best if the mutation matrix has 0 on the diagonal, since we will limit the number of mutation events
	protected ArrayList<double[][][]> computeAncientDemeEmissionLogProbabilities (double[] mutationRates, double[][] mutationMatrix) {
		
		assert mutationRates.length ==  this.refinedDemography.numIntervals();
		
		// WARNING: again we should use a different p
		// we make our own later
		// but we need the f
		ArrayList<double[][]> f = demoTransitions.f;

		// get some number
		int numAlleles = mutationMatrix.length;

		// compute mutation for each piece
		// [epoch,begin,nrMut,ending] dim=(numInt, numDeme, 2, 2*numDeme)
		ArrayList<double[][][]> intervalMutationProbs = new ArrayList<double[][][]>();
		// do some pre-computation for every interval
		for (int epoch = 0; epoch < this.refinedDemography.numIntervals(); epoch++) {
			// get the probabilities for having a mutation or not and doing certain transitions
			// should be of length 2, first thing is probs with no mutation, second is with one
			// the use of the absorption rates seems to be fine for now
			intervalMutationProbs.add(CoalescentODEs.computeMutationEvents (this.refinedDemography.epochList[epoch], mutationRates[epoch], refinedDemography.pulseMigrationMatrixList.get(epoch), refinedDemography.migrationMatrixList.get(epoch), this.trunk.getAbsorbRates(epoch), refinedDemography.expGrowthRates.get(epoch)));
		}

		// now make your own p, or something slightly different
		// we might have to do some more pre-processing here
		// transform the interval mutation probs into ending probs
		// [epoch,nrTotMut,ending] dim=(numInt, 2, 2*numDeme)
		// we have this.observedPresentDeme
		ArrayList<double[][]> mutationTransition = new ArrayList<double[][]>();
		
		// loop over all epochs
		for (int epoch = 0; epoch < this.refinedDemography.numIntervals(); epoch++){
			
			//make an array big enough
			// indexing should be [startDeme,EndDeme] dim=(numInt, 2*numInt)
			int thisNumDemes = refinedDemography.numAncientDemes(epoch);
			// get the previous one
			double[][] prevEnding = null;
			if (mutationTransition.size() > 0) {
				prevEnding = mutationTransition.get(mutationTransition.size()-1);
			}

			// we has zero or one mutation
			double[][] thisEpochP = new double[CoalescentODEs.MAX_NR_MUT_EVENTS+1][2*thisNumDemes];
			// just to be sure
			// add it to array list
			mutationTransition.add(thisEpochP);
			
			// get that one
			// [begin,nrMut,ending] dim=(numDeme, 2, 2*numDeme)
			double[][][] currTransition = intervalMutationProbs.get(epoch);
			double[][] thisF = f.get(epoch);

			// first one special
			if (epoch == 0) {
				assert (currTransition[observedPresentDeme].length == 2);
				// we could only have come from the oberseved deme
				for (int numMut=0; numMut<=CoalescentODEs.MAX_NR_MUT_EVENTS; numMut++) {
					assert (currTransition[observedPresentDeme][numMut].length == 2*thisNumDemes);
					for (int endDeme=0; endDeme<2*thisNumDemes; endDeme++) {
						thisEpochP[numMut][endDeme] = currTransition[observedPresentDeme][numMut][endDeme];
					}
				}
			}
			else {
				// get the right starting probabilities here
				// the new ones
				double[][] newStarts = new double[CoalescentODEs.MAX_NR_MUT_EVENTS+1][thisNumDemes];
				// now go through all them beginning demes
				// how many mutations?
				for (int numMut=0; numMut<=CoalescentODEs.MAX_NR_MUT_EVENTS; numMut++) {
					for (int startDeme=0; startDeme<thisNumDemes; startDeme++) {
						newStarts[numMut][startDeme] = 0d;
						// but actually we has to go through them members
						for (int memberAncientDeme : refinedDemography.getMemberDemesIndices (startDeme, epoch)) {
							// and put the right stuff here
							newStarts[numMut][startDeme] += prevEnding[numMut][memberAncientDeme];
						}
					}
				}
				
				// now compute the new probabilities
				// which end
				for (int endDeme=0; endDeme<2*thisNumDemes; endDeme++) {
					// how many mutations (in total, up to now)
					for (int numMut=0; numMut<=CoalescentODEs.MAX_NR_MUT_EVENTS; numMut++) {
						thisEpochP[numMut][endDeme] = 0d;
						// how could you have come here?
						// from which start
						for (int startDeme=0; startDeme<thisNumDemes; startDeme++) {
							assert (CoalescentODEs.MAX_NR_MUT_EVENTS == 1);
							// how many mutations along the way? (in total)
							if (numMut == 0) {
								// no mutation allowed
								thisEpochP[numMut][endDeme] += newStarts[0][startDeme] * currTransition[startDeme][0][endDeme];
							}
							else {
								// tertium non datur
								assert (numMut == 1);
								// either had a mutation already, and could not have one this time
								double prob = newStarts[1][startDeme] * thisF[startDeme][endDeme];
								// or we picked up a new one in this epoch
								prob += newStarts[0][startDeme] * currTransition[startDeme][1][endDeme];
								// save it
								thisEpochP[numMut][endDeme] += prob;
							}
						}
					}
				}
			}
			// this epoch should be good
		}
		// mutationP seems to be ok-ish now
		// now use it wisely		
		
		// now get some logarithmic emmisionProbabilties, somehow anyways
		// [epoch][ancient deme][trunk type][emission type]
		ArrayList<double[][][]> logJointEmission = new ArrayList<double[][][]>();
		
		// fill it
		// loop over epoch
		for (int epoch = 0; epoch < this.refinedDemography.numIntervals(); epoch++){
			// make storage and add it
			int numAncientDemes = this.refinedDemography.numAncientDemes(epoch);
			double[][][] thisEmission = new double[numAncientDemes][numAlleles][numAlleles];
			logJointEmission.add (thisEmission);
			
			// get mutation transition for this epoch
			double[][] thisMut = mutationTransition.get(epoch);
			
			
			// loop over ancient deme
			for (int ancientDeme=0; ancientDeme < thisEmission.length; ancientDeme++) {
				
				// loop over emission type
				for (int emissionType=0; emissionType<numAlleles; emissionType++) {
	
					// scary ghost demes (all ghosts zero, just to be sure)
					// this should deal with pulses too
					if (this.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
						// all the trunk types
						for (int trunkType=0; trunkType<numAlleles; trunkType++) {
							// log of zero
							thisEmission[ancientDeme][trunkType][emissionType] = Double.NEGATIVE_INFINITY;
						}
						// end next one
						continue;
					}
	

					// loop over trunk type
					for (int trunkType=0; trunkType<numAlleles; trunkType++) {

						double prob = 0d;
						// how many mutations along the way?
						for (int numMut=0; numMut<=CoalescentODEs.MAX_NR_MUT_EVENTS; numMut++) {
							assert (CoalescentODEs.MAX_NR_MUT_EVENTS == 1);
							if (numMut == 0) {
								// has to be the same or no chance
								if (trunkType == emissionType) {
									prob += thisMut[0][numAncientDemes + ancientDeme];
								}
							}
							else {
								// tertium non datur
								assert (numMut == 1);
								//can mutate
								prob += thisMut[1][numAncientDemes + ancientDeme] * mutationMatrix[trunkType][emissionType];
							}
						}
						
						assert (prob > 0 - ONELOCUS_EPSILON);
						// have to be nice
						if (prob < 0) prob = 0.0;
	
						// take the log
						thisEmission[ancientDeme][trunkType][emissionType] = Math.log(prob);

						
						assert(!Double.isNaN(thisEmission[ancientDeme][trunkType][emissionType]));
					}
				}
			}
		}
		
		// give it away now
		return logJointEmission;
		
	}	
	
	
	// [demoState]
	// returns conditional probability of no recombination given current demoState
	protected double[] computeDemoStateNoRecoTransition (ArrayList<double[]> epochAncientDemeNoRecoJoint) {
		
		double[] demoStateNoRecoTransition = new double[this.numDemoStates()];
		// fill with log 0s
		Arrays.fill(demoStateNoRecoTransition, Double.NEGATIVE_INFINITY);
		
		// add up the joint probabilities
		for (int epoch = 0; epoch < refinedDemography.numIntervals(); epoch++){
			for (int ancientDeme = 0; ancientDeme < refinedDemography.numAncientDemes(epoch); ancientDeme++){
				
				if (this.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
					assert (epochAncientDemeNoRecoJoint.get(epoch)[ancientDeme] == Double.NEGATIVE_INFINITY);
					continue;
				}
				
				for (int presentDeme = 0; presentDeme < this.trunk.getSampleSizes().length; presentDeme++) {
					
					int demoState = this.getDemoStateCollection().getDemoState(epoch, ancientDeme, presentDeme);
					
					// joint probability of event
					double value = epochAncientDemeNoRecoJoint.get(epoch)[ancientDeme];
					
					// multiply by number of lineages to get probability for the refined-interval, present-deme
					value += Math.log(this.trunk.fractionAncientDemeToPresentDeme(epoch, presentDeme, ancientDeme));
					
					// add to probability of state-interval, present-deme
					demoStateNoRecoTransition[demoState] = LogSum.computePairLogSum(demoStateNoRecoTransition[demoState], value);
				}
			}
		}
		
		// divide by the marginal to get conditional
		for (int demoState = 0; demoState < this.numDemoStates(); demoState++){
			int currDeme = this.getDemoStateCollection().getPresentDeme(demoState);
			assert(this.trunk.getSampleSizes()[currDeme] != 0 || (Double.NEGATIVE_INFINITY == this.demoStateLogMarginal[demoState]));
			assert(this.trunk.getSampleSizes()[currDeme] != 0 || (Double.NEGATIVE_INFINITY == demoStateNoRecoTransition[demoState]));

			if (Double.NEGATIVE_INFINITY == (this.demoStateLogMarginal[demoState])) {
				
				assert(Double.NEGATIVE_INFINITY == (demoStateNoRecoTransition[demoState]));
				continue;
			}
			
			demoStateNoRecoTransition[demoState] -= this.demoStateLogMarginal[demoState];
			
			assert(!Double.isNaN(demoStateNoRecoTransition[demoState]));
		}
		
		return demoStateNoRecoTransition;
	}

	// [prevDemoState][nextDemoState]
	// returns conditional probability of recombination to next demoState given previous demoState
	protected double[][] computeDemoStateRecoTransition (ArrayList<ArrayList<double[][]>> epochAncientDemeRecoJoint) {
		
		double[][] demoStateRecoTransition = new double[this.numDemoStates()][this.numDemoStates()];
		// fill with log 0s
		for (int i=0; i<demoStateRecoTransition.length; i++) {
			Arrays.fill(demoStateRecoTransition[i], Double.NEGATIVE_INFINITY);
		}
		
		// add up the joint probabilities for a single lineage residing in each interval/deme
		// loop through previuos refined intervals, ancient demes
		for (int prevEpoch = 0; prevEpoch < this.refinedDemography.numIntervals(); prevEpoch++){
			for (int prevAncientDeme = 0; prevAncientDeme < this.refinedDemography.numAncientDemes(prevEpoch); prevAncientDeme++){
				// loop through next refined intervals, ancient demes
				for (int nextEpoch = 0; nextEpoch < this.refinedDemography.numIntervals(); nextEpoch++){
					for (int nextAncientDeme = 0; nextAncientDeme < this.refinedDemography.numAncientDemes(nextEpoch); nextAncientDeme++){
					
						// skip ghost demes
						if (this.trunk.isNonAbsorbing(prevEpoch, prevAncientDeme, UberDemographyCore.ONELOCUS_EPSILON)|| this.trunk.isNonAbsorbing(nextEpoch, nextAncientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
							assert (epochAncientDemeRecoJoint.get(prevEpoch).get(nextEpoch)[prevAncientDeme][nextAncientDeme] == Double.NEGATIVE_INFINITY);
							continue;
						}
						
						// loop through the present demes at both loci
						for (int prevPresentDeme = 0; prevPresentDeme < this.trunk.getSampleSizes().length; prevPresentDeme++){
							for (int nextPresentDeme = 0; nextPresentDeme < this.trunk.getSampleSizes().length; nextPresentDeme++){
						
								// get the right intervalDeme index
								int prevDemoState = this.getDemoStateCollection().getDemoState(prevEpoch, prevAncientDeme, prevPresentDeme);
								int nextDemoState = this.getDemoStateCollection().getDemoState(nextEpoch, nextAncientDeme, nextPresentDeme);
								
								// joint probability of event
								double value = epochAncientDemeRecoJoint.get(prevEpoch).get(nextEpoch)[prevAncientDeme][nextAncientDeme];
								
								// multiply by present deme probabilities
								value += Math.log(this.trunk.fractionAncientDemeToPresentDeme(prevEpoch, prevPresentDeme, prevAncientDeme));
								value += Math.log(this.trunk.fractionAncientDemeToPresentDeme(nextEpoch, nextPresentDeme, nextAncientDeme));

								
								// add to probability of state-interval, present-deme
								demoStateRecoTransition[prevDemoState][nextDemoState] = LogSum.computePairLogSum(demoStateRecoTransition[prevDemoState][nextDemoState], value);
							}
						}
					}
				}
			}
		}
		
		// divide by the marginal at prev locus to get conditional
		for (int prevDemoState = 0; prevDemoState < this.numDemoStates(); prevDemoState++){
			int prevDeme = this.getDemoStateCollection().getPresentDeme(prevDemoState);

			
			for (int nextDemoState = 0; nextDemoState < this.numDemoStates(); nextDemoState++){

				assert(this.trunk.getSampleSizes()[prevDeme] != 0 || (Double.NEGATIVE_INFINITY == this.demoStateLogMarginal[prevDemoState]));
				assert(this.trunk.getSampleSizes()[prevDeme] != 0 || (Double.NEGATIVE_INFINITY == demoStateRecoTransition[prevDemoState][nextDemoState]));
				
				if (Double.NEGATIVE_INFINITY == this.demoStateLogMarginal[prevDemoState]) {
					
					assert(Double.NEGATIVE_INFINITY == demoStateRecoTransition[prevDemoState][nextDemoState]);
					continue;
				}
				
				demoStateRecoTransition[prevDemoState][nextDemoState] -= this.demoStateLogMarginal[prevDemoState];
				
				assert(!Double.isNaN(demoStateRecoTransition[prevDemoState][nextDemoState]));
			}
		}
		
		// should be fine and return it
		return demoStateRecoTransition;
	}

	
	// first joint, then logNoReco
	private static SimplePair<ArrayList<double[]>, ArrayList<double[]>> computeJointP (int observedPresentDeme, Demography refinedDemography, ArrayList<double[]> noRecoCache) {

		// also the jointP
		ArrayList<double[]> jointP = new ArrayList<double[]>();
		// now calculate the actual y probabilities
		ArrayList<double[]> noReco = new ArrayList<double[]>();
		
		double[] endingProbs = null;
		
		// loop over all epochs
		for (int epoch = 0; epoch < refinedDemography.numIntervals(); epoch++){
			
			//make an array big enough
			// indexing should be [startDeme,EndDeme] dim=(numInt, 2*numInt)
			int thisNumDemes = refinedDemography.numAncientDemes(epoch);

			double[] thisEpochP = null;
			
			double[] thisEpochNoReco = noRecoCache.get(epoch);
			assert (thisEpochNoReco.length == 2*thisNumDemes);
			
			// first one special
			if (epoch == 0) {
				// only one guy for p
				thisEpochP = new double[thisNumDemes];
				Arrays.fill(thisEpochP, 0d);
				thisEpochP[observedPresentDeme] = 1d;
				
				// compute the probs for noReco
				endingProbs = new double[2*thisNumDemes];
				for (int i=0; i<endingProbs.length; i++) {
					endingProbs[i] = thisEpochNoReco[i];
				}
			}
			else {
				
				// get some strating probs (which by chance happens to be what we want to put in P)
				thisEpochP = new double[thisNumDemes];
				Arrays.fill(thisEpochP, 0d);
				// go through them this start demes
				for (int thisStartDeme=0; thisStartDeme<thisEpochP.length; thisStartDeme++ ){
					// and go through the members in the previous deme
					for (int memberAncientDeme : refinedDemography.getMemberDemesIndices (thisStartDeme, epoch)) {
						assert (memberAncientDeme < refinedDemography.numAncientDemes(epoch-1));
						assert (2*memberAncientDeme < endingProbs.length);
						thisEpochP[thisStartDeme] += endingProbs[memberAncientDeme];
					}
				}
				
				
				endingProbs = new double[2*thisNumDemes];
				Arrays.fill (endingProbs, 0d);
				// and now some  transition to get to the end
				for (int endDeme=0; endDeme<2*thisNumDemes; endDeme++) {
					// from where did we come from?
					endingProbs[endDeme] = thisEpochNoReco[endDeme];
				}
				// this should be it
			}
			
			// add it to array list
			jointP.add(thisEpochP);
			
			// add something to endings too
			double[] realEndings = new double[thisNumDemes];
			for (int i=0; i<realEndings.length; i++) {
				realEndings[i] = Math.log(endingProbs[thisNumDemes + i]);
			}
			noReco.add(realEndings);
		}
		
		// seems to be ok-ish
		// give it away now
		return new SimplePair<ArrayList<double[]>, ArrayList<double[]>> (jointP, noReco);
	}

	// returns conditional probability of recombination to next demoState given previous demoState
	@Override
	protected RecoLogProbs computeRecoLogProbs (double recoRate) {
		
		// precompute some cashes
		SimplePair<ArrayList<double[][]>, ArrayList<double[]>> cachePair = preComputeRecoCaches (recoRate, observedPresentDeme, this.trunk);

		// precompute some caches
		ArrayList<double[][]> rCache = cachePair.first();
		// indexing noReco via [startDeme,endDeme] dim (numDemes, 2*numDemes)
		ArrayList<double[]> noRecoCache = cachePair.second();

		// now we have pre computed some stuff, so other stuff
		SimplePair<ArrayList<double[]>, ArrayList<double[]>> pair = computeJointP (observedPresentDeme, refinedDemography, noRecoCache);
		ArrayList<double[]> jointP = pair.first();

		// first the no recombination part
		ArrayList<double[]> ancientNoRecoLogJoint = pair.second();
				
		double[] noReco = computeDemoStateNoRecoTransition (ancientNoRecoLogJoint); 
		
		// then the no reco part
		ArrayList<ArrayList<double[][]>> ancientRecoJoint = computeAncientDemeRecoJoint (rCache, jointP);
		
		double[][] reco = computeDemoStateRecoTransition (ancientRecoJoint);
		
		// put together the return object
		RecoLogProbs demoStateTransition = new RecoLogProbs (reco, noReco);
		if (UberDemographyCore.RENORMALIZE_PROBS) ODECore.renormalizeLogTransitions (demoStateTransition, this.demoStateLogMarginal, this.renormalizeEpsilon);
				
		// should be fine and return it
		return demoStateTransition;
	}

	private SimplePair<ArrayList<double[][]>, ArrayList<double[]>> preComputeRecoCaches (double recoRate, int observedPresentDeme, TrunkProcess trunk) {
		ArrayList<double[][]> rCache = new ArrayList<double[][]>();
		// indexing noReco via [endDeme] dim (2*numDemes)
		ArrayList<double[]> noRecoCache = new ArrayList<double[]>();
		// let's fill it, do some ODEs in each epoch
		SimplePair<double[], double[][]> thisPair = null;
		
		for (int epoch = 0; epoch < trunk.refinedDemography.numIntervals(); epoch++) {
			
			// precompute R
			thisPair = CoalescentODEs.computeR(trunk, epoch, observedPresentDeme, thisPair, recoRate, smcPrime);

			// add things to the containers
			rCache.add (thisPair.second());
			noRecoCache.add (thisPair.first());
		}
		
		return new SimplePair<ArrayList<double[][]>, ArrayList<double[]>> (rCache, noRecoCache);
	}

	// [demoState][trunkType][emissionType]
	// returns conditional probability of emission type, given absorption into a lineage with trunk type
	@Override
	protected double[][][] computeLogEmission (double[] mutationRate, MyEigenDecomposition eigenMutation) {
		// [epoch][ancient deme][trunk type][emission type]
		// This gives marginal probability of absorbing inside (epoch, ancient deme), times conditional probability of emission type given trunk type
		// we only use the eigenstuff to get the original matrix out
		double[] refinedMutationRate = this.getRefinedTheta(mutationRate);
		
		ArrayList<double[][][]> refinedIntervalAncientDemeEmissionLogProbs = computeAncientDemeEmissionLogProbabilities (refinedMutationRate, eigenMutation.getOriginalMatrix());
		
		int numAlleles = eigenMutation.lambda.length;
		
		// to be returned
		double[][][] demoStateEmission = new double[this.numDemoStates()][numAlleles][numAlleles];
		for (int i = 0; i < demoStateEmission.length; i++){
			for (int j=0; j < demoStateEmission[i].length; j++){
				for (int k = 0; k < demoStateEmission[i][j].length; k++) {
					demoStateEmission[i][j][k] = Double.NEGATIVE_INFINITY;
				}
			}
		}
		
		// add up the joint probabilities for a single lineage residing in each interval/deme
		for (int epoch = 0; epoch < this.refinedDemography.numIntervals(); epoch++){
			for (int ancientDeme = 0; ancientDeme < this.refinedDemography.numAncientDemes(epoch); ancientDeme++){
				for (int presentDeme = 0; presentDeme < this.trunk.getSampleSizes().length; presentDeme++){
					
					for (int emissionType = 0; emissionType < numAlleles; emissionType++){
						for (int trunkType = 0; trunkType < numAlleles; trunkType++){
							
							if (this.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
								assert (refinedIntervalAncientDemeEmissionLogProbs.get(epoch)[ancientDeme][trunkType][emissionType] == Double.NEGATIVE_INFINITY);
								continue;
							}
							
							int demoState = this.getDemoStateCollection().getDemoState(epoch, ancientDeme, presentDeme);
							
							// probability absorption times conditional probability of emission, for refined-interval/ancient-deme
							double value = refinedIntervalAncientDemeEmissionLogProbs.get(epoch)[ancientDeme][trunkType][emissionType];
							
							// convert it to the present deme
							value += Math.log(this.trunk.fractionAncientDemeToPresentDeme(epoch, presentDeme, ancientDeme));
							
							// add it to the emission
							demoStateEmission[demoState][trunkType][emissionType] = LogSum.computePairLogSum(demoStateEmission[demoState][trunkType][emissionType], value);
						}
					}
					
				}			
			}
		}
		
		// divide by the marginal to get conditional
		for (int demoState = 0; demoState < this.numDemoStates(); demoState++){
			int currDeme = this.getDemoStateCollection().getPresentDeme(demoState);
			
			for (int emissionType = 0; emissionType < numAlleles; emissionType++){
				for (int trunkType = 0; trunkType < numAlleles; trunkType++){
					
					
					// being a ghost deme should imply some zero (joint) probabilities
					assert(this.trunk.getSampleSizes()[currDeme] != 0 || (Double.NEGATIVE_INFINITY == this.demoStateLogMarginal[demoState]));
					assert(this.trunk.getSampleSizes()[currDeme] != 0 || (Double.NEGATIVE_INFINITY == demoStateEmission[demoState][trunkType][emissionType]));

					if (Double.NEGATIVE_INFINITY == (this.demoStateLogMarginal[demoState])) {
						
						assert(Double.NEGATIVE_INFINITY == demoStateEmission[demoState][trunkType][emissionType]);
						continue;
					}
					
					// divide by the marginal probability of this interval deme
					demoStateEmission[demoState][trunkType][emissionType] -= (this.demoStateLogMarginal[demoState]);

					assert(!Double.isNaN(demoStateEmission[demoState][trunkType][emissionType]));
				}
			}
			
			// renomalize it
			if (this.demoStateLogMarginal[demoState] != Double.NEGATIVE_INFINITY) {
				if (UberDemographyCore.RENORMALIZE_PROBS) renormalizeLogStochasticMatrix(demoStateEmission[demoState], this.renormalizeEpsilon);
			}
		}
		
		return demoStateEmission;		
	}
	
	@Override
	public double[] getDemoStateLogMarginal() {
		return demoStateLogMarginal;
	}
	
	@Override
	public Interval getTimeInterval(int demoStateIdx) {
		return new Interval (this.getDemoStateCollection().startTime(demoStateIdx) , this.getDemoStateCollection().endTime(demoStateIdx));
	}

	// note: this function may need to be called before demoStateLogMarginal has been constructed
	@Override
	public int numDemoStates () {
		return this.getDemoStateCollection().numDemoStates();
	}
	
	@Override
	public void checkAncientEmissions (double[] mutationRate, MyEigenDecomposition eigenMutation) {
		ArrayList<double[][][]> ancientDemeEmissions = computeAncientDemeEmissionLogProbabilities (mutationRate, eigenMutation.getOriginalMatrix());
		
		int numAlleles = eigenMutation.lambda.length;
		
		for (int epoch=0; epoch<this.refinedDemography.numIntervals(); epoch++) {
			for (int ancientDeme=0; ancientDeme<this.refinedDemography.numAncientDemes(epoch); ancientDeme++) {
				
				for (int trunkType = 0; trunkType < numAlleles; trunkType++) {
					
					double sum = 0d;
					for (int emissionType = 0; emissionType < numAlleles; emissionType++) {
						sum += Math.exp(ancientDemeEmissions.get(epoch)[ancientDeme][trunkType][emissionType]);
					}
					
					sum /= Math.exp(this.epochAncientDemeMarginalLogProbabilities.get(epoch)[ancientDeme]);
					System.out.println("Ancient emission sum: " + sum);
				}
				
				
			}
		}
		
	}
	
	
	
	@Override
	public void checkAncientDemeProbs(double recoRate) {
		

		// precompute some cashes
		SimplePair<ArrayList<double[][]>, ArrayList<double[]>> cachePair = preComputeRecoCaches (recoRate, observedPresentDeme, this.trunk);
		ArrayList<double[][]> rCache = cachePair.first();
		// indexing noReco via [startDeme,endDeme] dim (numDemes, 2*numDemes)
		ArrayList<double[]> noRecoCache = cachePair.second();

		
		// now we have pre computed some stuff, so other stuff
		SimplePair<ArrayList<double[]>, ArrayList<double[]>> pair = computeJointP (observedPresentDeme, refinedDemography, noRecoCache);
		ArrayList<double[]> jointP = pair.first();

		// first the no recombination part
		ArrayList<double[]> ancientNoRecoJoint = pair.second();

		
		ArrayList<ArrayList<double[][]>> epochAncientDemeReco = computeAncientDemeRecoJoint (rCache, jointP);

		double[][] reco = computeDemoStateRecoTransition(epochAncientDemeReco);
		double[] noReco = computeDemoStateNoRecoTransition(ancientNoRecoJoint);
		
		// check no reco
		
		double noRecoSumOld = 0d;
		double noRecoSumNew = 0d;
		
		for (int epoch=0; epoch<this.refinedDemography.numIntervals(); epoch++) {
			for (int ancientDeme=0; ancientDeme<this.refinedDemography.numAncientDemes(epoch); ancientDeme++) {
				noRecoSumOld += Math.exp(ancientNoRecoJoint.get(epoch)[ancientDeme]);
			}
		}
		for (int demoState=0; demoState<this.numDemoStates(); demoState++) {
			double marginal = Math.exp (this.demoStateLogMarginal[demoState]);
			noRecoSumNew += Math.exp (noReco[demoState]) * marginal;
		}
		System.out.println ("old Y: " + noRecoSumOld);
		System.out.println ("sum Y diff: " + (noRecoSumOld - noRecoSumNew));
		
		// check reco
		
		double recoSumOld = 0d;
		double recoSumNew = 0d;
		
		for (int prevEpoch=0; prevEpoch<this.refinedDemography.numIntervals(); prevEpoch++) {
			for (int prevAncientDeme=0; prevAncientDeme<this.refinedDemography.numAncientDemes(prevEpoch); prevAncientDeme++) {
				for (int nextEpoch=0; nextEpoch<this.refinedDemography.numIntervals(); nextEpoch++) {
					for (int nextAncientDeme=0; nextAncientDeme<this.refinedDemography.numAncientDemes(nextEpoch); nextAncientDeme++) {
						recoSumOld += Math.exp(epochAncientDemeReco.get(prevEpoch).get(nextEpoch)[prevAncientDeme][nextAncientDeme]) ;

					}
				}
			}
		}
		for (int prevDemoState=0; prevDemoState<this.numDemoStates(); prevDemoState++) {
			for (int nextDemoState=0; nextDemoState<this.numDemoStates(); nextDemoState++) {
				double marginal = Math.exp (this.demoStateLogMarginal[prevDemoState]);
				recoSumNew += Math.exp (reco[prevDemoState][nextDemoState]) * marginal;
			}
		}
		System.out.println ("old Z: " + recoSumOld);
		System.out.println ("sum Z diff: " + (recoSumOld - recoSumNew));
		
		System.out.println ("sum old Y Z: " + (recoSumOld + noRecoSumOld));
		
		// check it
		for (int prevEpoch = 0; prevEpoch < this.refinedDemography.numIntervals(); prevEpoch++){
			for (int prevAncientDeme = 0; prevAncientDeme < this.refinedDemography.numAncientDemes(prevEpoch); prevAncientDeme++){
				// this has to sum to one				
				double transitionSum = Math.exp(ancientNoRecoJoint.get(prevEpoch)[prevAncientDeme]);
				
				for (int nextEpoch = 0; nextEpoch < this.refinedDemography.numIntervals(); nextEpoch++){
					for (int nextAncientDeme = 0; nextAncientDeme < this.refinedDemography.numAncientDemes(nextEpoch); nextAncientDeme++){
						double value = epochAncientDemeReco.get(prevEpoch).get(nextEpoch)[prevAncientDeme][nextAncientDeme];
						transitionSum += Math.exp(value);
					}
				}
				
				transitionSum /= Math.exp (epochAncientDemeMarginalLogProbabilities.get(prevEpoch)[prevAncientDeme]);
				
				System.out.println("transition: " + transitionSum);
			}
		}
		
	}
	
	@Override
	public void checkAncientMarginal() {
		double sum = 0d;
		for (double[] theseDemes : this.epochAncientDemeMarginalLogProbabilities) {
			for (int i=0; i<theseDemes.length; i++) {
				sum += Math.exp(theseDemes[i]);
			}
		}
		System.out.println ("marginal: " + sum);
	}


	
	@Override
	public void dump (double recoRate, double[] mutRate, double[][] mutMatrix, PrintStream outStream) {
		// show things
		
		// marginal
		outStream.println("=== Marginal probabilities:");
		double logSum = Double.NEGATIVE_INFINITY;
		for (double value : this.demoStateLogMarginal) {
			outStream.print (Math.exp (value) + "\t");
			logSum = LogSum.computePairLogSum(logSum, value);
		}
		outStream.println();
		outStream.println("=== logSum: " + logSum);
		
		// no reco & reco
		RecoLogProbs pair = this.getRecoLogProbs(recoRate);
		double[] logNoReco = pair.logNoReco;
		double[][] logReco = pair.logReco; 
		double[] logSums = new double[logNoReco.length];
		for (int i=0; i<logSums.length; i++) {
			logSums[i] = Double.NEGATIVE_INFINITY;
		}
		
		// no reco
		outStream.println();
		outStream.println("=== no Reco:");
		for (int i=0; i<logNoReco.length; i++) {
			outStream.print (Math.exp(logNoReco[i]) + "\t");
			logSums[i] = LogSum.computePairLogSum(logSums[i], logNoReco[i]);
		}
		outStream.println();
		
		// reco
		outStream.println();
		outStream.println("=== Reco (i times j):");
		for (int i=0; i<logReco.length; i++) {
			for (int j=0; j<logReco[i].length; j++) {
				outStream.print (Math.exp(logReco[i][j]) + "\t");
				logSums[i] = LogSum.computePairLogSum(logSums[i], logReco[i][j]);
			}
			if (i < logReco.length-1) outStream.println();
		}
		outStream.println();
		outStream.println("=== logSums: " + Arrays.toString(logSums));
		
		// some emission stuff
		double[][][] logEmission = this.getLogEmission (mutRate, new MyEigenDecomposition(mutMatrix, UberDemographyCore.ONELOCUS_EPSILON));
		outStream.println();
		outStream.println("=== Emission:");
		for (int demoState=0; demoState < logEmission.length; demoState++) {
			outStream.println();
			outStream.println("[time " + demoState + "]");
			logSums = new double[logEmission[demoState].length];
			for (int i=0; i<logSums.length; i++) {
				logSums[i] = Double.NEGATIVE_INFINITY;
			}
			outStream.println("=== Mutation (trunk times observed):");
			for (int i=0; i<logEmission[demoState].length; i++) {
				for (int j=0; j<logEmission[demoState][i].length; j++) {
					outStream.print (Math.exp(logEmission[demoState][i][j]) + "\t");
					logSums[i] = LogSum.computePairLogSum(logSums[i], logEmission[demoState][i][j]);
				}
				if (i < logEmission[demoState].length-1) outStream.println();
			}
			outStream.println();
			outStream.println("=== logSums: " + Arrays.toString(logSums));
		}
		

	}
	
	///  pair for p and Q
	private static class DemographyTransitions {
		// constructor
		public DemographyTransitions (ArrayList<ArrayList<double[][]>> p, ArrayList<ArrayList<double[][]>> Q, ArrayList<double[][]> f) {
			this.p = p;
			this.Q = Q;
			this.f = f;
		}
		
		// members
		// [first interval][second interval][first start deme][second start deme]
		final ArrayList<ArrayList<double[][]>> p;
		// [first interval][absorbing interval][first start deme][absorbing deme]
		final ArrayList<ArrayList<double[][]>> Q;
		// [interval][startDeme][endDeme]
		final ArrayList<double[][]> f;
	}
	
	// getter
	@Override
	public boolean isConditionalConfigEmpty () {
		return conditionalConfigEmpty;
	}
	
	@Override
	public int getObservedPresentDeme () {
		return observedPresentDeme;
	}
	
	@Override
	public TrunkProcess getTrunk() {
		return trunk;
	}

	private final boolean smcPrime;
	
	protected final TrunkProcess trunk;
	final private int observedPresentDeme;
	public final Demography refinedDemography;
	
	// to deal with fully or partially empty configurations
	private final boolean conditionalConfigEmpty;
	
	// some precomputed stuff
	protected final DemographyTransitions demoTransitions; 
	
	// initial [refined interval][ancient deme] probabilities, for a single lineage
	protected final ArrayList<double[]> epochAncientDemeMarginalLogProbabilities;
	
	protected final double[] demoStateLogMarginal;

}

