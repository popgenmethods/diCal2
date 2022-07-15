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
import java.util.HashMap;
import java.util.List;

import Jama.Matrix;
import Jampack.JampackException;
import Jampack.Z;
import edu.berkeley.diCal2.csd.DemoState.DemoStateCollection;
import edu.berkeley.diCal2.csd.TrunkProcess.SimpleTrunk;
import edu.berkeley.diCal2.csd.auxiliary.ComplexMath;
import edu.berkeley.diCal2.csd.auxiliary.MyEigenDecomposition;
import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.haplotype.FSAHaplotype;
import edu.berkeley.diCal2.utility.LogSum;
import edu.berkeley.diCal2.utility.RateMatrixTools;
import edu.berkeley.diCal2.utility.RealPartition.Interval;
import edu.berkeley.diCal2.utility.SumArray;

public class EigenCore extends UberDemographyCore {

	public EigenCore (int[] sampleSizes, int observedPresentDeme, Demography demo, Interval[] hiddenStateIntervals) {
		this (sampleSizes, new DemoStateCollection(demo, hiddenStateIntervals, null, UberDemographyCore.ONELOCUS_EPSILON), observedPresentDeme);
	}
	
	public EigenCore (int[] sampleSizes, DemoStateCollection demoStates, int observedPresentDeme) {
		this (new SimpleTrunk(sampleSizes, demoStates.refinedDemography), demoStates, observedPresentDeme, false);
	}
	
	public EigenCore (TrunkProcess trunk, DemoStateCollection demoStates, int observedPresentDeme, boolean cacheRecoTransitions) {
		this (trunk, demoStates, observedPresentDeme, cacheRecoTransitions, UberDemographyCore.DEFAULT_RENORMALIZATION_EPSILON);
	}
	
	/// the constructor, gets configuration, parameters and number of intervals
	public EigenCore (TrunkProcess trunk, DemoStateCollection demoStates, int observedPresentDeme, boolean cacheRecoTransitions, double renormalizeEpsilon) {
		super(demoStates, cacheRecoTransitions, renormalizeEpsilon);
		
		// remember some things
		this.trunk = trunk;
		this.observedPresentDeme = observedPresentDeme;

		// is it empty?
		this.conditionalConfigEmpty = trunk == null ? true : (SumArray.getSum(trunk.getSampleSizes()) == 0);
		
		// again, if the conditional config is empty, we don't have to do anything
		if (!this.conditionalConfigEmpty) {

			this.refinedDemography = this.trunk.getRefinedDemography();
			assert (this.refinedDemography == this.getDemoStateCollection().refinedDemography);
			
						// construct the migration parameters [give cake stuff too] {refined guy}
			this.migrationParams = new StructuredMigrationParameters (this.trunk);
	
			// locus independent part of the HMM
			// precompute the non-absoprtion and absorption probabilities
			this.demoTransitions = preComputeDemoTransitions (this.observedPresentDeme, this.migrationParams, this.renormalizeEpsilon);

			// initialize the marginal intervalDeme probabilities
			this.epochAncientDemeMarginalLogProbabilities = computeEpochAncientDemeLogMarginal (this.observedPresentDeme, this.migrationParams, this.demoTransitions.Q);
			
			for (double[] thisEpochMarginalProbs : epochAncientDemeMarginalLogProbabilities){
				for (double value : thisEpochMarginalProbs){
					assert (Math.exp(value) > 0d - ONELOCUS_EPSILON);
				}
			}
			
			// initialize the marginal intervalDeme probabilities
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
			this.migrationParams = null;
			this.refinedDemography = null;
		}	
	}

	private static <H extends FSAHaplotype> DemographyTransitions preComputeDemoTransitions (int observedPresentDeme, StructuredMigrationParameters migrationParams, double renormalizeEpsilon) {
 
		// for each interval we have a matrix (includes the absorbing states)
		// [interval][starting state][ending state]
		ArrayList<double[][]> f = new ArrayList<double[][]>();
		
		// go through intervals
		for (int e=0; e<migrationParams.numIntervals(); e++) {
			int numAncientDemes = migrationParams.numAncientDemes(e);
			double[][] currF;
			// is pulse?
			if (migrationParams.trunk.getRefinedDemography().isPulse(e)) {
				// its a pulse!
				// initialize all to zero
				currF = new double[2*numAncientDemes][2*numAncientDemes];
				// just copy the pulse matrix for the non-absorbing demes
				double[][] pulseMatrix = migrationParams.trunk.getRefinedDemography().pulseMigrationMatrixList.get(e);
				for (int i=0; i<numAncientDemes; i++) {
					for (int j=0; j<numAncientDemes; j++) {
						currF[i][j] = pulseMatrix[i][j];
					}
					// and fill in some 1s on the diagonal to make it a stochastic matrix
					currF[numAncientDemes+i][numAncientDemes+i] = 1d;
				}
			}
			else {
				// non-pulse, so do the usual eigenstuff
				// get the right eigenstuff
				MyEigenDecomposition eigenStuff = migrationParams.getEigenStuffList().get(e);

				
				// now fill the matrix of non-absorption probabilities
				assert (eigenStuff.V.nc == 2*numAncientDemes);
				assert (eigenStuff.V.nr == 2*numAncientDemes);
				assert (eigenStuff.lambda.length == 2*numAncientDemes);
				
				currF = eigenStuff.getRateMatrixExponential(migrationParams.getIntervalList()[e].endPoint - migrationParams.getIntervalList()[e].startPoint);
			}
			// matrix should sum to one
			
			// renormalize the f's so every row sums to one (they should do so anyways, but numerical stability ...)
			if (UberDemographyCore.RENORMALIZE_PROBS) UberDemographyCore.renormalizeStochasticMatrix(currF, renormalizeEpsilon);
			assert (RateMatrixTools.isProperStochaticMatrix (currF, UberDemographyCore.ONELOCUS_EPSILON));
			
			f.add (currF);
		}
		// so these matrices should be fine now
		
		
		int numIntervals = migrationParams.numIntervals();
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
					// initialize the P
					currP.add(Matrix.identity(migrationParams.numAncientDemes(firstInterval), migrationParams.numAncientDemes(firstInterval)).getArray());
				}
				else {
					// step forwward to next matrix
					// last in list should be the previous matrix
					double[][] prevMatrix = currP.get(currP.size()-1);
					assert (prevMatrix.length == migrationParams.numAncientDemes(firstInterval));
					assert (prevMatrix[0].length == migrationParams.numAncientDemes(secondInterval-1));
					// now step forward
					double[][] newMatrix = new double[migrationParams.numAncientDemes(firstInterval)][migrationParams.numAncientDemes(secondInterval)];
					// fill it
					// loop over first starting deme
					for (int firstAncientDeme=0; firstAncientDeme<newMatrix.length; firstAncientDeme++) {
						// loop over second starting deme
						for (int secondAncientDeme=0; secondAncientDeme<newMatrix[0].length; secondAncientDeme++) {
							// initialize
							newMatrix[firstAncientDeme][secondAncientDeme] = 0d;
							for (int memberAncientDeme : migrationParams.getDemography().getMemberDemesIndices (secondAncientDeme, secondInterval)) {
								assert (memberAncientDeme <= migrationParams.numAncientDemes(secondInterval-1));
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
					double[][] newMatrix = new double[migrationParams.numAncientDemes(firstInterval)][migrationParams.numAncientDemes(absorbingInterval)];
					// fill it
					// loop over first starting deme
					for (int startAncientDeme=0; startAncientDeme<newMatrix.length; startAncientDeme++) {
						// loop over second starting deme
						for (int absorbingAncientDeme=0; absorbingAncientDeme<newMatrix[0].length; absorbingAncientDeme++) {
							// initialize
							newMatrix[startAncientDeme][absorbingAncientDeme] = 0d;
							// starting deme in the absorbing interval
							for (int intermediateAncientDeme=0; intermediateAncientDeme<migrationParams.numAncientDemes(absorbingInterval); intermediateAncientDeme++) {
								// we want to go from startDeme to intermediate deme to absorbing deme
								// remeber the indexnig for the absorbing states in the extended matrices
								newMatrix[startAncientDeme][absorbingAncientDeme] += p.get(firstInterval).get(absorbingInterval)[startAncientDeme][intermediateAncientDeme] * f.get(absorbingInterval)[intermediateAncientDeme][migrationParams.numAncientDemes(absorbingInterval) + absorbingAncientDeme];
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
	private static <H extends FSAHaplotype> ArrayList<double[]> computeEpochAncientDemeLogMarginal (int observedPresentDeme, StructuredMigrationParameters migrationParams, ArrayList<ArrayList<double[][]>> q ) {

		// fill the initial probabilities from q
		// it's just the aborption stuff when starting in the first interval
		ArrayList<double[]> initialProb = new ArrayList<double[]>();
		
		// go through epochs
		for (int epoch = 0; epoch < migrationParams.numIntervals(); epoch++){
			
			//make an array big enough
			double[] thisEpochLogMarginal = new double[migrationParams.numAncientDemes(epoch)];
			// add it to array list
			initialProb.add(thisEpochLogMarginal);
		
			// go through them ancient demes
			for (int ancientDeme=0; ancientDeme<thisEpochLogMarginal.length; ancientDeme++ ){
				// scary
				if (!migrationParams.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
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

	private static Z H (double a, double b, Z u, Z lambda) {
		assert (u.re != Double.POSITIVE_INFINITY);
		assert (a >= 0d - ONELOCUS_EPSILON);
		assert (b > a);
		assert (a < b + ONELOCUS_EPSILON);
		// easy case
		if (Math.abs(a-b) < ONELOCUS_EPSILON) {
			return new Z(0d,0d);
		}
		
		// the real work
		Z result = null;
		// first decide on u
		if (u.re == Double.NEGATIVE_INFINITY) {
			// in this case it is just always zero
			result = new Z(0d, 0d);
		}
		else {
			// is lambda zero?
			if (Z.abs (lambda) < ONELOCUS_EPSILON) {
				// what about b?
				if (b == Double.POSITIVE_INFINITY) {
					// then we have infinity
					result = new Z(Double.POSITIVE_INFINITY, 0d);
				}
				else {
					// b not zero
					// compute exponential
					Z expU = ComplexMath.exponential (u);
					// and multiply it
					result = (new Z ()).Times (b - a, expU);
				}
			}
			else {
				// lambda not zero
				// what with b
				// we need the fraction in both cases
				Z fraction = null;
				try {
					fraction = (new Z()).Div(new Z(1d,0d), lambda);
				} catch (JampackException e1) {
					// no division by zero
					e1.printStackTrace();
					assert (false);
				}
				// also the Au guy
				Z exponentAu = (new Z()).Plus((new Z()).Times (a, lambda), u);
				Z expAu = ComplexMath.exponential (exponentAu);
				// prepare the bracket
				Z bracket = (new Z()).Times(-1d, expAu);
				
				// the value of the other exponential depends on b
				// if it is infinty, we do not alter the bracket
				// but we want to check the lambda
				assert (!(b == Double.POSITIVE_INFINITY) || (lambda.re < 0d + ONELOCUS_EPSILON));
				if (b != Double.POSITIVE_INFINITY) {
					// the completely normal case
					// first exponent
					Z exponentBu = (new Z()).Plus((new Z()).Times (b, lambda), u);
					// do the exponential
					Z expBu = ComplexMath.exponential (exponentBu);
					// and add it to the bracket
					bracket = (new Z()).Plus(bracket, expBu);
				}

				// get the result
				result = (new Z()).Times(fraction, bracket);
			}
		}
		// give it away now
		return result;
	}

	private static <H extends FSAHaplotype> double W (int startAncientDeme, int firstAbsorbInterval, int firstAbsorbAncientDeme, int secondAbsorbInterval, int secondAbsorbAncientDeme, double recombinationRate, StructuredMigrationParameters migrationParams, ArrayList<ArrayList<double[][]>> Q, ArrayList<Z[][][]> H, HashMap<List<Integer>, Double> rCache) {
		// some assertions
		// they can be everything
		assert ((0 <= firstAbsorbAncientDeme) && (firstAbsorbAncientDeme < migrationParams.numAncientDemes(firstAbsorbInterval)));
		assert ((0 <= secondAbsorbAncientDeme) && (secondAbsorbAncientDeme < migrationParams.numAncientDemes(secondAbsorbInterval)));

		// if they are in the wrong order, switch them
		if (secondAbsorbInterval < firstAbsorbInterval) {
			// call doubleU with changed roles
			return W (startAncientDeme, secondAbsorbInterval, secondAbsorbAncientDeme, firstAbsorbInterval, firstAbsorbAncientDeme, recombinationRate, migrationParams, Q, H, rCache);
		}
		
		// normal stuff
		int numAncientDemes = migrationParams.numAncientDemes(firstAbsorbInterval);
		assert ((0 <= startAncientDeme) && (startAncientDeme < numAncientDemes));
		if (firstAbsorbInterval == secondAbsorbInterval) {
			// in this case we just have R
			return R (startAncientDeme, firstAbsorbInterval, numAncientDemes + firstAbsorbAncientDeme, numAncientDemes + secondAbsorbAncientDeme, recombinationRate, migrationParams, H, rCache);
		}
		else {
			assert (firstAbsorbInterval < migrationParams.numIntervals()); 
			// we have firstAbsorbInterval < secondAbsorbInterval, so we need to multiply with Q
			double result = 0d;
			
			// loop over next start deme
			for (int nextStartAncientDeme=0; nextStartAncientDeme<migrationParams.numAncientDemes(firstAbsorbInterval+1); nextStartAncientDeme++) {
				double tmp = 0d;
				
				// loop over content of next start deme
				for (int memberAncientDeme : migrationParams.getDemography().getMemberDemesIndices (nextStartAncientDeme, firstAbsorbInterval+1)) {
					// is the index proper?
					assert (memberAncientDeme <= migrationParams.numAncientDemes(firstAbsorbInterval));
					// add some stuff
					tmp += R (startAncientDeme, firstAbsorbInterval, numAncientDemes + firstAbsorbAncientDeme, memberAncientDeme, recombinationRate, migrationParams, H, rCache);
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

	/// end demes can be non-absorbing and absorbing
	private static <H extends FSAHaplotype> double R (int startAncientDeme, int interval, int firstEndAncientDeme, int secondEndAncientDeme, double recombinationRate, StructuredMigrationParameters migrationParams, ArrayList<Z[][][]> H, HashMap<List<Integer>, Double> cache) {
		
		// if it is a pulse, there is no time for recombination, so no probability
		if (migrationParams.trunk.getRefinedDemography().isPulse(interval)) return 0d;
		
		// no pulse, so standard recombination business
		List<Integer> currKey = new ArrayList<Integer>();
		for (int i : (new Integer[] {startAncientDeme, interval, firstEndAncientDeme, secondEndAncientDeme})) {
			currKey.add(i);
		}
		
		if (cache.containsKey(currKey)) {
			return cache.get(currKey);
		}
		
		
		// some assertions
		int numAncientDemes = migrationParams.numAncientDemes(interval);
		assert ((0 <= startAncientDeme) && (startAncientDeme < numAncientDemes));
		// they can be everything
		assert ((0 <= firstEndAncientDeme) && (firstEndAncientDeme < 2*numAncientDemes));
		assert ((0 <= secondEndAncientDeme) && (secondEndAncientDeme < 2*numAncientDemes));
		
		// get the relevant eigenstuff
		MyEigenDecomposition eigenStuff = migrationParams.getEigenStuffList().get(interval);
		Z[][][] thisH = H.get(interval);
		
		Z preResult = null;
		
		// init result
		preResult = new Z(0d, 0d);
		// loop over recombinaton deme (that's eta)
		for (int recombinationAncientDeme=0; recombinationAncientDeme<numAncientDemes; recombinationAncientDeme++) {
			// now loop over eigen indices
			for (int k=0; k<eigenStuff.lambda.length; k++) {
				for (int m=0; m<eigenStuff.lambda.length; m++) {
					for (int n=0; n<eigenStuff.lambda.length; n++) {
						// get the right value
						Z summand = (new Z()).Times (eigenStuff.VW[k].get0(startAncientDeme, recombinationAncientDeme), eigenStuff.VW[m].get0(recombinationAncientDeme, firstEndAncientDeme));
						summand = (new Z()).Times (summand, eigenStuff.VW[n].get0(recombinationAncientDeme, secondEndAncientDeme));
						summand = (new Z()).Times (summand, thisH[k][m][n]);
						// add it
						preResult = (new Z()).Plus(preResult, summand);
					}
				}
			}
		}
		
		double toReturn = recombinationRate * preResult.re;
		cache.put(currKey, toReturn);
		
		// return a proper result
		assert (Math.abs (preResult.im) < ONELOCUS_EPSILON);
		// return multiplied with recombination rate
		assert (toReturn > 0 - ONELOCUS_EPSILON);
		return toReturn;
	}
	
	// [epoch][ancient deme]
	// returns the joint probability of being absorbed into the given (epoch,ancient deme) and not recombining between the two loci
	protected ArrayList<double[]> computeAncientDemeNoRecoJointLog (double recombinationRate) {

		// get p
		ArrayList<ArrayList<double[][]>> p = demoTransitions.p;

		// first precompute the H's (the integrals)
		// [interval][eigenIndex]
		ArrayList<Z[]> H = new ArrayList<Z[]>();
		// loop through intervals
		for (int interval=0; interval<migrationParams.numIntervals(); interval++) {
			if (this.refinedDemography.isPulse(interval)) {
				// pulse has no integral
				H.add(null);
			}
			else {
				// usual non-pulse integral
				// local eigenStuff
				MyEigenDecomposition eigenStuff = migrationParams.getEigenStuffList().get(interval);
				// get a new array
				Z[] currH = new Z[eigenStuff.lambda.length];
				// fill it
				for (int eigIdx=0; eigIdx<eigenStuff.lambda.length; eigIdx++) {
					Z complexArgLambda = (new Z()).Minus (eigenStuff.lambda[eigIdx], new Z(recombinationRate,0d));
					Z complexArgU = (new Z()).Times (- migrationParams.getIntervalList()[interval].startPoint, eigenStuff.lambda[eigIdx]);
					currH[eigIdx] = H  (migrationParams.getIntervalList()[interval].startPoint, migrationParams.getIntervalList()[interval].endPoint, complexArgU, complexArgLambda);
				}
				// remember currI
				H.add (currH);
			}
		}

		// now calculate the actual y probabilities
		ArrayList<double[]> logJointY = new ArrayList<double[]>();
		
		// loop over all epochs
		for (int epoch = 0; epoch < migrationParams.numIntervals(); epoch++){
			
			//make an array big enough
			double[] thisEpochY = new double[migrationParams.numAncientDemes(epoch)];
			// add it to array list
			logJointY.add(thisEpochY);
		
			// go through them ancient demes
			for (int ancientDeme=0; ancientDeme<thisEpochY.length; ancientDeme++ ){
		
		
				// scary ghost demes
				// this non-absorbing business should also take care of the pulse intervals
				if (migrationParams.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
					// log of zero
					thisEpochY[ancientDeme] = Double.NEGATIVE_INFINITY;
					// end next one
					continue;
				}
	
				// storage for computing y
				Z y = new Z (0d, 0d);
			
				// remember, it has to be absorbing, so + num demes
				int thisAbsorbDeme = migrationParams.numAncientDemes(epoch) + ancientDeme;
				MyEigenDecomposition eigenStuff = migrationParams.getEigenStuffList().get(epoch);
				// get all H's for the interval
				Z[] hList = H.get(epoch);
	
				// now first over eigen index
				for (int k=0; k<eigenStuff.lambda.length; k++) {
					// compute lambda times H
					Z lambdaH = (new Z()).Times (eigenStuff.lambda[k], hList[k]);
					// loop over starting demes of thisInterval
					for (int intervalStartAncientDeme=0; intervalStartAncientDeme < migrationParams.numAncientDemes(epoch); intervalStartAncientDeme++) {
						// calculate the corresponding summand
						Z summand = (new Z()).Times (p.get(0).get(epoch)[observedPresentDeme][intervalStartAncientDeme], (new Z()).Times (lambdaH, eigenStuff.VW[k].get0(intervalStartAncientDeme, thisAbsorbDeme)));
						// add the summand (actually, we would need only the real part, but to check it's nice)
						y = (new Z()).Plus(y, summand);
					}
				}
				
				// make sure y is real
				assert (Math.abs(y.im) < ONELOCUS_EPSILON);
				// should be positive
				assert (y.re > 0d - ONELOCUS_EPSILON);
				// be nice
				if (y.re < 0) y.re = 0d;
				
				// logarithmify the real part
				thisEpochY[ancientDeme] = Math.log(y.re);

				
				assert(!Double.isNaN(thisEpochY[ancientDeme]));
			}
		}
		
		// give it away now
		return logJointY;
	}

	// [prev epoch][next epoch][prev ancient deme][next ancient deme]
	// returns the joint probability of being absorbed into the given (prev epoch,prev ancient deme) at the first locus and recombining and eing absorbed into the given (next epoch,next ancient deme) at the second locus
	protected ArrayList<ArrayList<double[][]>> computeAncientDemeRecoJointLog (double recombinationRate) {

		// get some stuff
		ArrayList<ArrayList<double[][]>> p = demoTransitions.p;
		ArrayList<ArrayList<double[][]>> Q = demoTransitions.Q;
		
		// first precompute the H's (the ones we really need in our special case of equating intervals and epochs)
		// [interval][k][m][n]
		ArrayList<Z[][][]> H = new ArrayList<Z[][][]>();

		// loop through intervals
		for (int interval=0; interval<migrationParams.numIntervals(); interval++) {
			if (this.refinedDemography.isPulse(interval)) {
				// pulse has no integral
				H.add(null);
			}
			else {
				// usual non-pulse integral pre-computation
				// local eigenStuff
				MyEigenDecomposition eigenStuff = migrationParams.getEigenStuffList().get(interval);
				// get a new array
				Z[][][] currH = new Z[eigenStuff.lambda.length][eigenStuff.lambda.length][eigenStuff.lambda.length];
				// fill it (by looping through all three eigenindices)
				for (int k=0; k<eigenStuff.lambda.length; k++) {
					for (int m=0; m<eigenStuff.lambda.length; m++) {
						for (int n=0; n<eigenStuff.lambda.length; n++) {
							// remember to use the modified times function
							// first lambda
							Z tmp =  (new Z()).Minus (eigenStuff.lambda[k], eigenStuff.lambda[m]);
							tmp = (new Z()).Minus (tmp, eigenStuff.lambda[n]);
							Z complexArgLambda = (new Z()).Minus (tmp, new Z(recombinationRate,0d));
							// then u, special times
							tmp = (new Z()).Plus (eigenStuff.lambda[m], eigenStuff.lambda[n]);
							tmp = ComplexMath.times (migrationParams.getIntervalList()[interval].endPoint, tmp);
							Z complexArgU = (new Z()).Minus(tmp, ComplexMath.times(migrationParams.getIntervalList()[interval].startPoint, eigenStuff.lambda[k]));
							// calculate H
							currH[k][m][n] = H  (migrationParams.getIntervalList()[interval].startPoint, migrationParams.getIntervalList()[interval].endPoint, complexArgU, complexArgLambda);
						}
					}
				}
				// remember currI
				H.add (currH);
			}
		}
		// should be done for now

		// we have a function for R and W, and precomputed Q
		// maybe pre calculate R at some point
				
		// so now we can calculate z
		ArrayList<ArrayList<double[][]>> logJointZ = new ArrayList<ArrayList<double[][]>>();
		
		HashMap<List<Integer>, Double> rCache = new HashMap<List<Integer>, Double>();
		
		// go through previous epochs
		for (int prevAbsorbEpoch = 0; prevAbsorbEpoch < migrationParams.numIntervals(); prevAbsorbEpoch++){
			
			// make some storage
			ArrayList<double[][]> prevEpochZ = new ArrayList<double[][]>();
			logJointZ.add(prevEpochZ);
			
			// go through next epochs
			for (int nextAbsorbEpoch = 0; nextAbsorbEpoch < migrationParams.numIntervals(); nextAbsorbEpoch++){
			
				// make some storage
				double[][] currEpochsZ = new double[migrationParams.numAncientDemes(prevAbsorbEpoch)][migrationParams.numAncientDemes(nextAbsorbEpoch)];
				prevEpochZ.add(currEpochsZ);
				
				// go through ancient demes
				for (int prevAbsorbDeme=0; prevAbsorbDeme<currEpochsZ.length; prevAbsorbDeme++ ){				
					for (int nextAbsorbDeme=0; nextAbsorbDeme<currEpochsZ[prevAbsorbDeme].length; nextAbsorbDeme++ ){				
		
						// scary ghost demes (all ghosts zero)
						// this should also take care of pulse intervals
						if (migrationParams.trunk.isNonAbsorbing(prevAbsorbEpoch, prevAbsorbDeme, UberDemographyCore.ONELOCUS_EPSILON) || migrationParams.trunk.isNonAbsorbing(nextAbsorbEpoch, nextAbsorbDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
							// log of zero
							currEpochsZ[prevAbsorbDeme][nextAbsorbDeme] = Double.NEGATIVE_INFINITY;
							// end next one
							continue;
						}
		
						// which is the lowest epoch?
						int minAbsorbEpoch = Math.min(prevAbsorbEpoch, nextAbsorbEpoch);
						assert (minAbsorbEpoch <= migrationParams.numIntervals());
						
						// init
						double preZ = 0d;
						// compute it
						// first loop over recombination happening in every interval up to minAbsorbInterval (not including)
						for (int recoInterval=0; recoInterval<minAbsorbEpoch; recoInterval++) {
							assert (recoInterval < migrationParams.numIntervals());
							double summand = 0d;
							// loop over reco interval start deme
							for (int recoIntervalStartDeme=0; recoIntervalStartDeme<migrationParams.numAncientDemes(recoInterval); recoIntervalStartDeme++) {
								// loop over start deme of next interval on previuos locus
								for (int prevLocusNextIntervalStartDeme=0; prevLocusNextIntervalStartDeme<migrationParams.numAncientDemes(recoInterval+1); prevLocusNextIntervalStartDeme++) {
									// now loop over the ending deme of reco interval at previuos locus
									for (int prevLocusMemberDeme : migrationParams.getDemography().getMemberDemesIndices (prevLocusNextIntervalStartDeme, recoInterval+1)) {
										// loop over start deme of next interval on current locus
										for (int nextLocusNextIntervalStartDeme=0; nextLocusNextIntervalStartDeme<migrationParams.numAncientDemes(recoInterval+1); nextLocusNextIntervalStartDeme++) {
											// now loop over the ending deme of reco interval at current locus
											for (int nextLocusMemberDeme : migrationParams.getDemography().getMemberDemesIndices (nextLocusNextIntervalStartDeme, recoInterval+1)) {
												// and prepare adding something
												// p [first interval][second interval][first start deme][second start deme]
												double tmp = p.get(0).get(recoInterval)[observedPresentDeme][recoIntervalStartDeme];
												tmp *= R (recoIntervalStartDeme, recoInterval, prevLocusMemberDeme, nextLocusMemberDeme, recombinationRate, migrationParams, H, rCache);
												// Q [first interval][absorbing interval][first start deme][absorbing deme]
												tmp *= Q.get(recoInterval+1).get(prevAbsorbEpoch)[prevLocusNextIntervalStartDeme][prevAbsorbDeme];
												tmp *= Q.get(recoInterval+1).get(nextAbsorbEpoch)[nextLocusNextIntervalStartDeme][nextAbsorbDeme];
												
												// add it
												summand += tmp;
											}
										}
									}						
								}
							}
							// add the summand to our result
							preZ += summand;
						}
						
						// and then add the term for recombination in same interval as minAbsorbInterval
						double summand = 0d;
						// loop over reco interval start deme
						for (int recoIntervalStartAncientDeme=0; recoIntervalStartAncientDeme<migrationParams.numAncientDemes(minAbsorbEpoch); recoIntervalStartAncientDeme++) {
							// p [first interval][second interval][first start deme][second start deme]
							double tmp = p.get(0).get(minAbsorbEpoch)[observedPresentDeme][recoIntervalStartAncientDeme];
							tmp *= W (recoIntervalStartAncientDeme, prevAbsorbEpoch, prevAbsorbDeme, nextAbsorbEpoch,nextAbsorbDeme, recombinationRate, migrationParams, Q, H, rCache);
							// remember tmp
							summand += tmp;
						}
						// add the summand
						preZ += summand;
						
						// store the normalized log
						assert (preZ > 0 - ONELOCUS_EPSILON);
						// have to be nice
						if (preZ < 0) preZ = 0.0;
						
						// take the log
						currEpochsZ[prevAbsorbDeme][nextAbsorbDeme] = Math.log(preZ);
												
						assert(!Double.isNaN(currEpochsZ[prevAbsorbDeme][nextAbsorbDeme]));
					}
				}
			}
		}
		
		//give back the calculated transitions
		return logJointZ;
	}
	
	// [epoch][ancient deme][trunk type][emission type]
	// This gives marginal probability of absorbing inside (epoch, ancient deme), times conditional probability of emission type given trunk type
	protected ArrayList<double[][][]> computeAncientDemeEmissionLogProbabilities (double[] mutationRate, MyEigenDecomposition eigenMutation) {

		assert(mutationRate.length == migrationParams.numIntervals());
		
		// get p
		ArrayList<ArrayList<double[][]>> p = demoTransitions.p;

		// first precompute the H's (the integrals)
		// [interval][migrationEigenIndex][mutationEigenIndex]
		ArrayList<Z[][]> H = new ArrayList<Z[][]>();
		// loop over intervals
		for (int absorbInterval=0; absorbInterval<migrationParams.numIntervals(); absorbInterval++) {
			// no emission from pulse interval (no integral pre-computation)
			if (this.refinedDemography.isPulse(absorbInterval)) {
				H.add (null);
			}
			else {
				// usual non-pulse integral pre-computation
				// get some eigenstuff
				MyEigenDecomposition eigenMigration = migrationParams.getEigenStuffList().get(absorbInterval);
				// get an H
				Z[][] currH = new Z[eigenMigration.lambda.length][eigenMutation.lambda.length];
				
				// fill the H
				// loop over eigenMigration index
				for (int k=0; k<eigenMigration.lambda.length; k++) {
					// loop over eigenMutation index
					for (int j=0; j<eigenMutation.lambda.length; j++) {
						// and get the H
						Z tmp = (new Z()).Times (mutationRate[absorbInterval], (new Z()).Minus(eigenMutation.lambda[j], new Z(1d,0d)));
						Z complexArgLambda = (new Z()).Plus (tmp, eigenMigration.lambda[k]);
						Z complexArgU = (new Z()).Times (- migrationParams.getIntervalList()[absorbInterval].startPoint, eigenMigration.lambda[k]);
						currH[k][j] = H  (migrationParams.getIntervalList()[absorbInterval].startPoint, migrationParams.getIntervalList()[absorbInterval].endPoint, complexArgU, complexArgLambda);
					}
				}
				
				// remember H
				H.add(currH);
			}
		}
		// done
		
			
		int numAlleles = eigenMutation.lambda.length;
		// now get some logarithmic emmisionProbabilties
		// [epoch][ancient deme][trunk type][emission type]
		ArrayList<double[][][]> logJointEmission = new ArrayList<double[][][]>();
		
		// fill it
		// loop over epoch
		for (int epoch = 0; epoch < migrationParams.numIntervals(); epoch++){
			// make storage and add it
			double[][][] thisEmission = new double[migrationParams.numAncientDemes(epoch)][numAlleles][numAlleles];
			logJointEmission.add (thisEmission);
			
			// loop over ancient deme
			for (int ancientDeme=0; ancientDeme < thisEmission.length; ancientDeme++) {
				
				// loop over emission type
				for (int emissionType=0; emissionType<numAlleles; emissionType++) {
	
					// scary ghost demes (all ghosts zero, just to be sure)
					// this should deal with pulse stuff too
					if (migrationParams.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
						// all the trunk types
						for (int trunkType=0; trunkType<numAlleles; trunkType++) {
							// log of zero
							thisEmission[ancientDeme][trunkType][emissionType] = Double.NEGATIVE_INFINITY;
						}
						// end next one
						continue;
					}
	
					// get some stuff
					MyEigenDecomposition eigenMigration = migrationParams.getEigenStuffList().get(epoch);
					int numAncientDemes = migrationParams.numAncientDemes(epoch);

					// loop over trunk type
					for (int trunkType=0; trunkType<numAlleles; trunkType++) {
						// init to win it
						Z preResult = new Z(0d,0d);
						// loop over start deme of absorption interval
						for (int absorbIntervalStartAncientDeme=0; absorbIntervalStartAncientDeme<numAncientDemes; absorbIntervalStartAncientDeme++) {
							Z preSum =new Z(0d,0d);
							// loop over eigenMigration index
							for (int k=0; k<eigenMigration.lambda.length; k++) {
								// loop over eigenMutation index
								for (int j=0; j<eigenMutation.lambda.length; j++) {
									// now add something
									// we have to use our special times
									Z tmp = ComplexMath.ourTimes(eigenMigration.lambda[k], H.get(epoch)[k][j]);
									// normal times here
									tmp = (new Z()).Times (tmp, eigenMutation.VW[j].get0(trunkType, emissionType));
									// remember, the absorbing states have to be indexed with + numAncientDemes
									tmp = (new Z()).Times (tmp, eigenMigration.VW[k].get0(absorbIntervalStartAncientDeme, numAncientDemes + ancientDeme));
									// add it to pre sum
									preSum = (new Z()).Plus(preSum, tmp);
								}
							}
							
							//multiply by p
							// p [first interval][second interval][first start deme][second start deme]
							preSum = (new Z()).Times (p.get(0).get(epoch)[observedPresentDeme][absorbIntervalStartAncientDeme], preSum);
							// and add to preResult
							preResult = (new Z()).Plus(preResult, preSum);
						}
						
						// store normalized log probability
						assert (Math.abs(preResult.im) < ONELOCUS_EPSILON);
						assert (preResult.re > 0 - ONELOCUS_EPSILON);
						// have to be nice
						if (preResult.re < 0) preResult.re = 0.0;
	
						// take the log
						thisEmission[ancientDeme][trunkType][emissionType] = Math.log(preResult.re);

						
						assert(!Double.isNaN(thisEmission[ancientDeme][trunkType][emissionType]));
					}
				}
			}
		}
		
		// give it away now
		return logJointEmission;
		
	}	
	
	// [demoState]
	// returns probability of being absorbed into the demoState (summed over all states living in the demoState)
	private double[] computeDemoStateLogMarginal() {
		
		double[] demoStateLogMarginal = new double[this.numDemoStates()];
		// fill with log 0s
		Arrays.fill(demoStateLogMarginal, Double.NEGATIVE_INFINITY);
		
		// add up the joint probabilities
		for (int epoch = 0; epoch < migrationParams.numIntervals(); epoch++){
			for (int ancientDeme = 0; ancientDeme < migrationParams.numAncientDemes(epoch); ancientDeme++){
					
				if (this.migrationParams.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
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
		
		if (UberDemographyCore.RENORMALIZE_PROBS) EigenCore.renormalizeLogStochasticVector(demoStateLogMarginal, this.renormalizeEpsilon);
		
		return demoStateLogMarginal;

	}
	
	// [demoState]
	// returns conditional probability of no recombination given current demoState
	private double[] computeLogNoReco (double recoRate) {
		
		ArrayList<double[]> epochAncientDemeNoReco = computeAncientDemeNoRecoJointLog (recoRate);

		double[] demoStateNoRecoTransition = new double[this.numDemoStates()];
		// fill with log 0s
		Arrays.fill(demoStateNoRecoTransition, Double.NEGATIVE_INFINITY);
		
		// add up the joint probabilities
		for (int epoch = 0; epoch < migrationParams.numIntervals(); epoch++){
			for (int ancientDeme = 0; ancientDeme < migrationParams.numAncientDemes(epoch); ancientDeme++){
				
				if (this.migrationParams.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
					assert (epochAncientDemeNoReco.get(epoch)[ancientDeme] == Double.NEGATIVE_INFINITY);
					continue;
				}
				
				for (int presentDeme = 0; presentDeme < this.trunk.getSampleSizes().length; presentDeme++) {
					
					int demoState = this.getDemoStateCollection().getDemoState(epoch, ancientDeme, presentDeme);
					
					// joint probability of event
					double value = epochAncientDemeNoReco.get(epoch)[ancientDeme];
					
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
			assert(this.trunk.sampleSizes[currDeme] != 0 || (Double.NEGATIVE_INFINITY == this.demoStateLogMarginal[demoState]));
			assert(this.trunk.sampleSizes[currDeme] != 0 || (Double.NEGATIVE_INFINITY == demoStateNoRecoTransition[demoState]));

			
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
	private double[][] computeLogReco (double recoRate) {
		
		ArrayList<ArrayList<double[][]>> epochAncientDemeReco = computeAncientDemeRecoJointLog (recoRate);
				
		double[][] demoStateRecoTransition = new double[this.numDemoStates()][this.numDemoStates()];
		// fill with log 0s
		for (int i=0; i<demoStateRecoTransition.length; i++) {
			Arrays.fill(demoStateRecoTransition[i], Double.NEGATIVE_INFINITY);
		}
		
		// add up the joint probabilities for a single lineage residing in each interval/deme
		// loop through previuos refined intervals, ancient demes
		for (int prevEpoch = 0; prevEpoch < migrationParams.numIntervals(); prevEpoch++){
			for (int prevAncientDeme = 0; prevAncientDeme < migrationParams.numAncientDemes(prevEpoch); prevAncientDeme++){
				// loop through next refined intervals, ancient demes
				for (int nextEpoch = 0; nextEpoch < migrationParams.numIntervals(); nextEpoch++){
					for (int nextAncientDeme = 0; nextAncientDeme < migrationParams.numAncientDemes(nextEpoch); nextAncientDeme++){
					
						// skip ghost demes
						if (this.migrationParams.trunk.isNonAbsorbing(prevEpoch, prevAncientDeme, UberDemographyCore.ONELOCUS_EPSILON)|| this.migrationParams.trunk.isNonAbsorbing(nextEpoch, nextAncientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
							assert (epochAncientDemeReco.get(prevEpoch).get(nextEpoch)[prevAncientDeme][nextAncientDeme] == Double.NEGATIVE_INFINITY);
							continue;
						}

						for (int prevPresentDeme = 0; prevPresentDeme < this.trunk.getSampleSizes().length; prevPresentDeme++){
							for (int nextPresentDeme = 0; nextPresentDeme < this.trunk.getSampleSizes().length; nextPresentDeme++){

								int prevDemoState = this.getDemoStateCollection().getDemoState(prevEpoch, prevAncientDeme, prevPresentDeme);
								int nextDemoState = this.getDemoStateCollection().getDemoState(nextEpoch, nextAncientDeme, nextPresentDeme);
								
								// joint probability of event
								double value = epochAncientDemeReco.get(prevEpoch).get(nextEpoch)[prevAncientDeme][nextAncientDeme];
								
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

				assert(this.trunk.sampleSizes[prevDeme] != 0 || (Double.NEGATIVE_INFINITY == this.demoStateLogMarginal[prevDemoState]));
				assert(this.trunk.sampleSizes[prevDeme] != 0 || (Double.NEGATIVE_INFINITY == demoStateRecoTransition[prevDemoState][nextDemoState]));
				
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
	
	@Override
	protected RecoLogProbs computeRecoLogProbs(double recoRate) {
				
		RecoLogProbs toReturn = new RecoLogProbs (this.computeLogReco(recoRate), this.computeLogNoReco(recoRate));
		if (UberDemographyCore.RENORMALIZE_PROBS) EigenCore.renormalizeLogTransitions(toReturn, this.demoStateLogMarginal, this.renormalizeEpsilon);
		
		return toReturn;
	}
	
	// [demoState][trunkType][emissionType]
	// returns conditional probability of emission type, given absorption into a lineage with trunk type
	@Override
	protected double[][][] computeLogEmission (double[] mutationRate, MyEigenDecomposition eigenMutation) {
		// If there is only one mutation rate, just use that everywhere
		
		double[] refinedMutationRate = this.getRefinedTheta(mutationRate);
		// If there is more than one mutation rate, we must correctly partition it
		
		ArrayList<double[][][]> refinedIntervalAncientDemeEmissionLogProbs = computeAncientDemeEmissionLogProbabilities(refinedMutationRate, eigenMutation);
		
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
		for (int epoch = 0; epoch < migrationParams.numIntervals(); epoch++){
			for (int ancientDeme = 0; ancientDeme < migrationParams.numAncientDemes(epoch); ancientDeme++){
					
				for (int presentDeme = 0; presentDeme < this.trunk.getSampleSizes().length; presentDeme++){
					
					for (int emissionType = 0; emissionType < numAlleles; emissionType++){
						for (int trunkType = 0; trunkType < numAlleles; trunkType++){
							
							if (this.migrationParams.trunk.isNonAbsorbing(epoch, ancientDeme, UberDemographyCore.ONELOCUS_EPSILON)) {
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
					assert(this.trunk.sampleSizes[currDeme] != 0 || (Double.NEGATIVE_INFINITY == this.demoStateLogMarginal[demoState]));
					assert(this.trunk.sampleSizes[currDeme] != 0 || (Double.NEGATIVE_INFINITY == demoStateEmission[demoState][trunkType][emissionType]));

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
				if (UberDemographyCore.RENORMALIZE_PROBS) EigenCore.renormalizeLogStochasticMatrix(demoStateEmission[demoState], this.renormalizeEpsilon);
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
		ArrayList<double[][][]> ancientDemeEmissions = computeAncientDemeEmissionLogProbabilities(mutationRate, eigenMutation);
		
		int numAlleles = eigenMutation.lambda.length;
		
		for (int epoch=0; epoch<this.migrationParams.numIntervals(); epoch++) {
			for (int ancientDeme=0; ancientDeme<this.migrationParams.numAncientDemes(epoch); ancientDeme++) {
				
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
		ArrayList<ArrayList<double[][]>> epochAncientDemeReco = computeAncientDemeRecoJointLog (recoRate);
		ArrayList<double[]> epochAncientDemeNoReco = computeAncientDemeNoRecoJointLog (recoRate);

		double[][] reco = computeLogReco(recoRate);
		double[] noReco = computeLogNoReco(recoRate);
		
		// check no reco
		
		double noRecoSumOld = 0d;
		double noRecoSumNew = 0d;
		
		for (int epoch=0; epoch<this.migrationParams.numIntervals(); epoch++) {
			for (int ancientDeme=0; ancientDeme<this.migrationParams.numAncientDemes(epoch); ancientDeme++) {
				noRecoSumOld += Math.exp(epochAncientDemeNoReco.get(epoch)[ancientDeme]);
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
		
		for (int prevEpoch=0; prevEpoch<this.migrationParams.numIntervals(); prevEpoch++) {
			for (int prevAncientDeme=0; prevAncientDeme<this.migrationParams.numAncientDemes(prevEpoch); prevAncientDeme++) {
				for (int nextEpoch=0; nextEpoch<this.migrationParams.numIntervals(); nextEpoch++) {
					for (int nextAncientDeme=0; nextAncientDeme<this.migrationParams.numAncientDemes(nextEpoch); nextAncientDeme++) {
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
		
		// check things
		for (int prevEpoch = 0; prevEpoch < this.migrationParams.numIntervals(); prevEpoch++){
			for (int prevAncientDeme = 0; prevAncientDeme < this.migrationParams.numAncientDemes(prevEpoch); prevAncientDeme++){
				// this has to sum to one				
				double transitionSum = Math.exp(epochAncientDemeNoReco.get(prevEpoch)[prevAncientDeme]);
				
				for (int nextEpoch = 0; nextEpoch < this.migrationParams.numIntervals(); nextEpoch++){
					for (int nextAncientDeme = 0; nextAncientDeme < this.migrationParams.numAncientDemes(nextEpoch); nextAncientDeme++){
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
		// show who you are
		
		// but only if you are not empty
		if (this.conditionalConfigEmpty) {
			System.out.println("[NO CORE, SINCE CONDITIONAL CONFIG EMPTY]");
			return;
		}
		
		// marginal
		outStream.println("=== Marginal probabilities:");
		double logSum = Double.NEGATIVE_INFINITY;
		for (double value : this.demoStateLogMarginal) {
			outStream.print (Math.exp (value) + "\t");
			logSum = LogSum.computePairLogSum(logSum, value);
		}
		outStream.println();
		outStream.println("=== logSum: " + logSum);
		
		// no reco & reco stuff
		double[] logNoReco = this.computeLogNoReco(recoRate);
		double[][] logReco = this.computeLogReco(recoRate); 
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
		return this.observedPresentDeme;
	}
	
	@Override
	public TrunkProcess getTrunk() {
		return this.trunk;
	}
	
	private final TrunkProcess trunk;
	protected final StructuredMigrationParameters migrationParams;
	final private int observedPresentDeme;
	public final Demography refinedDemography;
	
	// to deal with fully or partially empty configurations
	private final boolean conditionalConfigEmpty;
	
	// some precomputed stuff
	private final DemographyTransitions demoTransitions; 
	
	// initial [refined interval][ancient deme] probabilities, for a single lineage
	protected final ArrayList<double[]> epochAncientDemeMarginalLogProbabilities;
	
	private final double[] demoStateLogMarginal;
	
}
