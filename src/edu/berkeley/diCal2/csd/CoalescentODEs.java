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
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.HighamHall54Integrator;

import edu.berkeley.diCal2.csd.CoalescentODEs.TimeDependentDouble.BackwardExponentialDoublePlus;
import edu.berkeley.diCal2.csd.CoalescentODEs.TimeDependentDouble.ConstantDouble;
import edu.berkeley.diCal2.csd.CoalescentODEs.TimeDependentDouble.DoubleBackwardExponentialDoublePlus;
import edu.berkeley.diCal2.utility.Pair;
import edu.berkeley.diCal2.utility.RateMatrixTools;
import edu.berkeley.diCal2.utility.RealPartition.Interval;
import edu.berkeley.diCal2.utility.SimplePair;

public class CoalescentODEs {

	public static final int MAX_NR_MUT_EVENTS = 1;
	public static final double INTEGRATION_LENGTH = 5.0d;
	public static final double ABS_TOLERANCE = 1.0e-14;
	public static final double REL_TOLERANCE = 1.0e-14;
	private static final double MIN_STEP = 1.0e-10;
	private static final double MAX_TIME = 2000;
	private static final double FACTOR_START_MAX = 10;
	
	public static double[][] getMarginalTransitionMatrix (Interval epoch, double[][] migrationRates, double[] absorptionRates, double[] expGrowthRates) {
		// just check the dimensions
		assert (migrationRates.length == absorptionRates.length);
		assert (RateMatrixTools.isSquare(migrationRates));
		assert (absorptionRates.length == expGrowthRates.length);
		
		// get some numbers
		double startTime = epoch.startPoint;
		double endTime = epoch.endPoint;
		
		// and start with the ode
		int numDemes = expGrowthRates.length;
		// get the states right
		MarginalStatesDict stateDict = new MarginalStatesDict(numDemes);
		// build a matrix
		SparseMatrix<TimeDependentDouble> rateMatrix = new SparseMatrix<TimeDependentDouble>(stateDict.numStates());
		// and fill it
		for (int srcDeme=0; srcDeme<numDemes; srcDeme++) {
			for (int dstDeme=0; dstDeme<numDemes; dstDeme++) {
				// get them non-absorb indices
				int srcIdx = stateDict.getIndex(false, srcDeme);
				int dstIdx = stateDict.getIndex(false, dstDeme);
				// now see what to do
				if (srcDeme == dstDeme) {
					// we need another index
					int absorbIdx = stateDict.getIndex(true, dstDeme);
					// we have some negative diagonal, but also the absorption rate
					double constantPart = migrationRates[srcDeme][dstDeme];
					assert (constantPart <= 0d);
					// the given rates are backward in time growth rates, so since we want the inverse, we have to take the negative rate
					// no inverse of rates, cause is already an absorpton rate
					TimeDependentDouble diagonalValue = new BackwardExponentialDoublePlus (constantPart, - absorptionRates[dstDeme], endTime, - expGrowthRates[dstDeme]);
					rateMatrix.set (srcIdx, dstIdx, diagonalValue);
					// and the absorption rate
					TimeDependentDouble absRate = new BackwardExponentialDoublePlus (0d, absorptionRates[dstDeme], endTime, - expGrowthRates[dstDeme]);
					rateMatrix.set (srcIdx, absorbIdx, absRate);					
				}
				else {
					// just normal transition
					rateMatrix.set (srcIdx, dstIdx, new ConstantDouble(migrationRates[srcDeme][dstDeme]));
				}
			}
		}
		assert (isRateMatrix(rateMatrix, (startTime + endTime)/2));
		
		
		// init return container
		double[][] marginalMatrix = new double[2*numDemes][2*numDemes];
		// just be sure
		for (int i=0; i<marginalMatrix.length; i++) {
			Arrays.fill (marginalMatrix[i], 0d);
		}
		
		// run ode (a few times)
		for (int startDeme=0; startDeme<numDemes; startDeme++) {
			// prepare the initial vector
			double[] y_start = new double[stateDict.numStates()];
			Arrays.fill(y_start, 0d);
			y_start[startDeme] = 1d;
			// and run it
			double[] y_end = solveODEHighamHall (y_start, rateMatrix, startTime, endTime);
			assert (Math.abs (sum(y_end) - 1d) < UberDemographyCore.ONELOCUS_EPSILON);
			
			// get stuff out
			assert (marginalMatrix[startDeme].length == y_end.length);
			for (int j=0; j<y_end.length; j++) {
				marginalMatrix[startDeme][j] = y_end[j];
			}
		}
		// for completeness sake
		for (int i=numDemes; i<2*numDemes; i++) {
			marginalMatrix[i][i] = 1d;
		}
		
		
		// things has to sum to one
		assert (RateMatrixTools.isProperStochaticMatrix(marginalMatrix, UberDemographyCore.ONELOCUS_EPSILON));
		// and some more checking
		double tmpOneLocusEpsilon = UberDemographyCore.ONELOCUS_EPSILON * 100;
		TreeSet<Integer> realDemes = new TreeSet<Integer>();
		for (int i=0; i<absorptionRates.length; i++) {
			if (absorptionRates[i] > tmpOneLocusEpsilon) {
				realDemes.add(i);
			}
		}
		for (int i=0; i<numDemes; i++) {
			// not real demes have to be zero
			if (epoch.endPoint != Double.POSITIVE_INFINITY) {
				for (int j=0; j<marginalMatrix[i].length; j++) {
					if ((j >= numDemes) && !realDemes.contains(j-numDemes)) {
						assert (marginalMatrix[i][j] <  tmpOneLocusEpsilon);
					}
				}
			}
			else {
				for (int j=0; j<2*numDemes; j++) {
					if ((j < numDemes) || !realDemes.contains(j-numDemes)) {
						assert (marginalMatrix[i][j] <  tmpOneLocusEpsilon);
					}
				}
			}
		}
			
		// return it
		return marginalMatrix;
	}

	
	public static double[] solveODEHighamHall (double[] y_start,	SparseMatrix<TimeDependentDouble> rateMatrix, double startTime, double endTime) {
		assert (y_start.length == rateMatrix.getDimension());
		// pass it through
		return solveODEHighamHall(y_start, new FancyODE (rateMatrix), startTime, endTime);
	}

	
	// this solver starts at startTime
	public static double[] solveODEHighamHall (double[] y_start, FirstOrderDifferentialEquations ode, double startTime, double endTime) {

//		int times = 6;
		
		// make a copy of the starting values, also gives initial stat
		double[] values = Arrays.copyOf (y_start, y_start.length);
		
		double minStep = MIN_STEP;
		double maxStep = (endTime == Double.POSITIVE_INFINITY ? INTEGRATION_LENGTH : endTime - startTime);
		double absTolerance = ABS_TOLERANCE;
		double relTolerance = REL_TOLERANCE;
		HighamHall54Integrator integrator = new RobustHighamHall(minStep, maxStep, absTolerance, relTolerance);
		
		// run the integrator
		integrator.integrate (ode, startTime, values, endTime == Double.POSITIVE_INFINITY ? startTime + INTEGRATION_LENGTH : endTime, values);
		
		// This will keep integrating until equilibrium is reached (derivatives are very small)
		if(endTime == Double.POSITIVE_INFINITY){

			// this should not be run forever
			double maxTime = Math.max (FACTOR_START_MAX*startTime, MAX_TIME);

			double[] currDeriv = new double[values.length];
			ode.computeDerivatives(startTime + INTEGRATION_LENGTH, values, currDeriv);
			boolean reRun = false;
			for(double deriv : currDeriv){
				if (deriv > UberDemographyCore.ONELOCUS_EPSILON){
					reRun = true;
				}
			}
			while(reRun){
				startTime = startTime + INTEGRATION_LENGTH;
				integrator.integrate(ode, startTime, values, startTime + INTEGRATION_LENGTH, values);
				ode.computeDerivatives(startTime + INTEGRATION_LENGTH, values, currDeriv);
				reRun = false;
				double maxDeriv = 0d;
				for(double deriv : currDeriv){
					if (Math.abs(deriv) > maxDeriv){
						maxDeriv = Math.abs(deriv);
					}
				}
				if (maxDeriv > .1*UberDemographyCore.ONELOCUS_EPSILON){	
//				if (maxDeriv > .01*UberDemographyCore.ONELOCUS_EPSILON){	
//				if (maxDeriv > .01*UberDemographyCore.ONELOCUS_EPSILON && (startTime < maxTime)){
					if (startTime > maxTime) {
						throw new RuntimeException ("The ODE solver computing the HMM probabilities is running too long (maxDeriv: " + maxDeriv + "). Try changing the paramters.");
					}
					reRun = true;
				}
//				else {
//					assert (maxDeriv < UberDemographyCore.ONELOCUS_EPSILON);
//				}
//				else {
//					MetaOptimization.synchronizedPrintln("startTime " + startTime + ", maxDeriv " +  maxDeriv + ", epsilon " + .01*UberDemographyCore.ONELOCUS_EPSILON);
//				}
//				if (startTime > maxTime && (times-- > 0)) {
//					MetaOptimization.synchronizedPrintln ("too long:  " + Arrays.toString(values) + ", maxDeriv " +  maxDeriv + ", epsilon " + .01*UberDemographyCore.ONELOCUS_EPSILON);
//				}
			}
		}
		
		// copy END into container to return it
		double[] y_final = new double[y_start.length];
		for (int i=0; i<y_final.length; i++) {
			y_final[i] = values[i];
			assert y_final[i] > -100*UberDemographyCore.ONELOCUS_EPSILON;
			if (y_final[i] < 0d) {
				y_final[i] = 0d;
			}
		}
		
		// give it away now
		return y_final;
	}

	/*
	 * This is an ODE solver, that works for ODE's that are _guaranteed_ to never be negative.
	 * If the ODE is negative, or approaches 0 from above too quickly -- then the ODE will take infinitely many tiny steps
	 * (so don't use it!!!)
	 */
	public static class PositiveODESolver { 
		
		private static final double REL_CHANGE_PER_STEP = 0.005;
		private static final double ABS_CHANGE_FROM_ZERO = 0.005;
		private static final double RETRY_STEP_FACTOR = .9;
		
		public static double[] solve(double[] y_start, FirstOrderDifferentialEquations ode, double startTime, double endTime) {
			
			assert endTime > startTime;
			assert endTime - startTime < Double.POSITIVE_INFINITY;
			
			double currTime = startTime;
			
			double[] y_curr = Arrays.copyOf(y_start, y_start.length);
			double[] y_dot1 = new double[y_curr.length];
			
			//Runge-Kutta things
			double[] y_dot2 = new double[y_curr.length];
			double[] y_2 = new double[y_curr.length];
			double[] y_dot3 = new double[y_curr.length];
			double[] y_3 = new double[y_curr.length];
			double[] y_dot4 = new double[y_curr.length];
			double[] y_4 = new double[y_curr.length];
			
			while (currTime < endTime) {
				ode.computeDerivatives(currTime, y_curr, y_dot1);
				
				double prevTime = currTime;
				
				//Approximate step size using only crude y_dot estimate
				double stepSize = Math.min(getStepSize(y_curr, y_dot1), .015);
				
				if (currTime + stepSize > endTime) {
					stepSize = currTime - prevTime;
					if(stepSize < UberDemographyCore.ONELOCUS_EPSILON){
						break;
					}
				}
				boolean tryStep = true;
				while(tryStep){
					tryStep = false;
					//improve y_dot estimate via Runge-Kutta bullshit

					for(int i = 0; i < y_curr.length; i++){
						y_2[i] = y_curr[i] + stepSize / 2d * y_dot1[i];
					}
					ode.computeDerivatives(currTime + stepSize / 2d, y_2, y_dot2);
					for(int i = 0; i < y_curr.length; i++){
						y_3[i] = y_curr[i] + stepSize / 2d * y_dot2[i];
					}
					ode.computeDerivatives(currTime + stepSize / 2d, y_3, y_dot3);
					for(int i = 0; i < y_curr.length; i++){
						y_4[i] = y_curr[i] + stepSize * y_dot3[i];
					}
					ode.computeDerivatives(currTime + stepSize, y_4, y_dot4);
				
					for(int i = 0; i < y_curr.length; i++){
						if(stepSize / 6d * (y_dot1[i] + 2 * y_dot2[i] + 2 * y_dot3[i] + y_dot4[i]) > 0 && Math.abs(stepSize / 6d * (y_dot1[i] + 2 * y_dot2[i] + 2 * y_dot3[i] + y_dot4[i])) > ABS_CHANGE_FROM_ZERO){
							tryStep = true;
						}
						if(stepSize / 6d * (y_dot1[i] + 2 * y_dot2[i] + 2 * y_dot3[i] + y_dot4[i]) < 0 && Math.abs((stepSize / 6d * (y_dot1[i] + 2 * y_dot2[i] + 2 * y_dot3[i] + y_dot4[i])) / y_curr[i]) > REL_CHANGE_PER_STEP){
							tryStep = true;
						}
					}
					//If error is too large, rescale the step size and try again
					if(tryStep == true){
						stepSize = RETRY_STEP_FACTOR * stepSize;
					}
				}
				currTime += stepSize;
				
				//Perform step
				for (int i=0; i<y_curr.length; i++) {
					y_curr[i] += stepSize / 6d * (y_dot1[i] + 2 * y_dot2[i] + 2 * y_dot3[i] + y_dot4[i]);
					assert(y_curr[i] > -UberDemographyCore.ONELOCUS_EPSILON && y_curr[i] < 1 + UberDemographyCore.ONELOCUS_EPSILON);
				}
			}
			
			return y_curr;
		}
		
		private static double getStepSize(double[] y_curr, double[] y_dot) {
			double stepSize = Double.POSITIVE_INFINITY;
			
			for (int i = 0; i < y_curr.length; i++) {
				assert y_curr[i] >= 0d;
				
				// if y_curr is 0, then y_dot must be positive
				assert y_curr[i] > 0d || y_dot[i] >= 0d;
				
				double currStepSize;
				if (y_curr[i] < ABS_CHANGE_FROM_ZERO && y_dot[i] > 0d) {
					// y_curr is close to 0 and increasing; want to take a step that takes us far enough away from zero
					// (otherwise, if y_curr is small and y_dot is big, then taking step size based on relative change would require very small steps)
					currStepSize = ABS_CHANGE_FROM_ZERO / y_dot[i];					
				} else {
					// either y_curr is far enough from 0, or we are decreasing.
					// Either way, we don't need to worry about taking tiny steps (y_dot / y_curr shouldn't be too large)

					double relChangePerUnitStep = Math.abs(y_dot[i] / y_curr[i]);
					if (Double.isNaN(relChangePerUnitStep)) relChangePerUnitStep = 0d;
					
					currStepSize = REL_CHANGE_PER_STEP / relChangePerUnitStep;
				}
				
				stepSize = Math.min(stepSize, currStepSize);
			}
			
			// changePerStep / stepSize = changePerUnitTime
			assert stepSize > 0d;
			return stepSize;
		}
	}
	
	// ode object
	public static class FancyODE implements FirstOrderDifferentialEquations {

		final private SparseMatrix<TimeDependentDouble> rateMatrix;
		
		public FancyODE (SparseMatrix<TimeDependentDouble> rateMatrix) {
			
			// remember values
			this.rateMatrix = rateMatrix;
		}
		
		@Override
		public void computeDerivatives (double t, double[] y, double[] yDot) {

			assert !Double.isNaN(t);
			
			Arrays.fill(yDot, 0d);
			// go through the matrix entries
			for ( Entry<Pair<Integer, Integer>, TimeDependentDouble> matrixEntry : this.rateMatrix.getEntrySet()) {
				// get the indices
				int srcIdx = matrixEntry.getKey().first();
				int dstIdx = matrixEntry.getKey().second();
				// and the value
				double value = matrixEntry.getValue().getValue(t);
				// and do something with it
				// it flows from src to dst, with strength value
				yDot[dstIdx] += value * y[srcIdx];
			}
			// should be good
		}

		@Override
		public int getDimension () {
			return this.rateMatrix.getDimension();
		}
		
		public String dumpRateMatrix (double t) {
			return CoalescentODEs.dumpSparseMatrix (this.rateMatrix, t);
		}
	}


	// reco is indexed by [firstEndAncientDeme, secondEndAncientDeme] dim=(2*numDemes, 2*numDemes)
	// indexing noReco via [endDeme] dim=(2*numDemes)
	public static SimplePair<double[], double[][]> computeR (TrunkProcess trunk, int epochIdx, int additionalDemeIdx, SimplePair<double[], double[][]> prevEpochR, double recombinationRate, boolean smcPrime) {

		Interval epoch = trunk.refinedDemography.epochList[epochIdx];
		double[][] pulseMigRates = trunk.refinedDemography.pulseMigrationMatrixList.get(epochIdx);
		double[][] migrationRates = trunk.refinedDemography.migrationMatrixList.get(epochIdx);
		double[] absorptionRates = trunk.getAbsorbRates(epochIdx);
		double[] expGrowthRates = trunk.refinedDemography.expGrowthRates.get(epochIdx);
		double[] initPopSizes = new double[trunk.refinedDemography.numAncientDemes(epochIdx)];
		for (int deme = 0; deme < trunk.refinedDemography.numAncientDemes(epochIdx); deme++) {
			initPopSizes[deme] = trunk.refinedDemography.getPopSize(deme, epochIdx, trunk.refinedDemography.epochList[epochIdx].endPoint);
		}
		int numDemes = pulseMigRates != null ? pulseMigRates.length : expGrowthRates.length;

		
		double[] startNoReco = new double[numDemes];
		double[][] startReco = new double[numDemes][numDemes];
		SimplePair<double[], double[][]> currEpochStartR = new SimplePair<double[], double[][]>(startNoReco, startReco);
		if (epochIdx == 0) {
			assert prevEpochR == null;
			startNoReco[additionalDemeIdx] = 1.0;
		} else {
			assert prevEpochR != null;
			
			// no reco
			for (int deme = 0; deme < numDemes; deme++) {
				for (int memberDemeIdx : trunk.refinedDemography.getMemberDemesIndices(deme, epochIdx)) {
					startNoReco[deme] += prevEpochR.first()[memberDemeIdx];
				}
			}
			
			// reco
			for (int leftDeme = 0; leftDeme < numDemes; leftDeme++) {
				for (int leftMemberDemeIdx : trunk.refinedDemography.getMemberDemesIndices(leftDeme, epochIdx)) {
					for (int rightDeme = 0; rightDeme < numDemes; rightDeme++) {
						for (int rightMemberDemeIdx : trunk.refinedDemography.getMemberDemesIndices(rightDeme, epochIdx)) {
							startReco[leftDeme][rightDeme] += prevEpochR.second()[leftMemberDemeIdx][rightMemberDemeIdx];							
						}
					}
				}
			}

		}
		
		
		if (pulseMigRates != null) {
			assert absorptionRates == null;
			assert migrationRates == null;
			assert expGrowthRates == null;
			
			double[] noReco = new double[2*numDemes];
			
			for (int i = 0; i < numDemes; i++) {
				assert pulseMigRates[i].length == numDemes;
				
				for (int j = 0; j < numDemes; j++) {
					noReco[j] = currEpochStartR.first()[i] * pulseMigRates[i][j];
				}
			}
			
			double[][] reco = new double[2*numDemes][2*numDemes];
			for (int leftStart = 0; leftStart < numDemes; leftStart++) {
				for (int rightStart = 0; rightStart < numDemes; rightStart++) {
					for (int leftEnd = 0; leftEnd < numDemes; leftEnd++) {
						for (int rightEnd = 0; rightEnd < numDemes; rightEnd++) {
							reco[leftEnd][rightEnd] += currEpochStartR.second()[leftStart][rightStart] * pulseMigRates[leftStart][leftEnd] * pulseMigRates[rightStart][rightEnd];
						}
					}
				}
			}
			
			return new SimplePair<double[], double[][]>(noReco, reco);
			
		}
		
		
		// just check the dimensions
		assert (migrationRates.length == absorptionRates.length);
		assert (RateMatrixTools.isSquare(migrationRates));
		assert (absorptionRates.length == expGrowthRates.length);

		// get some numbers
		// start time is zero
		double startTime = epoch.startPoint;
		double endTime = epoch.endPoint;

		// them indices
		RecoStatesDict recoStates = new RecoStatesDict (numDemes);
		
		// make y_start
		assert recoStates.indexToPair.size() == recoStates.numStates();
		double[] y_start = new double[recoStates.numStates()];
		for (int i = 0; i < y_start.length; i++) {
			Pair<Integer, Integer> currState = recoStates.indexToPair.get(i);
			int firstDeme = currState.first();
			int secondDeme = currState.second();
			if (firstDeme >= numDemes || secondDeme >= numDemes) continue;
			
			// if together
			if (secondDeme == -1) {
				y_start[i] = currEpochStartR.first()[firstDeme];
			} else {
				// if separate
				y_start[i] = currEpochStartR.second()[firstDeme][secondDeme];
			}
		}
		
		// build a rate matrix
		SparseMatrix<TimeDependentDouble> rateMatrix = new SparseMatrix<TimeDependentDouble> (recoStates.numStates());
		// fill it
		// loop over first deme
		for (int firstDeme=0; firstDeme<numDemes; firstDeme++) {
			
			// what happens if they are together
			for (int dstDeme=0; dstDeme<numDemes; dstDeme++) {
				// get them non-absorb indices
				int srcIdx = recoStates.getIndexTogether (firstDeme);
				int dstIdx = recoStates.getIndexTogether (dstDeme);
				// now see what to do
				if (firstDeme == dstDeme) {
					// we need another index
					int absorbIdx = recoStates.getIndexTogether (dstDeme + numDemes);
					// we have some negative diagonal, but also the absorption rate
					double constantPart = migrationRates[firstDeme][dstDeme] - recombinationRate;
					assert (constantPart <= 0d);
					// the given rates are backward in time growth rates, so since we want the inverse, we have to take the negative rate
					// no inverse of rates, cause is already an absorpton rate
					TimeDependentDouble diagonalValue = new BackwardExponentialDoublePlus (constantPart, - absorptionRates[dstDeme], endTime, - expGrowthRates[dstDeme]);
					rateMatrix.set (srcIdx, dstIdx, diagonalValue);
					// and the absorption rate
					TimeDependentDouble absRate = new BackwardExponentialDoublePlus (0, absorptionRates[dstDeme], endTime, - expGrowthRates[dstDeme]);
					rateMatrix.set (srcIdx, absorbIdx, absRate);
					// also possibly recombining
					TimeDependentDouble recoValue = new ConstantDouble(recombinationRate);
					int recoIdx = recoStates.getIndexApart (firstDeme, firstDeme);
					rateMatrix.set (srcIdx, recoIdx, recoValue);
				}
				else {
					// just normal transition
					rateMatrix.set (srcIdx, dstIdx, new ConstantDouble(migrationRates[firstDeme][dstDeme]));
				}
			}
			
			// what happens if they are apart (both free)
			// loop over second deme
			for (int secondDeme=0; secondDeme<numDemes; secondDeme++) {
				
				// both are free
				
				// where do we come from
				int srcIdx = recoStates.getIndexApart(firstDeme, secondDeme);
				
				// where do we go
				// possible back recombination
				if (smcPrime && firstDeme == secondDeme) {
					TimeDependentDouble backCoal = new BackwardExponentialDoublePlus (0d, 1.0 / initPopSizes[firstDeme], endTime, - expGrowthRates[firstDeme]);
					int dstIdx = recoStates.getIndexTogether (firstDeme);
					rateMatrix.set (srcIdx, dstIdx, backCoal);
				}
				
				// fill the diagonal
				double diagMig = migrationRates[firstDeme][firstDeme] + migrationRates[secondDeme][secondDeme];
				assert (diagMig <= 0d);
				TimeDependentDouble diagValue;
				if (smcPrime && firstDeme == secondDeme) {
					diagValue = new BackwardExponentialDoublePlus(diagMig, - 2*absorptionRates[firstDeme] - 1.0/initPopSizes[firstDeme], endTime, -expGrowthRates[firstDeme]);
				} else {					
					diagValue = new DoubleBackwardExponentialDoublePlus (diagMig, - absorptionRates[firstDeme], - absorptionRates[secondDeme], endTime, - expGrowthRates[firstDeme], - expGrowthRates[secondDeme]);
				}
				rateMatrix.set (srcIdx, srcIdx, diagValue);
				
				// off-diagonal
				// loop over target deme (can be at both first or second)
				for (int targetDeme=0; targetDeme<numDemes; targetDeme++) {
					// we already did the diagonal
					// WARNING: but there is still some absorption missing
					
					// where could the first deme go?
					if (targetDeme == firstDeme) {
						// could get absorbed
						// get the absorption index
						int absorbIdx = recoStates.getIndexApart (targetDeme + numDemes, secondDeme);
						TimeDependentDouble absRate = new BackwardExponentialDoublePlus (0d, absorptionRates[firstDeme], endTime, - expGrowthRates[firstDeme]);
						rateMatrix.set (srcIdx, absorbIdx, absRate);
					} else {
						// or migrate there
						int dstIdx = recoStates.getIndexApart (targetDeme, secondDeme);
						TimeDependentDouble firstRate = new ConstantDouble (migrationRates[firstDeme][targetDeme]);
						rateMatrix.set (srcIdx, dstIdx, firstRate);
					}
					
					// where could the second deme go?
					if (targetDeme == secondDeme) {
						// could get absorbed
						int absorbIdx = recoStates.getIndexApart (firstDeme, targetDeme + numDemes);
						TimeDependentDouble absRate = new BackwardExponentialDoublePlus (0d, absorptionRates[secondDeme], endTime, - expGrowthRates[secondDeme]);
						rateMatrix.set (srcIdx, absorbIdx, absRate);
					}
					else {
						// or migrate there
						int dstIdx = recoStates.getIndexApart (firstDeme, targetDeme);
						TimeDependentDouble secondRate = new ConstantDouble (migrationRates[secondDeme][targetDeme]);
						rateMatrix.set (srcIdx, dstIdx, secondRate);
					}
				}
				
				// one of them is fixed
				// which one?
				for (int i=0; i<2; i++) {
					boolean firstFixed = (i == 0);
					// get the source index
					int firstSrcDeme = (firstFixed ? firstDeme + numDemes : firstDeme);
					int secondSrcDeme = (firstFixed ?  secondDeme : secondDeme + numDemes);
					int oneFixedSrcIdx = recoStates.getIndexApart (firstSrcDeme, secondSrcDeme);
					int realSrcDeme = (firstFixed ? secondDeme : firstDeme);
					
					// now just marginal dynamics for the other one
					for (int realDstDeme=0; realDstDeme<numDemes; realDstDeme++) {
						// get non-absorb indices
						int firstDstDeme = (firstFixed ? firstSrcDeme : realDstDeme);
						int secondDstDeme = (firstFixed ?  realDstDeme : secondSrcDeme);
						int oneFixedDstIdx = recoStates.getIndexApart (firstDstDeme, secondDstDeme);
						
						// now see what to do
						if (realSrcDeme == realDstDeme) {
							// we need another index
							int firstAbsorb = (firstFixed ? firstSrcDeme : realDstDeme + numDemes);
							int secondAbsorb = (firstFixed ? realDstDeme +  numDemes : secondSrcDeme);
							int absorbIdx = recoStates.getIndexApart (firstAbsorb, secondAbsorb);
							// we have some negative diagonal, but also the absorption rate
							double constantPart = migrationRates[realSrcDeme][realDstDeme];
							assert (constantPart <= 0d);
							// the given rates are backward in time growth rates, so since we want the inverse, we have to take the negative rate
							// no inverse of rates, cause is already an absorption rate
							TimeDependentDouble diagonalValue = new BackwardExponentialDoublePlus (constantPart, - absorptionRates[realDstDeme], endTime, - expGrowthRates[realDstDeme]);
							assert (oneFixedSrcIdx == oneFixedDstIdx);
							rateMatrix.set (oneFixedSrcIdx, oneFixedDstIdx, diagonalValue);
							// and the absorption rate
							TimeDependentDouble absRate = new BackwardExponentialDoublePlus (0d, absorptionRates[realDstDeme], endTime, - expGrowthRates[realDstDeme]);
							rateMatrix.set (oneFixedSrcIdx, absorbIdx, absRate);					
						}
						else {
							// just normal transition
							rateMatrix.set (oneFixedSrcIdx, oneFixedDstIdx, new ConstantDouble(migrationRates[realSrcDeme][realDstDeme]));
						}
					}
				}
				// should be done for this first/second-pair
			}
		}
		assert (isRateMatrix(rateMatrix, (startTime + endTime)/2));
		
		// run ode
		// intialize return containers
		// [endDeme]
		double[] noReco = new double[2*numDemes];
		// [firstEndDeme][secondEndDeme]
		double[][] reco = new double[2*numDemes][2*numDemes];
		// and then run the ode to fill them
		// for each startDeme in this interval

//		MetaOptimization.synchronizedPrintln (Arrays.toString (y_start));
//		MetaOptimization.synchronizedPrintln (dumpSparseMatrix (rateMatrix, startTime));
		
		double initMass = sum(y_start);
		double[] y_end = solveODEHighamHall (y_start, rateMatrix, startTime, endTime);
		assert (Math.abs (sum(y_end) - initMass) < UberDemographyCore.ONELOCUS_EPSILON);
		
//		MetaOptimization.synchronizedPrintln (Arrays.toString (y_end));

		// and then fill the matrices
		// first no reco
		Arrays.fill(noReco, 0d);
		// with real values
		for (int endDeme=0; endDeme<2*numDemes; endDeme++) {
			noReco[endDeme] = y_end[recoStates.getIndexTogether(endDeme)];
		}
		
		// then reco
		for (int firstEnd=0; firstEnd<2*numDemes; firstEnd++){
			Arrays.fill(reco[firstEnd], 0d);
			for (int secondEnd=0; secondEnd<2*numDemes; secondEnd++) {
				reco[firstEnd][secondEnd] = y_end[recoStates.getIndexApart(firstEnd, secondEnd)];
			}
		}
		// that's it
		
		// check it
		// split reco no reco probability equal by 1/2
		TreeSet<Integer> realDemes = new TreeSet<Integer>();
		for (int i=0; i<absorptionRates.length; i++) {
			if (absorptionRates[i] > 0 + UberDemographyCore.ONELOCUS_EPSILON) {
				realDemes.add(i);
			}
		}
		// first no reco
		double tmpOneLocusEpsilon = UberDemographyCore.ONELOCUS_EPSILON * 100;
		if (epoch.endPoint != Double.POSITIVE_INFINITY) {
			for (int j=0; j<noReco.length; j++) {
				boolean nonZeroJ = (j<numDemes) || realDemes.contains(j - numDemes);
				assert ( (nonZeroJ && (noReco[j] > tmpOneLocusEpsilon)) || (noReco[j] < tmpOneLocusEpsilon));
			}
		}
		else {
			for (int j=numDemes; j<2*numDemes; j++) {
				assert ((realDemes.contains(j-numDemes) && (noReco[j] > tmpOneLocusEpsilon)) || (noReco[j] < tmpOneLocusEpsilon));
			}
		}
		// this stuff also has to sum to 1/2
		if (epoch.endPoint != Double.POSITIVE_INFINITY) {
			for (int j=0; j<reco.length; j++) {
				boolean nonZeroJ = (j<numDemes) || realDemes.contains(j - numDemes);
				for (int k=0; k<reco[j].length; k++) {
					boolean nonZeroK = (k<numDemes) || realDemes.contains(k - numDemes);
					assert ( (nonZeroJ && nonZeroK && (reco[j][k] > tmpOneLocusEpsilon)) || (reco[j][k] < tmpOneLocusEpsilon));
				}
			}
		}
		else {
			for (int j=0; j<2*numDemes; j++) {
				boolean nonZeroJ = realDemes.contains(j - numDemes);
				for (int k=0; k<2*numDemes; k++) {
					boolean nonZeroK = realDemes.contains(k - numDemes);
					assert ( (nonZeroJ && nonZeroK && (reco[j][k] > tmpOneLocusEpsilon)) || (reco[j][k] < tmpOneLocusEpsilon));
				}
			}
		}
		
		// return it in a pair
		return new SimplePair<double[], double[][]>(noReco, reco);
	}

	private static String dumpSparseMatrix (SparseMatrix<TimeDependentDouble> rateMatrix, double t) {
		// first make empty matrix
		double[][] dumpMatrix = new double[rateMatrix.getDimension()][rateMatrix.getDimension()];
		for (int i=0; i<dumpMatrix.length; i++) {
			for (int j=0; j<dumpMatrix[i].length; j++) {
				dumpMatrix[i][j] = 0d;
			}
		}
		
		// fill the entries
		for ( Entry<Pair<Integer, Integer>, TimeDependentDouble> matrixEntry : rateMatrix.getEntrySet()) {
			// get the indices
			int srcIdx = matrixEntry.getKey().first();
			int dstIdx = matrixEntry.getKey().second();
			// and the value
			dumpMatrix[srcIdx][dstIdx] = matrixEntry.getValue().getValue(t);
		}

		// and then write the dump matrix
		String toReturn = "{";
		for (int i=0; i<dumpMatrix.length; i++) {
			toReturn += "{";
			for (int j=0; j<dumpMatrix[i].length; j++) {
				toReturn += dumpMatrix[i][j];
				if (j == dumpMatrix[i].length-1) {
					toReturn += "}";
				}
				else {
					toReturn += ", ";
				}
			}
			if (i == dumpMatrix.length-1) {
				toReturn += "}";
			}
			else {
				toReturn += ",\n";
			}
		}
		return toReturn;
	}


	// [begin,nrMut,ending] dim=(numDeme, 2, 2*numDeme)
	public static double[][][] computeMutationEvents (Interval epoch, double mutRate, double[][] pulseMigRates, double[][] migrationRates, double[] absorptionRates, double[] expGrowthRates) {
		
		assert (MAX_NR_MUT_EVENTS == 1);
		
		if (pulseMigRates != null) {
			assert migrationRates == null;
			assert absorptionRates == null;
			assert expGrowthRates == null;
			
			int numDemes = pulseMigRates.length;
			
			double[][][] toReturn = new double[numDemes][2][2*numDemes];
			for (int i = 0; i < numDemes; i++) {
				assert pulseMigRates[i].length == numDemes;
				for (int j = 0; j < numDemes; j++) {
					toReturn[i][0][j] = pulseMigRates[i][j];
				}
			}
			return toReturn;
		}
		
		// just check the dimensions
		assert (migrationRates.length == absorptionRates.length);
		assert (RateMatrixTools.isSquare(migrationRates));
		assert (absorptionRates.length == expGrowthRates.length);

		// get some numbers
		// start time is zero
		double startTime = epoch.startPoint;
		double endTime = epoch.endPoint;
		
		// them indices
		int numDemes = expGrowthRates.length;
		MutStatesDict mutStates = new MutStatesDict (numDemes);
		// build a matrix
		SparseMatrix<TimeDependentDouble> rateMatrix = new SparseMatrix<TimeDependentDouble> (mutStates.numStates());
		// fill it
		// loop nr of mutations
		for (int nrMut=0; nrMut<MAX_NR_MUT_EVENTS+1; nrMut++) {
			// then loop over source deme
			for (int srcDeme=0; srcDeme<numDemes; srcDeme++) {
				// also over destination deme
				for (int dstDeme=0; dstDeme<numDemes; dstDeme++) {
					// get non-absorb indices
					int srcIdx = mutStates.getIndex (nrMut, srcDeme);
					int dstIdx = mutStates.getIndex (nrMut, dstDeme);
					// now see what to do
					if (srcDeme == dstDeme) {
						// we need another index
						int absorbIdx = mutStates.getIndex (nrMut, dstDeme + numDemes);
						
						// we have some negative diagonal, but also the absorption rate
						double constantPart = migrationRates[srcDeme][dstDeme];
						// subtract some mutation in case we have non yet
						if (nrMut == 0) {
							constantPart -= mutRate;
						}
						assert (constantPart <= 0d);
						// the given rates are backward in time growth rates, so since we want the inverse, we have to take the negative rate
						// no inverse of rates, cause is already an absorpton rate
						TimeDependentDouble diagonalValue = new BackwardExponentialDoublePlus (constantPart, - absorptionRates[dstDeme], endTime, - expGrowthRates[dstDeme]);
						assert (srcIdx == dstIdx);
						rateMatrix.set (srcIdx, dstIdx, diagonalValue);
						
						// and the absorption rate
						TimeDependentDouble absRate = new BackwardExponentialDoublePlus (0d, absorptionRates[dstDeme], endTime, - expGrowthRates[dstDeme]);
						rateMatrix.set (srcIdx, absorbIdx, absRate);
						
						// and also mutation, if we don't have one yet
						if (nrMut == 0) {
							int mutIdx = mutStates.getIndex (nrMut + 1, dstDeme);
							TimeDependentDouble timeMutRate = new ConstantDouble(mutRate);
							rateMatrix.set (srcIdx, mutIdx, timeMutRate);
						}
					}
					else {
						// just normal transition
						rateMatrix.set (srcIdx, dstIdx, new ConstantDouble(migrationRates[srcDeme][dstDeme]));
					}
				}
				// should be done for this mut/deme-pair
			}
		}
		// matrix might be fine
		assert (isRateMatrix(rateMatrix, (startTime + endTime)/2));

		
		// run ode
		// intialize return containers
		double[][][] mut = new double[numDemes][MAX_NR_MUT_EVENTS+1][2*numDemes];
		// and then run the ode to fill them
		// for each startDeme in this interval
		for (int startDeme=0; startDeme<numDemes; startDeme++) {
			// prepare stuff
			double[] y_start = new double[mutStates.numStates()];
			Arrays.fill(y_start, 0d);
			y_start[mutStates.getIndex (0, startDeme)] = 1d;
			// and run it
			double[] y_end = solveODEHighamHall (y_start, rateMatrix, startTime, endTime);
			assert (Math.abs (sum(y_end) - 1d) < UberDemographyCore.ONELOCUS_EPSILON);
			
			// and then fill matrices
			// then reco
			for (int nrMut=0; nrMut<mut[startDeme].length; nrMut++){
				Arrays.fill(mut[startDeme][nrMut], 0d);
				for (int endDeme=0; endDeme<2*numDemes; endDeme++) {
					mut[startDeme][nrMut][endDeme] = y_end[mutStates.getIndex(nrMut, endDeme)];
				}
			}
		}
		// maybe that's it

		// check it
		TreeSet<Integer> realDemes = new TreeSet<Integer>();
		for (int i=0; i<absorptionRates.length; i++) {
			if (absorptionRates[i] > 0 + UberDemographyCore.ONELOCUS_EPSILON) {
				realDemes.add(i);
			}
		}
		// real things
		double[] theSumz = new double[numDemes];
		double tmpOneLocusEpsilon = UberDemographyCore.ONELOCUS_EPSILON * 100;
		// split the probs evenly between zero or one mutation
		for (int i=0; i<mut.length; i++) {
			// no mut part
			if (epoch.endPoint != Double.POSITIVE_INFINITY) {
				for (int j=0; j<mut[i][0].length; j++) {
					boolean nonZeroJ = (j<numDemes) || realDemes.contains(j - numDemes);
					assert ( (nonZeroJ && (mut[i][0][j] > tmpOneLocusEpsilon)) || (mut[i][0][j] < tmpOneLocusEpsilon));
					theSumz[i] += mut[i][0][j];
				}
			}
			else {
				for (int j=0; j<mut[i][0].length; j++) {
					boolean nonZeroJ = realDemes.contains(j - numDemes);
					assert ( (nonZeroJ && (mut[i][0][j] > tmpOneLocusEpsilon)) || (mut[i][0][j] < tmpOneLocusEpsilon));
					theSumz[i] += mut[i][0][j];
				}
			}

			// mut part
			if (epoch.endPoint != Double.POSITIVE_INFINITY) {
				for (int j=0; j<mut[i][1].length; j++) {
					boolean nonZeroJ = (j<numDemes) || realDemes.contains(j - numDemes);
					assert ( (nonZeroJ && (mut[i][1][j] > tmpOneLocusEpsilon)) || (mut[i][1][j] < tmpOneLocusEpsilon));
					theSumz[i] += mut[i][1][j];
				}
			}
			else {
				for (int j=0; j<mut[i][1].length; j++) {
					boolean nonZeroJ = realDemes.contains(j - numDemes);
					assert ( (nonZeroJ && (mut[i][1][j] > tmpOneLocusEpsilon)) || (mut[i][1][j] < tmpOneLocusEpsilon));
					theSumz[i] += mut[i][1][j];
				}
			}
		}
		
		// return the mutation event matrix
		return mut;
	}
	
	public static class SparseMatrix<T> {
		private TreeMap<Pair<Integer, Integer>,T> entries;
		private int dimension;
		
		SparseMatrix (int dimension) {
			this.dimension = dimension;
			this.entries = new TreeMap<Pair<Integer,Integer>, T>();
		}
		
		public static double getAbsMaxRate (SparseMatrix<TimeDependentDouble> rateMatrix, double time) {
			
			// go though all values
			double maxRate = 0d;
			for (TimeDependentDouble value : rateMatrix.entries.values()) {
				// adjust rate maybe
				maxRate = Math.max (maxRate, value.getValue(time));
			}
			
			return maxRate;
		}

		public T get (int i, int j) {
			assert ((0 <= i) && (i < this.dimension));
			assert ((0 <= j) && (j < this.dimension));
			Pair<Integer,Integer> key = new Pair<Integer,Integer>(i,j);
			assert (this.entries.containsKey(key));
			return this.entries.get(key);
		}
		
		public void set (int i, int j, T value) {
			assert ((0 <= i) && (i < this.dimension));
			assert ((0 <= j) && (j < this.dimension));
			Pair<Integer,Integer> pairIdx = new Pair<Integer,Integer>(i,j);
			assert (!this.entries.containsKey(pairIdx) || this.entries.get(pairIdx).equals(value));
			this.entries.put (pairIdx, value);
		}
				
		public int getDimension () {
			return this.dimension;
		}
		
		public Set<Entry<Pair<Integer, Integer>, T>> getEntrySet () {
			return this.entries.entrySet();
		}
	}
	
	
	public static boolean isRateMatrix (SparseMatrix<TimeDependentDouble> m, double t) {
		// get some vector
		double[] sums = new double[m.getDimension()];
		Arrays.fill(sums, 0d);
		double maxAbsEntry = 0d;
		
		// and put all guys where they belong
		for (Entry<Pair<Integer, Integer>, TimeDependentDouble> e : m.getEntrySet()) {
			maxAbsEntry = Math.max(Math.abs(e.getValue().getValue(t)) , maxAbsEntry);
			sums[e.getKey().first()] += e.getValue().getValue(t);
		}
		
		// and see where we are at
		for (double value : sums) {
			if (Math.abs(value) > UberDemographyCore.ONELOCUS_EPSILON * maxAbsEntry) return false;
		}
		// we made it
		return true;
	}
	
	public static abstract class TimeDependentDouble {
		public abstract double getValue (double t);
		
		public static class ConstantDouble extends TimeDependentDouble {

			private double d;
			
			public ConstantDouble (double d) {
				this.d = d;
			}
			
			@Override
			public double getValue (double t) {
				// just ignore the time
				return this.d;
			}
		}
		
		public static class MinusSumDouble extends TimeDependentDouble {
			private final List<TimeDependentDouble> vars;
			
			public MinusSumDouble() {
				this.vars = new ArrayList<TimeDependentDouble>();
			}
			
			public void add(TimeDependentDouble d) {
				vars.add(d);
			}
			
			@Override
			public double getValue (double t) {
				double ret = 0d;
				for (TimeDependentDouble v : vars) ret -= v.getValue(t);
				return ret;
			}
		}
		
		public static class BackwardExponentialDoublePlus extends TimeDependentDouble {
			
			private final double b0;
			private final double d0;
			private final double rate;
			private final double t0;
			
			public BackwardExponentialDoublePlus (double b0, double d0, double t0, double rate) {
				this.b0 = b0;
				this.d0 = d0;
				this.t0 = t0;
				this.rate = rate;
			}

			@Override
			public double getValue (double t) {
				// outside is ok
				if (this.rate == 0d) {
					return this.b0 + this.d0;
				}
				else {
					return this.b0 + this.d0 * Math.exp((t0 - t)*this.rate);
				}
			}
		}
		
		public static class DoubleBackwardExponentialDoublePlus extends TimeDependentDouble {
			
			private final double b0;
			private final double d0;
			private final double rate0;
			private final double d1;
			private final double rate1;
			private final double t0;
			
			public DoubleBackwardExponentialDoublePlus (double b0, double d0, double d1, double t0, double rate0, double rate1) {
				this.b0 = b0;
				this.d0 = d0;
				this.d1 = d1;
				this.t0 = t0;
				this.rate0 = rate0;
				this.rate1 = rate1;
			}

			@Override
			public double getValue (double t) {
				// outside is ok
				// base
				double value = this.b0;
				// exponential zero
				if (this.rate0 == 0d) {
					value += this.d0;
				}
				else {
					value += this.d0 * Math.exp((t0 - t)*this.rate0);
				}
				// exponential one
				if (this.rate1 == 0d) {
					value += this.d1;
				}
				else {
					value += this.d1 * Math.exp((t0 - t)*this.rate1);
				}
				// give it away
				return value;
			}
		}

	}
	
	
	public static class MarginalStatesDict {
		// false is not yet absorbed, true means absorbed
		private final ArrayList<Pair<Boolean, Integer>> indexToPair;
		private final TreeMap<Pair<Boolean, Integer>, Integer> pairToIndex;
		
		public MarginalStatesDict (int numDemes) {
			// enumerate all pairs
			this.indexToPair = new ArrayList<Pair<Boolean,Integer>>();
			for (int abs=0; abs<2; abs++) {
				for (int deme=0; deme<numDemes; deme++) {
					// false is not yet absorbed, true means absorbed
					this.indexToPair.add (new Pair<Boolean, Integer> (abs != 0, deme));
				}
			}
			assert (this.indexToPair.size() == 2*numDemes);
			
			// then construct the reverse map
			this.pairToIndex = new TreeMap<Pair<Boolean,Integer>, Integer>();
			for (int i=0; i<this.indexToPair.size(); i++) {
				this.pairToIndex.put(this.indexToPair.get(i), i);
			}
		}
		
		public int getIndex (boolean absored, int deme) {
			return this.pairToIndex.get (new Pair<Boolean, Integer> (absored, deme));
		}
		
		public int numStates () {
			return this.indexToPair.size();
		}
	}
	
	public static class RecoStatesDict {
		// if the second one is minus one, this means they are still together
		private final ArrayList<Pair<Integer, Integer>> indexToPair;
		private final TreeMap<Pair<Integer, Integer>, Integer> pairToIndex;
		
		public RecoStatesDict (int numDemes) {
			// enumerate all pairs
			this.indexToPair = new ArrayList<Pair<Integer,Integer>>();
			for (int firstDeme=0; firstDeme<2*numDemes; firstDeme++) {
				for (int secondDeme=-1; secondDeme<2*numDemes; secondDeme++) {
					// -1 for second means together
					this.indexToPair.add (new Pair<Integer, Integer> (firstDeme, secondDeme));
				}
			}
			assert (this.indexToPair.size() == (2*numDemes)*(2*numDemes + 1));
			
			// then construct the reverse map
			this.pairToIndex = new TreeMap<Pair<Integer,Integer>, Integer>();
			for (int i=0; i<this.indexToPair.size(); i++) {
				this.pairToIndex.put(this.indexToPair.get(i), i);
			}
		}
		
		public int getIndexTogether (int deme) {
			return this.pairToIndex.get (new Pair<Integer, Integer> (deme, -1));
		}
		
		public int getIndexApart (int firstDeme, int secondDeme) {
			return this.pairToIndex.get (new Pair<Integer, Integer> (firstDeme, secondDeme));
		}

		public Pair<Integer, Integer> getState (int index) {
			return this.indexToPair.get(index);
		}
		
		public int numStates () {
			return this.indexToPair.size();
		}
	}

	public static class MutStatesDict {
		// if the second one is minus one, this means they are still together
		private final ArrayList<Pair<Integer, Integer>> indexToPair;
		private final TreeMap<Pair<Integer, Integer>, Integer> pairToIndex;
		
		public MutStatesDict (int numDemes) {
			// enumerate all pairs
			this.indexToPair = new ArrayList<Pair<Integer,Integer>>();
			for (int nrMut=0; nrMut<MAX_NR_MUT_EVENTS+1; nrMut++) {
				for (int deme=0; deme<2*numDemes; deme++) {
					// first is number mutation events, second is deme
					this.indexToPair.add (new Pair<Integer, Integer> (nrMut, deme));
				}
			}
			assert (this.indexToPair.size() == 2*(2*numDemes));
			
			// then construct the reverse map
			this.pairToIndex = new TreeMap<Pair<Integer,Integer>, Integer>();
			for (int i=0; i<this.indexToPair.size(); i++) {
				this.pairToIndex.put(this.indexToPair.get(i), i);
			}
		}
		
		public int getIndex (int nrMut, int deme) {
			assert (0 <= nrMut);
			assert (nrMut <= MAX_NR_MUT_EVENTS);
			return this.pairToIndex.get (new Pair<Integer, Integer> (nrMut, deme));
		}

		public Pair<Integer, Integer> getState (int index) {
			return this.indexToPair.get(index);
		}
		
		public int numStates () {
			return this.indexToPair.size();
		}
	}

	
	public static double sum (double[] array) {
		double sum=0d;
		for (double value : array) {
			sum += value;
		}
		return sum;
	}
	

}
