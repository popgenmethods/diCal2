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
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import Jama.Matrix;
import edu.berkeley.diCal2.csd.DemoState.DemoStateCollection;
import edu.berkeley.diCal2.csd.auxiliary.MyEigenDecomposition;
import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.utility.AncestralProcessProbs;
import edu.berkeley.diCal2.utility.RateMatrixTools;
import edu.berkeley.diCal2.utility.RealPartition.Interval;

public abstract class TrunkProcess {

	public enum CakeStyle {
		BEGINNING, MIDDLE, END, AVERAGE
	}
	
	public static CakeStyle getCakeStyle(String cakeStyleName) {
		if (cakeStyleName.equals("beginning")) {
			return CakeStyle.BEGINNING;
		}
		else if (cakeStyleName.equals("middle")) {
			return CakeStyle.MIDDLE;
		}
		else if (cakeStyleName.equals("end")) {
			return CakeStyle.END;
		} else if (cakeStyleName.equals("average")) {
			return CakeStyle.AVERAGE;
		}
		return null;
	}
	
	protected final Demography refinedDemography;
	protected final int[] sampleSizes;
	
	public TrunkProcess(int [] sampleSizes, Demography refinedDemography) {
		this.sampleSizes = Arrays.copyOf(sampleSizes, sampleSizes.length);
		this.refinedDemography = refinedDemography;
	}
	
	public abstract double fractionAncientDemeToPresentDeme(int refinedInterval, int presentDeme, int ancientDeme);
	public abstract double getAbsorbRate(int refinedInterval, int ancientDeme);
	
	public boolean isNonAbsorbing (int epochIdx, int ancientDemeIdx, double EPSILON){
		if (this.refinedDemography.isPulse(epochIdx)) {
			return true;
		}
		else {
			// times popSIze, becuase they can get big
			return this.getAbsorbRate(epochIdx, ancientDemeIdx) * this.refinedDemography.popSizes.get(epochIdx)[ancientDemeIdx] < EPSILON;
		}		
	}

	
	public double[] getAbsorbRates(int refinedInterval) {
		
		if (this.refinedDemography.isPulse(refinedInterval)) return null;
		
		int numDemes = this.refinedDemography.treePartitionList.get(refinedInterval).size();
		
		double[] toReturn = new double[numDemes];
		
		for (int ancientDeme = 0; ancientDeme < toReturn.length; ancientDeme++) {
			toReturn[ancientDeme] = this.getAbsorbRate(refinedInterval, ancientDeme);
		}
		
		return toReturn;
	}
	
	// GETTERS
	
	public Demography getRefinedDemography() {
		return refinedDemography;
	}

	public int[] getSampleSizes(){
		return this.sampleSizes;
	}
	
	// NESTED FACTORY CLASS
	public static class TrunkProcessFactory {

		public final TrunkStyle trunkStyle;
		public final CakeStyle cakeStyle;
				
		public TrunkProcessFactory(TrunkStyle trunkStyle, CakeStyle cakeStyle) {
			this.trunkStyle = trunkStyle;
			this.cakeStyle = cakeStyle;
		}

		public TrunkProcessFactory(String trunkStyleName, CakeStyle cakeStyle) {
			this (TrunkProcessFactory.getTrunkStyle(trunkStyleName), cakeStyle);
		}
		
		public synchronized TrunkProcess getTrunk (int[] sampleSizes, DemoStateCollection demoStates) {
			return this.getTrunk(sampleSizes, demoStates != null ? demoStates.refinedDemography : null);
		}
		
		public boolean halfMigrationRate () {
			if (halfMigrationTrunks.contains(this.trunkStyle)) {
				return true;
			}
			return false;
		}
		private TrunkProcess getTrunk (int[] sampleSizes, Demography refinedDemography) {
			int totalSampleSize= 0;
			for (int n : sampleSizes){
				totalSampleSize += n;
			}
			if (totalSampleSize == 0) return null;
			
			if (this.trunkStyle == TrunkStyle.SimpleTrunk) {
				return new SimpleTrunk (sampleSizes, refinedDemography);
			}
			if (this.trunkStyle == TrunkStyle.OldCakeTrunk) {
				return new ConstPopSizeCake (sampleSizes, refinedDemography, this.cakeStyle);
			}
			if (this.trunkStyle == TrunkStyle.MeanCakeTrunk) {
				return new MeanPopSizeCake(sampleSizes, refinedDemography, this.cakeStyle);
			}
			if (this.trunkStyle == TrunkStyle.MultiCakeTrunk) {
				return new MultiCakeTrunk(sampleSizes, refinedDemography, this.cakeStyle);
			}
			if (this.trunkStyle == TrunkStyle.MultiCakeUpdatingTrunk) {
				return new MultiCakeUpdatingTrunk(sampleSizes, refinedDemography, this.cakeStyle);
			}
			if (this.trunkStyle == TrunkStyle.MigratingMultiCakeTrunk) {
				return new SimpleMigratingTrunk(sampleSizes, refinedDemography, this.cakeStyle, new MultiCakeTrunk(sampleSizes, refinedDemography, this.cakeStyle));
			}
			if (this.trunkStyle == TrunkStyle.MigratingEthan) {
				return new EthanTrunk(sampleSizes, refinedDemography, this.cakeStyle);
			}
			if (this.trunkStyle == TrunkStyle.MigMultiCakeUpdatingTrunk) {
				return new SimpleMigratingTrunk(sampleSizes, refinedDemography, this.cakeStyle, new MultiCakeUpdatingTrunk(sampleSizes, refinedDemography, this.cakeStyle));
			}
			if (this.trunkStyle == TrunkStyle.RecursiveTrunk){
				// should be using the RecursiveTrunkFactory subclass
				assert (false);
			}
			if (this.trunkStyle == TrunkStyle.ExactCake) {
				assert (this.cakeStyle == null);
				return new SinglePopExactTrunk(sampleSizes, refinedDemography);
			}
			
			
			// should not happen
			assert (false);
			return null;
		}
		
		public enum TrunkStyle { SimpleTrunk, OldCakeTrunk, MeanCakeTrunk, MultiCakeTrunk, MultiCakeUpdatingTrunk, MigratingMultiCakeTrunk, MigMultiCakeUpdatingTrunk, RecursiveTrunk, ExactCake, MigratingEthan };
		private static Set<TrunkStyle> halfMigrationTrunks = new HashSet<TrunkStyle>(Arrays.asList(new TrunkStyle[] {TrunkStyle.MigratingMultiCakeTrunk, TrunkStyle.MigMultiCakeUpdatingTrunk, TrunkStyle.RecursiveTrunk, TrunkStyle.MigratingEthan}));
		
		public static TrunkStyle getTrunkStyle(String trunkStyleName) {
			if (trunkStyleName.equals("simple")) {
				return TrunkStyle.SimpleTrunk;
			} else if (trunkStyleName.equals("oldCake")) {
				return TrunkStyle.OldCakeTrunk;
			} else if (trunkStyleName.equals("meanCake")) {
				return TrunkStyle.MeanCakeTrunk;
			} else if (trunkStyleName.equals("multiCake")) {
				return TrunkStyle.MultiCakeTrunk;
			} else if (trunkStyleName.equals("multiCakeUpdating")) {
				return TrunkStyle.MultiCakeUpdatingTrunk;
			} else if (trunkStyleName.equals("migratingMultiCake")) {
				return TrunkStyle.MigratingMultiCakeTrunk;
			} else if (trunkStyleName.equals("migMultiCakeUpdating")) {
				return TrunkStyle.MigMultiCakeUpdatingTrunk;
			} else if (trunkStyleName.equals("recursive")) {
				return TrunkStyle.RecursiveTrunk;
			} else if (trunkStyleName.equals ("exactCake")) {
				return TrunkStyle.ExactCake;
			} else if (trunkStyleName.equals ("migratingEthan")) {
				return TrunkStyle.MigratingEthan;
			}
			
			return null;
		}
	}

	
	// NESTED SUBCLASSES
	// each lineage goes back in time forever and never moves or changes mass
	public static class SimpleTrunk extends TrunkProcess {

		public SimpleTrunk (int[] sampleSizes, Demography refinedDemography) {
			super (sampleSizes, refinedDemography);
		}		
		
		@Override
		public double fractionAncientDemeToPresentDeme (int refinedInterval, int presentDeme, int ancientDeme) {
			if (this.refinedDemography.treePartitionList.get(refinedInterval).get(ancientDeme).contains(presentDeme)) {
				
				int presentDemeSize = sampleSizes[presentDeme];
				int ancientDemeSize = this.getAncientDemeSize(refinedInterval, ancientDeme);
				
				if (ancientDemeSize == 0) {
					assert (presentDemeSize == 0);
					return 0d;
				}
				
				return (double) presentDemeSize / (double) ancientDemeSize;
			}
			
			return 0d;
		}

		@Override
		public double getAbsorbRate(int refinedInterval, int ancientDeme) {
			return this.getAncientDemeSize(refinedInterval, ancientDeme) / this.refinedDemography.popSizes.get(refinedInterval)[ancientDeme];
		}
		
		
		protected int getAncientDemeSize(int refinedInterval, int ancientDeme) {
			int ancientDemeSize = 0;
			for (int deme : this.refinedDemography.treePartitionList.get(refinedInterval).get(ancientDeme)) {
				ancientDemeSize += this.sampleSizes[deme];
			}
			return ancientDemeSize;
		}
	}
	
	// each lineage goes back in time forever and never migrates.
	// all lineages lose mass at the same rate, no matter what deme they started in
	// the rate at which mass is loss is given by the "effective population size" of the demography at that time
	public static abstract class SimpleCakeTrunk extends SimpleTrunk {
		final double[] cakeFactors;
		
		public SimpleCakeTrunk(int[] sampleSizes, Demography refinedDemo, CakeStyle cakeStyle) {
			super (sampleSizes, refinedDemo);
		
			assert (refinedDemo == this.getRefinedDemography());
			
			int totalSampleSize= 0;
			for (int n : sampleSizes){
				totalSampleSize += n;
			}
			
			
			if (refinedDemo == null) {
				assert(totalSampleSize == 0);
				cakeFactors = null;
				return;
			}
						
			double[] popSizes = new double[refinedDemo.epochList.length];
			for (int i = 0; i < popSizes.length; i++) {
				popSizes[i] = effectivePopSize(refinedDemo.popSizes.get(i), refinedDemo.migrationMatrixList.get(i));
			}
			
			cakeFactors = TrunkProcess.getCakeFactors (totalSampleSize, this.refinedDemography.epochList, popSizes, cakeStyle);
		}
		
		@Override
		public double getAbsorbRate(int refinedInterval, int ancientDeme) {
			return cakeFactors[refinedInterval] * this.getAncientDemeSize(refinedInterval, ancientDeme) / this.refinedDemography.popSizes.get(refinedInterval)[ancientDeme];
		}
		
		protected abstract double effectivePopSize (double[] popSizes, double[][] migMatrix);
	}
	
	// Effective pop size = 1
	public static class ConstPopSizeCake extends SimpleCakeTrunk {		

		public ConstPopSizeCake(int[] sampleSizes, Demography refinedDemo, CakeStyle cakeStyle) {
			super(sampleSizes, refinedDemo, cakeStyle);
		}
		
		@Override
		protected double effectivePopSize (double[] popSizes, double[][] migMatrix) {
			return 1d;
		}
	}
	
	// Effective pop size = mean of pop sizes
	public static class MeanPopSizeCake extends SimpleCakeTrunk {		

		public MeanPopSizeCake(int[] sampleSizes, Demography refinedDemo, CakeStyle cakeStyle) {
			super(sampleSizes, refinedDemo, cakeStyle);
		}
		
		@Override
		protected double effectivePopSize (double[] popSizes, double[][] migMatrix) {
			
			if (popSizes == null) return 1d;
			
			double sum = 0d;
			for (int i = 0; i < popSizes.length; i++) {
				sum += popSizes[i];
			}
			return sum / popSizes.length;
		}
	}	
	
	
	// Each lineage never migrates
	// Each present deme loses absorption mass according to an independent cake
	public static class MultiCakeTrunk extends SimpleTrunk {

		// [presentDeme][refinedInterval]
		final double[][] cakeFactors;
		
		public MultiCakeTrunk(int[] sampleSizes, Demography refinedDemo, CakeStyle cakeStyle) {
		
			super (sampleSizes, refinedDemo);
			int totalSampleSize= 0;
			for (int n : sampleSizes){
				totalSampleSize += n;
			}
			if (this.refinedDemography == null) {
				assert (totalSampleSize == 0);
				cakeFactors = null;
				return;
			}
			
			
			// [presentDeme][refinedInterval]
			double[][] popSizes = new double[sampleSizes.length][this.refinedDemography.epochList.length];
			
			// fill popsizes
			for (int epoch = 0; epoch < this.refinedDemography.epochList.length; epoch++) {
				for (int ancientDeme = 0; ancientDeme < this.refinedDemography.treePartitionList.get(epoch).size(); ancientDeme++) {
					for (int presentDeme : this.refinedDemography.treePartitionList.get(epoch).get(ancientDeme)) {
						popSizes[presentDeme][epoch] = this.refinedDemography.isPulse(epoch) ? 1d : this.refinedDemography.popSizes.get(epoch)[ancientDeme];
					}
				}
			}
			
			this.cakeFactors = new double[sampleSizes.length][];
			for (int presentDeme = 0; presentDeme < sampleSizes.length; presentDeme++) {
				cakeFactors[presentDeme] = TrunkProcess.getCakeFactors (sampleSizes[presentDeme], this.refinedDemography.epochList, popSizes[presentDeme], cakeStyle);
			}
		}

		@Override
		public double getAbsorbRate(int refinedInterval, int ancientDeme) {
			
			double absorbRate = 0d;
			for (int presentDeme : this.refinedDemography.treePartitionList.get(refinedInterval).get(ancientDeme)) {
				absorbRate += this.cakeFactors[presentDeme][refinedInterval] * sampleSizes[presentDeme];
			}
			absorbRate /= this.refinedDemography.popSizes.get(refinedInterval)[ancientDeme];
			
			return absorbRate;
		}

	}
	
	public static class MultiCakeUpdatingTrunk extends SimpleTrunk {

		// [refinedInterval][ancientDeme]
		final private double[][] trunkSize;

		public MultiCakeUpdatingTrunk (int[] sampleSizes, Demography refinedDemo, CakeStyle cakeStyle) {

			super (sampleSizes, refinedDemo);
			
			int totalSampleSize= 0;
			for (int n : sampleSizes){
				totalSampleSize += n;
			}
			if (this.refinedDemography == null) {
				// should only be null for an empty trunk
				assert (totalSampleSize == 0);
				this.trunkSize = null;
				return;
			}
			
			// [refinedInterval][ancientDeme]
			this.trunkSize = new double[this.refinedDemography.epochList.length][];
			
			// get some representative times for each epoch. effectively the points where we want the trunk size
			double[] trunkTimes = null;
			// only fill it for the right ones, where it makes sense
			// average has no representative time, which kind of seems to make sense
			if (hasRepresentativeTimes(cakeStyle)) {
				trunkTimes = getRepresentativeTimes (this.refinedDemography.epochList, cakeStyle);
			}
			
			// some variables
			double[] lastEndingSizes = null;
			// calculate them iteratively
			for (int epochIdx = 0; epochIdx < this.trunkSize.length; epochIdx++) {
				
				// get the sizes at the start of this epoch
				double[] startingSizes = new double[this.refinedDemography.treePartitionList.get(epochIdx).size()];
				// by summing previous stuff in suitable ways
				for (int ancientDeme=0; ancientDeme <startingSizes.length; ancientDeme++) {
					// inside or beginning?
					if (epochIdx > 0) {
						// sum over where the previous guys were
						for (int memberDeme : refinedDemography.getMemberDemesIndices(ancientDeme, epochIdx)) {
							startingSizes[ancientDeme] += lastEndingSizes[memberDeme];
						}
					}
					else {
						// in first epoch this is just the initial guys
						assert (epochIdx == 0);
						assert (startingSizes.length == sampleSizes.length);
						startingSizes[ancientDeme] = sampleSizes[ancientDeme];
					}
				}
				
				// get the time of observation in this interval
				Double thisTime = null;
				if (trunkTimes != null) {
					thisTime = trunkTimes[epochIdx];
				}
				
				// now go through every ancient deme, and do some shrinking trunk there
				// pulse or not?
				if (this.refinedDemography.isPulse(epochIdx)) {
					// just copy
					lastEndingSizes = startingSizes;
					// and the cakeFactors just the same
					// no cakeFactors should be necessary
				}
				else {
					// no pulse, so this is the real deal
					// make some room for the stuff
					this.trunkSize[epochIdx] = new double[startingSizes.length];
					// also remember in each deme where you stopped
					lastEndingSizes = new double[startingSizes.length];
					
					// loop
					for (int ancientDeme=0; ancientDeme<startingSizes.length; ancientDeme++) {
						
						// find the right t_0 for the expected number of lineages to match with the starting size
						// get a close starting point
						int newInitialN = (int)Math.ceil(startingSizes[ancientDeme]);
						double[] preComputedFractions = preComputeExpectedNumberFraction (newInitialN);
						double t_0 = findTimeZero (newInitialN, startingSizes[ancientDeme], preComputedFractions);
						
						// we need some acceleration, cause trunk shrinks faster in small populations
						double accelerationFactor = 1/this.refinedDemography.popSizes.get(epochIdx)[ancientDeme];
						
						// now that we know where to start, see how much we have left at the end of the interval
						double endTime = t_0 + accelerationFactor * (this.refinedDemography.epochList[epochIdx].endPoint - this.refinedDemography.epochList[epochIdx].startPoint);
						
						// and store the end, cause we need to know it for the next time
						lastEndingSizes[ancientDeme] = getExpectedNumberLineages(newInitialN, endTime, preComputedFractions);
						
						// and we might also want an intermediate time
						if (cakeStyle == CakeStyle.BEGINNING) {
							// should be beginning
							assert (thisTime == this.refinedDemography.epochList[epochIdx].startPoint);
							// beginning trunk
							this.trunkSize[epochIdx][ancientDeme] = startingSizes[ancientDeme];
						}
						else if (cakeStyle == CakeStyle.END) {
							// should be end
							assert (thisTime == this.refinedDemography.epochList[epochIdx].endPoint);
							// ending trunk
							this.trunkSize[epochIdx][ancientDeme] = lastEndingSizes[ancientDeme];
						}
						else if (cakeStyle == CakeStyle.MIDDLE) {
							assert (thisTime > this.refinedDemography.epochList[epochIdx].startPoint);
							assert (thisTime < this.refinedDemography.epochList[epochIdx].endPoint);
							// get the time elapsing after t_0, accoutnig for the accelerator
							double adjustedTime = t_0 + accelerationFactor * (thisTime - this.refinedDemography.epochList[epochIdx].startPoint);
							
							// and then get the right factor
							this.trunkSize[epochIdx][ancientDeme] = getExpectedNumberLineages(newInitialN, adjustedTime, preComputedFractions);
						}
						else if (cakeStyle == CakeStyle.AVERAGE) {
							assert (thisTime == null);
							// and do some averaging
							this.trunkSize[epochIdx][ancientDeme] = (startingSizes[ancientDeme] + lastEndingSizes[ancientDeme])/2d;
						}
						else {
							// no cakestyles left
							assert (false);
						}
					}
				}
				
				// done with this epoch				
			}
		}

		private static double findTimeZero(int numIndividuals, double n_0, double[] preComputedFractions) {

			// get some start
			// this has to be below
			double lowerBound = 0d;
			assert (getExpectedNumberLineages(numIndividuals, lowerBound, preComputedFractions) >= n_0);
			// this has to be above
			double upperBound = Math.pow(10, 3);
			assert (getExpectedNumberLineages(numIndividuals, upperBound, preComputedFractions) <= n_0);
			
			// then update the bounds until you are close enough
			while (true) {
				// compute the pivot
				double pivot = (upperBound + lowerBound)/2;
				// compute value for pivot
				double n_p = getExpectedNumberLineages (numIndividuals, pivot, preComputedFractions);
				
				// now update the right bound
				if (n_p < n_0) {
					// pivot value is below n_0, so pivot is the new upper bound
					upperBound = pivot;
				}
				else {
					// pivot value is above or equal to n_0, so pivot is new lower bound
					lowerBound = pivot;
				}
				
				// if relative error below bound, we are done
				// WARING: a very small upper bound might inflate the relative error, so we need to check for that too
				if (((upperBound - lowerBound)/upperBound < 0.0001) || (upperBound < UberDemographyCore.ONELOCUS_EPSILON)) break;
			}
			
			// return the mean of lower and upper
			return (lowerBound + upperBound)/2d;
		}

		@Override
		public double getAbsorbRate(int refinedInterval, int ancientDeme) {
			assert (!this.refinedDemography.isPulse(refinedInterval));
			// absorption rate is trunkSize / popSize
			return this.trunkSize[refinedInterval][ancientDeme] / this.refinedDemography.popSizes.get(refinedInterval)[ancientDeme];
		}
	}

	public static class EthanTrunk extends MigratingTrunk {

		// [refinedInterval][ancientDeme]
		final private double[][] trunkSize;

		
		public EthanTrunk (int[] sampleSizes, Demography refinedDemo, CakeStyle cakeStyle) {
			// super
			super (sampleSizes, refinedDemo, cakeStyle);
			
			// don't want this one anymore
			assert (cakeStyle != CakeStyle.MIDDLE);
			
			int totalSampleSize= 0;
			for (int n : sampleSizes){
				totalSampleSize += n;
			}
			if (this.refinedDemography == null) {
				// should only be null for an empty trunk
				assert (totalSampleSize == 0);
				this.trunkSize = null;
				return;
			}
			
			// [refinedInterval][ancientDeme]
			this.trunkSize = new double[this.refinedDemography.epochList.length][];
			
			// get some representative times for each epoch. effectively the points where we want the trunk size
			double[] trunkTimes = null;
			// only fill it for the right cakes, where it makes sense
			// average has no representative time, which kind of seems to make sense
			if (hasRepresentativeTimes(cakeStyle)) {
				trunkTimes = getRepresentativeTimes (this.refinedDemography.epochList, cakeStyle);
			}
			
			// some variables
			double[] lastEndingSizes = null;
			// calculate them iteratively
			for (int epochIdx = 0; epochIdx < this.trunkSize.length; epochIdx++) {
				
				// get the sizes at the start of this epoch
				double[] startingSizes = new double[this.refinedDemography.treePartitionList.get(epochIdx).size()];
				// by summing previous stuff in suitable ways
				for (int ancientDeme=0; ancientDeme <startingSizes.length; ancientDeme++) {
					// inside or beginning?
					if (epochIdx > 0) {
						// sum over where the previous guys were
						for (int memberDeme : refinedDemography.getMemberDemesIndices(ancientDeme, epochIdx)) {
							startingSizes[ancientDeme] += lastEndingSizes[memberDeme];
						}
					}
					else {
						// in first epoch this is just the initial guys
						assert (epochIdx == 0);
						assert (startingSizes.length == sampleSizes.length);
						startingSizes[ancientDeme] = sampleSizes[ancientDeme];
					}
				}
				
				// get the time of observation in this interval
				Double thisTime = null;
				if (trunkTimes != null) {
					thisTime = trunkTimes[epochIdx];
				}
				
				
				// now go through every ancient deme, and do some shrinking trunk there
				// pulse or not?
				if (this.refinedDemography.isPulse(epochIdx)) {
					
					// actually, it is like making a one step transition using the pulse matrix
					double[][] pulseMatrix = this.refinedDemography.pulseMigrationMatrixList.get(epochIdx);
					assert (RateMatrixTools.isProperStochaticMatrix (pulseMatrix, UberDemographyCore.ONELOCUS_EPSILON));
					lastEndingSizes = new double[startingSizes.length];
					for (int j=0; j<lastEndingSizes.length; j++) {
						lastEndingSizes[j] = 0d;
						for (int i=0; i<startingSizes.length; i++) {
							lastEndingSizes[j] +=  startingSizes[i] * pulseMatrix[i][j];
						}
					}
					// should be good
				}
				else {
					
					// evolve the ODE starting from startingSizes for given amount of time
					double startTime = this.refinedDemography.epochList[epochIdx].startPoint;
					// make a suitable end time
					double endTime = this.refinedDemography.epochList[epochIdx].endPoint;
					assert (startTime >= 0d);
					assert (startTime < endTime);
					assert (endTime <= Double.POSITIVE_INFINITY);
					
					// final preparation
					lastEndingSizes = new double[startingSizes.length];
					// get them exponential rates if necessary
					double[] expRates = null;
					if (this.refinedDemography.expGrowthRates != null) {
						expRates = this.refinedDemography.expGrowthRates.get(epochIdx);
						// no exponential growth in last interval
						if (epochIdx >= this.trunkSize.length-1) {
							for (double rate : expRates) {
								assert (rate == 0);
							}
						}
					}
					
					// then finally do the ODE
					// I guess we need to fill a matrix of sorts (this one is not really sparse, but it's small, so that's ok)
					EthanODE ode = new EthanODE (this.refinedDemography.migrationMatrixList.get(epochIdx), this.refinedDemography.popSizes.get(epochIdx), expRates, endTime);
					
					// then run the ODE
					lastEndingSizes = CoalescentODEs.solveODEHighamHall(startingSizes, ode, startTime, endTime);
					
					// be safe
					assert (lastEndingSizes.length == startingSizes.length);
					
					// do something for the different cake styles
					if (cakeStyle == CakeStyle.BEGINNING) {
						// trunk is at beginning
						this.trunkSize[epochIdx] = Arrays.copyOf (startingSizes, startingSizes.length);
					}
					else if (cakeStyle == CakeStyle.END) {
						// sample at end
						this.trunkSize[epochIdx] = Arrays.copyOf (lastEndingSizes, lastEndingSizes.length);
					}
					else if (cakeStyle == CakeStyle.AVERAGE) {
						assert (thisTime == null);
						// make some room for the stuff
						this.trunkSize[epochIdx] = new double[startingSizes.length];
						for (int i=0; i<this.trunkSize[epochIdx].length; i++) {
							// and do some averaging
							this.trunkSize[epochIdx][i] = (startingSizes[i] + lastEndingSizes[i])/2d;
						}
					}
					else {
						// no more middle
						assert (false);
					}
					
					// should be done
					
					// trunk sizes are filled
				}
				
				// done with this epoch
			}
		}

		private static class EthanODE implements FirstOrderDifferentialEquations {

			final private double[][] migrationRates;
			final private double[] popSizes;
			final private double[] expGrowthRates;
			final private double endTime;
			
			public EthanODE (double[][] migrationRates, double[] popSizes, double[] expGrowthRates, double endTime) {
				// be safe
				assert (popSizes.length == migrationRates.length);
				assert ((expGrowthRates == null) || (expGrowthRates.length == popSizes.length));
				for (int i=0; i<migrationRates.length; i++) {
					assert (popSizes.length == migrationRates[i].length);
					assert (RateMatrixTools.isProperRateMatrix(migrationRates, UberDemographyCore.ONELOCUS_EPSILON));
				}
				
				// remember values
				this.migrationRates = migrationRates;
				this.popSizes = popSizes;
				this.expGrowthRates = expGrowthRates;
				this.endTime = endTime;
			}
			

			@Override
			public void computeDerivatives (double t, double[] y, double[] yDot) {

				// actual derivative for the ode 
				for (int l=0; l<yDot.length; l++) {
					// first the migration part (negative sum on the diagonal, so everything fine)
					double migTerm = 0d;
					for (int i=0; i<y.length; i++) {
						migTerm += this.migrationRates[i][l] * y[i];
					}
					yDot[l] = migTerm;
							
					// then maybe coalescence part
					// I think the variance would make this correct in real life
					if (y[l] > 1) {
						// add it (the negative one)
						if (this.expGrowthRates == null) {
							yDot[l] += - y[l]*(y[l]-1)/(2*this.popSizes[l]);
						}
						else {
							double growingSize = this.popSizes[l];
							// no growth in last interval should be checked outside
							if (this.expGrowthRates[l] != 0d) {
								growingSize *= Math.exp(this.expGrowthRates[l] * (this.endTime - t));
							}
							yDot[l] += - y[l]*(y[l]-1)/(2*growingSize);
							
						}
					}
				}
				
				// should be good
			}

			@Override
			public int getDimension () {
				return popSizes.length;
			}
		}
		
		@Override
		public double getAbsorbRate(int refinedInterval, int ancientDeme) {
			assert (!this.refinedDemography.isPulse(refinedInterval));
			// absorption rate is trunkSize / popSize
			return this.trunkSize[refinedInterval][ancientDeme] / this.refinedDemography.popSizes.get(refinedInterval)[ancientDeme];
		}
	}

	private static boolean hasRepresentativeTimes(CakeStyle cakeStyle) {
		return (new TreeSet<CakeStyle>(Arrays.asList(new CakeStyle[]{CakeStyle.BEGINNING, CakeStyle.MIDDLE, CakeStyle.END}))).contains(cakeStyle);
	}
	
	private static double[] getRepresentativeTimes (Interval[] intervals, CakeStyle cakeStyle) {
		
		assert (hasRepresentativeTimes(cakeStyle));
		
		// get a representative point in each interval of constant size
		double[] toReturn = new double[intervals.length];
		for (int i = 0; i < intervals.length; i++) {
			toReturn[i] = getRepresentativePoint(intervals[i], cakeStyle);
		}
		return toReturn;
	}
	
	private static Double getRepresentativePoint (Interval i, CakeStyle cakeStyle) {
		assert (i.startPoint <= i.endPoint);
		if (cakeStyle == CakeStyle.BEGINNING) {
			return i.startPoint;
		}
		else if (cakeStyle == CakeStyle.END) {
			return i.endPoint;
		}
		else if (cakeStyle == CakeStyle.MIDDLE) {
			// get the 50 % quantile of the trunctated exponential with given rate
			// maybe later
			// for now just the middle
			if (i.endPoint != Double.POSITIVE_INFINITY) {
				return (i.endPoint + i.startPoint)/2;
			}
			else {
				// the last guy, just do start point + 1 =)
				return i.startPoint + 1;
			}
		} else {
			// if it is average, we don't know what to do here
			assert (false);
			return null;
		}
	}
	
	public static double[] getCakeFactors (int n, Interval[] intervals, double[] popSizes, CakeStyle cakeStyle) {
		if (cakeStyle != CakeStyle.AVERAGE){ 
			double[] times = getRepresentativeTimes(intervals, cakeStyle);
			return getCakeFactors(n, times, intervals, popSizes);
		} else {
			double[] beginFactors = getCakeFactors(n, intervals, popSizes, CakeStyle.BEGINNING);
			double[] endFactors = getCakeFactors(n, intervals, popSizes, CakeStyle.END);
			assert (beginFactors.length == endFactors.length);
			
			double[] cakeFactors = new double[beginFactors.length];
			for (int i = 0; i < cakeFactors.length; i++) {
				cakeFactors[i] = (beginFactors[i] + endFactors[i]) / 2d;
			}
			
			return cakeFactors;
		}
	}
	
	// ancestral times is where we want to know what the factor is
	// intervals, popSizes are the intervals and population sizes of the demography (single population)
	public static double[] getCakeFactors (int n, double[] ancestralTimes, Interval[] intervals, double[] popSizes) {
		double[] rescaledAncestralTimes = getRescaledAncestralTimes(ancestralTimes, intervals, popSizes);
		
		return TrunkProcess.getCakeFactors(n, rescaledAncestralTimes);
	}
	
	// rescales ancestralTimes to equivalent times under a constant population size coalescent model
	private static double[] getRescaledAncestralTimes (double[] ancestralTimes, Interval[] intervals, double[] popSizes) {
		assert (intervals.length == popSizes.length);
		
		// get the time where the popizes might changes
		double lastEnd = 0d;
		double[] popSizeChangeTimes = new double[intervals.length - 1];
		for (int i=0; i<intervals.length; i++) {
			assert (lastEnd == intervals[i].startPoint);
			lastEnd = intervals[i].endPoint;
			if (lastEnd != Double.POSITIVE_INFINITY) {
				assert (i < intervals.length - 1);
				popSizeChangeTimes[i] = lastEnd;
			}
			else {
				assert (i == intervals.length - 1);
			}
		}
				
		
		// transform changeTimes
		// times in smaller populations come earlier, in large populations later
		double[] transformedChangeTimes = new double[popSizeChangeTimes.length];
		double currRescaledTime = 0d;
		double currTime = 0d;
		for (int i=0; i<popSizeChangeTimes.length; i++) {
			double nextTime = popSizeChangeTimes[i];
			assert (nextTime >= 0d);
			
			double currIncrement = nextTime - currTime;
			currTime = nextTime;
			
			assert(popSizes[i] > 0d);
			double rescaledIncrement = currIncrement / popSizes[i];
			currRescaledTime += rescaledIncrement;

			transformedChangeTimes[i] = currRescaledTime;
		}
		
		// the real sample times
		// have to be rescaled like the change times
		double[] rescaledAncestralTimes = new double[ancestralTimes.length];
		for (int i = 0; i < ancestralTimes.length; i++) {
			double origSamplingTime = ancestralTimes[i];
			assert (origSamplingTime >= 0d);
			
			// find ceiling in orig size change times
			int floorIdx = findFloor (popSizeChangeTimes, origSamplingTime);
			
			double origIncrement = origSamplingTime - (floorIdx == -1 ? 0d : popSizeChangeTimes[floorIdx]);
			double rescaledIncrement = origIncrement / popSizes[floorIdx + 1];
			rescaledAncestralTimes[i] = rescaledIncrement + (floorIdx == -1 ? 0d : transformedChangeTimes[floorIdx]);
		}
		
		return rescaledAncestralTimes;
	}
	
	
	private static int findFloor (double[] array, double number) {
		int i;
		for (i = array.length - 1; i >= 0; i--) {
			if (array[i] <= number) break;
		}
		return i;
	}
	
	public static double[] getCakeFactors (int n, double[] times) {
		
		assert (n >= 0);
		
		// put them here
		double[] cakeFactors = new double[times.length];

		if (n == 0) {
			Arrays.fill(cakeFactors, 1d);
			return cakeFactors;
		}
		
		// the fraction
		double[] theFraction = preComputeExpectedNumberFraction (n);
		
		// now calculate the expectations
		for (int t=0; t<times.length; t++) {
			// divide expected number of lineages by n
			cakeFactors[t] = getExpectedNumberLineages (n, times[t], theFraction) /n;
		}

		// return them
		return cakeFactors;
	}

	
	private static double getExpectedNumberLineages (int n, double time, double[] theFraction) {
		
		double toReturn;
		// maybe special cases?
		if (time == 0d) {
			toReturn = n;
		}
		else if (time == Double.POSITIVE_INFINITY) {
			toReturn = 1d;
		}
		else {
			// normal case
			
			// compute the sum (for the expected value first)
			toReturn = 0d;
			for (int i=1; i<=n; i++) {
				toReturn += Math.exp (-i*(i-1)/2d * time) * (2*i - 1) * theFraction[i];
			}
		}

		return toReturn;
	}

	private static double[] preComputeExpectedNumberFraction (int n) {
		double[] toReturn = new double[n+1];
		
		// calculate it iteratively
		toReturn[0] = 1;
		for (int i=1; i<toReturn.length; i++) {
			toReturn[i] = (n-i+1) / (double)(n+i-1) * toReturn[i-1];
		}
		
		// give it away
		return toReturn;
	}

	public static abstract class MigratingTrunk extends TrunkProcess {

		// [epoch][presentDeme][ancientDeme]
		// prob that lineage in presentDeme is in ancientDeme at chosen point of epoch
		protected final double[][][] migrationTransitionProbs;
		
		// [epoch][presentDeme][ancientDeme]
		private final double[][][] fractionAncientDemeToPresentDeme;
		
		private final CakeStyle cakeStyle;
		
		public MigratingTrunk(int[] sampleSizes, Demography refinedDemography, CakeStyle cakeStyle) {
			super(sampleSizes, refinedDemography);
			
			this.cakeStyle = cakeStyle;
			this.migrationTransitionProbs = new double[this.refinedDemography.epochList.length][][];
			
			double[] trunkPoints = null;
			if (hasRepresentativeTimes(this.cakeStyle)) {
				trunkPoints = TrunkProcess.getRepresentativeTimes(this.refinedDemography.epochList, this.cakeStyle);
			}
			
			ArrayList<MyEigenDecomposition> eigenStuffList = new ArrayList<MyEigenDecomposition>();
			// go through and get them
			for (double[][] thisMigMatrix : this.refinedDemography.migrationMatrixList) {
				eigenStuffList.add(thisMigMatrix == null ? null : new MyEigenDecomposition(thisMigMatrix, UberDemographyCore.ONELOCUS_EPSILON));
			}
			
			// [presentDeme][ancientDeme]
			double[][] prevEndTransition = null;
			
			for (int epoch = 0; epoch < this.refinedDemography.epochList.length; epoch++) {
				assert (trunkPoints == null || trunkPoints[epoch] <= this.refinedDemography.epochList[epoch].endPoint);
				assert (trunkPoints == null || trunkPoints[epoch] >= this.refinedDemography.epochList[epoch].startPoint);
				
				int numAncientDemes = this.refinedDemography.treePartitionList.get(epoch).size();
				
				double[][] startTransition;
				
				// update currTransition up to the startPoint
				if (epoch == 0) {
					
					assert (sampleSizes.length == numAncientDemes);
					
					startTransition = new double[sampleSizes.length][sampleSizes.length];
					for (int deme = 0; deme < sampleSizes.length; deme++) {
						startTransition[deme][deme] = 1d;
					}
					
				} else {
					startTransition = new double[sampleSizes.length][numAncientDemes];
					
					for (int presentDeme = 0; presentDeme < sampleSizes.length; presentDeme++) {						
						for (int ancientDeme = 0; ancientDeme < numAncientDemes; ancientDeme++) {
							for (int memberDeme : this.refinedDemography.getMemberDemesIndices(ancientDeme, epoch)) {
								startTransition[presentDeme][ancientDeme] += prevEndTransition[presentDeme][memberDeme];
							}
						}
					}
				}
				
				double[][] endTransitionUpdate;
				if (this.refinedDemography.isPulse(epoch)) {
					endTransitionUpdate = this.refinedDemography.pulseMigrationMatrixList.get(epoch);
				} else {
					endTransitionUpdate = eigenStuffList.get(epoch).getRateMatrixExponential (this.refinedDemography.epochList[epoch].endPoint - this.refinedDemography.epochList[epoch].startPoint);
				}
				double[][] endTransition = (new Matrix(startTransition)).times(new Matrix(endTransitionUpdate)).getArray();
				
				assert (startTransition.length == endTransition.length);
				
				if (this.refinedDemography.isPulse(epoch)) {
					this.migrationTransitionProbs[epoch] = null;
				} else if (this.cakeStyle == CakeStyle.BEGINNING){
					assert (trunkPoints[epoch] == this.refinedDemography.epochList[epoch].startPoint);
					
					this.migrationTransitionProbs[epoch] = startTransition;
				} else if (this.cakeStyle == CakeStyle.END) {
					assert (trunkPoints[epoch] == this.refinedDemography.epochList[epoch].endPoint);
					
					this.migrationTransitionProbs[epoch] = endTransition;
				} else if (this.cakeStyle == CakeStyle.MIDDLE) {
					
					assert (trunkPoints[epoch] <= this.refinedDemography.epochList[epoch].endPoint);
					assert (trunkPoints[epoch] >= this.refinedDemography.epochList[epoch].startPoint);
					
					double[][] transitionUpdate = eigenStuffList.get(epoch).getRateMatrixExponential (trunkPoints[epoch] - this.refinedDemography.epochList[epoch].startPoint);
					
					this.migrationTransitionProbs[epoch] = (new Matrix(startTransition)).times(new Matrix(transitionUpdate)).getArray();
				} else if (this.cakeStyle == CakeStyle.AVERAGE) {
					
					this.migrationTransitionProbs[epoch] = new double[startTransition.length][];
					for (int presentDeme=0; presentDeme < this.migrationTransitionProbs[epoch].length; presentDeme++) {
						assert (startTransition[presentDeme].length == endTransition[presentDeme].length);
						
						this.migrationTransitionProbs[epoch][presentDeme] = new double[startTransition[presentDeme].length];
						for (int ancientDeme = 0; ancientDeme < this.migrationTransitionProbs[epoch][presentDeme].length; ancientDeme++) {
							this.migrationTransitionProbs[epoch][presentDeme][ancientDeme] = (startTransition[presentDeme][ancientDeme] + endTransition[presentDeme][ancientDeme]) / 2d;
						}
					}
					
				} else {
					assert (false);
				}
				
				prevEndTransition = endTransition;
				
				// done with this epoch
			}
			
			// precompute fraction ancient deme to present deme
			this.fractionAncientDemeToPresentDeme = new double[this.refinedDemography.epochList.length][sampleSizes.length][];
			for (int epoch = 0; epoch < this.refinedDemography.epochList.length; epoch++) {
				
				if (this.refinedDemography.isPulse(epoch)) {
					this.fractionAncientDemeToPresentDeme[epoch] = null;
					continue;
				}
				
				// fill fractions (up to proportionality)
				for (int presentDeme = 0; presentDeme < sampleSizes.length; presentDeme++) {
					// each present deme could have come from a bunch of ancient demes
					this.fractionAncientDemeToPresentDeme[epoch][presentDeme] = new double[this.refinedDemography.popSizes.get(epoch).length];
					for (int ancientDeme = 0; ancientDeme < this.refinedDemography.popSizes.get(epoch).length; ancientDeme++) {
						// proportional to number guys in present deme * prob to get to this ancient deme from the present deme
						this.fractionAncientDemeToPresentDeme[epoch][presentDeme][ancientDeme] = sampleSizes[presentDeme] * this.migrationTransitionProbs[epoch][presentDeme][ancientDeme];
					}
				}
				
				// renormalize so fractions sum to 1
				for (int ancientDeme = 0; ancientDeme < this.refinedDemography.popSizes.get(epoch).length; ancientDeme++) {
					double sum = 0d;
					for (int presentDeme = 0; presentDeme < sampleSizes.length; presentDeme++) {
						sum += this.fractionAncientDemeToPresentDeme[epoch][presentDeme][ancientDeme];
					}
					
					if (sum == 0d) continue;
					
					// now divide
					for (int presentDeme = 0; presentDeme < sampleSizes.length; presentDeme++) {
						this.fractionAncientDemeToPresentDeme[epoch][presentDeme][ancientDeme] /= sum;
					}
				}
				
			}
			
			// DONE
		}
		
		
		// proportional to presentDemeSampleSize * migrationTransitionProbs[interval][presentDeme][ancientDeme]
		@Override
		public double fractionAncientDemeToPresentDeme(int refinedInterval, int presentDeme, int ancientDeme) {
			return this.fractionAncientDemeToPresentDeme[refinedInterval][presentDeme][ancientDeme];
		}

	}
	
	
	public static class SimpleMigratingTrunk extends MigratingTrunk {
		final SimpleTrunk nonMigratingTrunk;
		
		public SimpleMigratingTrunk(int[] sampleSizes, Demography refinedDemography, CakeStyle cakeStyle, SimpleTrunk nonMigratingTrunk) {
			super(sampleSizes, refinedDemography, cakeStyle);
			this.nonMigratingTrunk = nonMigratingTrunk;
		}

		@Override
		public double getAbsorbRate(int refinedInterval, int ancientDeme) {
			
			double sum = 0d;
			
			for (int otherAncientDeme = 0; otherAncientDeme < this.refinedDemography.popSizes.get(refinedInterval).length; otherAncientDeme++) {
				for (int presentDeme : this.refinedDemography.treePartitionList.get(refinedInterval).get(otherAncientDeme)) {
					double presentDemeMass = this.nonMigratingTrunk.getAbsorbRate(refinedInterval, otherAncientDeme) * this.nonMigratingTrunk.fractionAncientDemeToPresentDeme(refinedInterval, presentDeme, otherAncientDeme);	
					sum +=  presentDemeMass * this.migrationTransitionProbs[refinedInterval][presentDeme][ancientDeme];
				}
			}
			
			return sum;
		}
	}	
	
	public static class SinglePopExactTrunk extends SimpleTrunk {
		final double[] cakeFactors;
		
		public SinglePopExactTrunk(int[] sampleSizes, Demography refinedDemo) {
			super (sampleSizes, refinedDemo);
		
			assert (refinedDemo == this.getRefinedDemography());
			int totalSampleSize= 0;
			for (int n : sampleSizes){
				totalSampleSize += n;
			}
			if (refinedDemo == null) {
				assert (totalSampleSize == 0);
				cakeFactors = null;
				return;
			}
						
			double[] popSizes = new double[refinedDemo.epochList.length];
			for (int i = 0; i < popSizes.length; i++) {
				popSizes[i] = effectivePopSize(refinedDemo.popSizes.get(i), refinedDemo.migrationMatrixList.get(i));
			}
			
			cakeFactors = exactCoalescentCakeFactors (totalSampleSize, this.refinedDemography.epochList, popSizes);
		}
		
		@Override
		public double getAbsorbRate(int refinedInterval, int ancientDeme) {
			return cakeFactors[refinedInterval] * this.getAncientDemeSize(refinedInterval, ancientDeme) / this.refinedDemography.popSizes.get(refinedInterval)[ancientDeme];
		}
		
		protected double effectivePopSize (double[] popSizes, double[][] migMatrix) {
			
			if (popSizes == null) return 1d;
			
			double sum = 0d;
			for (int i = 0; i < popSizes.length; i++) {
				sum += popSizes[i];
			}
			return sum / popSizes.length;
		}
	}
	
	// returns cake factors so that the absorption probabilities in each epoch matches the true coalescent
	// for the final epoch, takes the expected time until absorption and finds exponential rate that matches it
	private static double[] exactCoalescentCakeFactors (int nTrunk, Interval[] epochs, double[] popSizes) {
		double[] endPoints = getRepresentativeTimes(epochs, CakeStyle.END);
		double[] rescaledEndPoints = getRescaledAncestralTimes(endPoints, epochs, popSizes);
		
		double[] cakeFactors = new double[epochs.length];
		// do all the epochs except the last one
		double lastAbsorbProb = 0d;
		double lastRescaledTime = 0d;
		for (int epochIdx = 0; epochIdx < epochs.length - 1; epochIdx++) {
			
			double epochLength = (epochs[epochIdx].endPoint - epochs[epochIdx].startPoint);
			
			// probability of absorption before the end of the epoch
			double currAbsorbProb = exactCoalescentAbsorptionProb(nTrunk, rescaledEndPoints[epochIdx]);
			
			// conditional probability of absorption, given no absorption before epoch
			double conditionalAbsorbProb = (currAbsorbProb - lastAbsorbProb) / (1d - lastAbsorbProb);
			
			// invert 1 - e^(-totalAbsorbRate) = conditionalAbsorbProb
			double totalAbsorbRate = - Math.log(1d - conditionalAbsorbProb);
			
			// set cake factor
			if (epochLength > 0d) {
				cakeFactors[epochIdx] = totalAbsorbRate / nTrunk /  epochLength * popSizes[epochIdx];
			} else {
				assert conditionalAbsorbProb < UberDemographyCore.ONELOCUS_EPSILON * 1000;
				cakeFactors[epochIdx] = 1d;
			}
			
			assert (cakeFactors[epochIdx] > 0d - UberDemographyCore.ONELOCUS_EPSILON * 1000);
			assert (cakeFactors[epochIdx] <= 1d + UberDemographyCore.ONELOCUS_EPSILON * 1000);
			
			cakeFactors[epochIdx] = Math.max(0d, cakeFactors[epochIdx]);
			cakeFactors[epochIdx] = Math.min(1d, cakeFactors[epochIdx]);

			
			// update lastAbsorbProb
			lastAbsorbProb = currAbsorbProb;
			lastRescaledTime = rescaledEndPoints[epochIdx];
		}
		
		// do the last epoch
		int lastEpochIdx = epochs.length - 1;
		double lastEpochExpectTime = 0d;
		
		// kTrunk = trunk size at beginning of last epoch
		// kTrunkExpectTime = expected time to absorption given kTrunk
		double kTrunkExpectTime = 0d;
		for (int kTrunk = 1; kTrunk <= nTrunk; kTrunk++) {
			// update kTrunkExpectTime
			kTrunkExpectTime = 2d / ((kTrunk+1) * kTrunk) + (kTrunk + 1 - 2) / (double) (kTrunk + 1) * kTrunkExpectTime;
			
			// get the probability of having not been absorbed, and having kTrunk at beginning of last epoch
			double kTrunkProb = probNoAbsorb(nTrunk, nTrunk - kTrunk) * AncestralProcessProbs.getProb(nTrunk + 1, kTrunk + 1, lastRescaledTime);
			
			lastEpochExpectTime += kTrunkExpectTime * kTrunkProb;
		}
		
		// divide to get conditional expectation, given not having been absorbed by beginning of last epoch
		lastEpochExpectTime /= (1d - lastAbsorbProb);
		cakeFactors[lastEpochIdx] = 1d / lastEpochExpectTime / nTrunk;
		
		return cakeFactors;
	}
	
	// the probability of conditional lineage not having been absorbed, given numCoal coalescent events hitting the total sample of size nTrunk + 1
	private static double probNoAbsorb (int nTrunk, int numCoal) {
		double survivalProb = 1d;
		for (int sampleSize = nTrunk + 1; sampleSize > nTrunk + 1 - numCoal; sampleSize--) {
			assert (sampleSize > 1);
			// when we have sampleSize lineages left and a coalescence happens, multiply by probability new lineage was not picked
			survivalProb *= (double) (sampleSize - 2) / (double) sampleSize;
		}			
		
		assert (survivalProb >= 0d);
		assert (survivalProb <= 1d);
		
		return survivalProb;
	}
	
	
	// returns the exact probability that the new guy has coalesced with a trunk lineage before time t 
	private static double exactCoalescentAbsorptionProb(int nTrunk, double t) {
		assert (nTrunk > 0);
		
		double absorptionProb = 0d;
		// sum over the number of coalescent events that have hit the whole sample (both trunk and new guy)
		for (int numCoal = 0; numCoal <= nTrunk; numCoal++) {
			// the probability that the new guy never coalesced (i.e., 'survived'), given that numCoal coalescent events happened
			double survivalProb = probNoAbsorb(nTrunk, numCoal);
			
			absorptionProb += AncestralProcessProbs.getProb(nTrunk + 1, nTrunk + 1 - numCoal, t) * (1d - survivalProb);
		}
		
		if (absorptionProb > - 1e-18 && absorptionProb < 0d){
			absorptionProb = 0d;	//round if there are minor numerical errors.
		}
		
		assert (absorptionProb >= 0d);
		assert (absorptionProb <= 1d);
		
		return absorptionProb;
	}
}
