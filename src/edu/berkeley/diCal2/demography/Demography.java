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

package edu.berkeley.diCal2.demography;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeSet;

import edu.berkeley.diCal2.utility.RealPartition;
import edu.berkeley.diCal2.utility.RealPartition.Interval;

public class Demography {
	
	// the epochs in this demography
	public Interval[] epochList;
	// partitions giving the tree structure of this demography
	public ArrayList<ArrayList<TreeSet<Integer>>> treePartitionList;
	// the migration rates in this demography (a list of them for each epoch)
	public ArrayList<double[][]> migrationMatrixList;
	// and the sizes of the different demes in each epoch
	public ArrayList<double[]> popSizes;
	// and the list of pulseMigrationMatrices
	public ArrayList<double[][]> pulseMigrationMatrixList;
	// and maybe some exponential rates
	public ArrayList<double[]> expGrowthRates;

	
	// empty constructor somehow necessary
	public Demography () {
		// und alle so: yeah
	}
	
	// the explicit constructor
	public Demography (Interval[] parEpochList, ArrayList<ArrayList<TreeSet<Integer>>> parTreePartitionList, ArrayList<double[]> parPopSizes, ArrayList<double[][]> parPulseMigrationMatrixList, ArrayList<double[][]> parMigrationMatrixList) {
		
		// make deep copies of stuff
		this.epochList = deepCopyIntervals (parEpochList);
		this.treePartitionList = deepCopyPartitionList (parTreePartitionList);
		this.migrationMatrixList = deepCopyMatrixList (parMigrationMatrixList);
		this.pulseMigrationMatrixList = deepCopyMatrixList (parPulseMigrationMatrixList);
		this.popSizes = deepCopyPopSizes (parPopSizes);
		// nothing here
		this.expGrowthRates = null;
		// should be done
		
		for (int i = 0; i < this.epochList.length; i++) {
			// contains assertions we want to check
			this.isPulse(i);
		}
	}
	
	// with exp rates
	public Demography (Interval[] parEpochList, ArrayList<ArrayList<TreeSet<Integer>>> parTreePartitionList, ArrayList<double[]> parPopSizes, ArrayList<double[]> parExpRates, ArrayList<double[][]> parPulseMigrationMatrixList, ArrayList<double[][]> parMigrationMatrixList) {
		// first the old stuff
		this (parEpochList, parTreePartitionList, parPopSizes, parPulseMigrationMatrixList, parMigrationMatrixList);
		
		// and then the exp rates, if we have some
		if (parExpRates != null) {
			// we should be ok with using the popsize deep-copy
			this.expGrowthRates = deepCopyPopSizes (parExpRates);
		}
		
	}
	
	protected static ArrayList<double[]> deepCopyPopSizes (ArrayList<double[]> parPopSizes) {
		ArrayList<double[]> returnSizes = new ArrayList<double[]>();
		for (double[] addPopSize : parPopSizes) {
			
			double[] toAdd = addPopSize == null ? null : Arrays.copyOf (addPopSize, addPopSize.length);
			returnSizes.add (toAdd);
		}
		return returnSizes;
	}

	protected static ArrayList<double[][]> deepCopyMatrixList (ArrayList<double[][]> parMigrationMatrixList) {
		ArrayList<double[][]> returnList = new ArrayList<double[][]>();
		for (double[][] origMatrix : parMigrationMatrixList) {
			if (origMatrix == null) {
				returnList.add (null);
			}
			else {
				double[][] addMatrix = new double[origMatrix.length][];
				for (int i=0; i<origMatrix.length; i++) {
					assert (origMatrix[i].length == origMatrix.length);
					addMatrix[i] = Arrays.copyOf (origMatrix[i], origMatrix[i].length);
				}
				returnList.add (addMatrix);
			}
		}
		return returnList;
	}

	protected static ArrayList<ArrayList<TreeSet<Integer>>> deepCopyPartitionList (ArrayList<ArrayList<TreeSet<Integer>>> parTreePartitionList) {
		ArrayList<ArrayList<TreeSet<Integer>>> returnList = new ArrayList<ArrayList<TreeSet<Integer>>>();
		for (ArrayList<TreeSet<Integer>> origPartition : parTreePartitionList) {
			ArrayList<TreeSet<Integer>> addPartition = new ArrayList<TreeSet<Integer>>();
			for (TreeSet<Integer> origBlock : origPartition) {
				TreeSet<Integer> addBlock = new TreeSet<Integer>();
				for (int addInt : origBlock) {
					addBlock.add (addInt);
				}
				addPartition.add (addBlock);
			}
			returnList.add (addPartition);
		}

		return returnList;
	}

	protected static Interval[] deepCopyIntervals (Interval[] parEpochList) {
		Interval[] returnList = new Interval[parEpochList.length];
		for (int i=0; i<returnList.length; i++) {
			returnList[i] = new Interval (parEpochList[i].startPoint, parEpochList[i].endPoint);
		}
		return returnList;
	}

	// deep-copy constructor
	public Demography (Demography demo) {
		this (demo.epochList, demo.treePartitionList, demo.popSizes, demo.expGrowthRates, demo.pulseMigrationMatrixList, demo.migrationMatrixList);
	}
	
	public Demography (Demography demo, Interval[] stateIntervals, ArrayList<Integer> refinedIntervalToStateInterval, double EPSILON) {
		this (demo, stateIntervals, refinedIntervalToStateInterval, null, EPSILON);
	}
	
	// refines the given demography (making a deep-copy of it)
	public Demography (Demography demo, Interval[] stateIntervals, ArrayList<Integer> refinedIntervalToStateInterval, ArrayList<Integer> refinedIntervalToDemoInterval, double EPSILON) {

		if (stateIntervals == null) {
			stateIntervals = new Interval[] {new Interval(0d, Double.POSITIVE_INFINITY)};
		}
		
		
		assert (RealPartition.isConsecutive(stateIntervals));
		// they has to be really bigger, otherwise the remainder doesn't quite work
		assert (RealPartition.allBiggerThanEpsilon(stateIntervals, 10*EPSILON));
		
		
		if (refinedIntervalToStateInterval == null) refinedIntervalToStateInterval = new ArrayList<Integer>();
		else refinedIntervalToStateInterval.clear();
		
		if (refinedIntervalToDemoInterval == null) refinedIntervalToDemoInterval = new ArrayList<Integer>();
		else refinedIntervalToDemoInterval.clear();
		
		// old things
		Interval[] demoIntervals = demo.epochList;
		// also refineByIntervals
		assert (RealPartition.isConsecutive(demoIntervals));
		// they has to be really bigger, otherwise the remainder doesn't quite work
		assert (demo.checkIntervalLengths(10*EPSILON));
		
		// some checking
		assert (demoIntervals[0].startPoint == stateIntervals[0].startPoint);
		assert (demoIntervals[0].startPoint == 0d);
		assert (demoIntervals[demoIntervals.length-1].endPoint == stateIntervals[stateIntervals.length-1].endPoint);
		assert (demoIntervals[demoIntervals.length-1].endPoint == Double.POSITIVE_INFINITY);

		
		ArrayList<Interval> newIntervals = new ArrayList<Interval>();
		
		// get the new intervals and map them to the original ones
		double currStart = 0d;
		double currEnd = 0d;
		
		int curDemoIntervalIdx = 0;
		int curStateIntervalIdx = 0;
		
		assert (demoIntervals.length > 0);
		assert (stateIntervals.length > 0);
		
		// now go through
		while (currStart < Double.POSITIVE_INFINITY) {
			double demoIntervalEnd = demoIntervals[curDemoIntervalIdx].endPoint;
			double refineByIntervalEnd = stateIntervals[curStateIntervalIdx].endPoint;
			
			refinedIntervalToDemoInterval.add(curDemoIntervalIdx);
			refinedIntervalToStateInterval.add(curStateIntervalIdx);

			// Note: an evil demon might make the endpoints the same
			if (demoIntervalEnd <= refineByIntervalEnd - EPSILON){
				// demoInterval gives end, so advance demo by one
				currEnd = demoIntervalEnd;
				curDemoIntervalIdx++;
			}
			else if (refineByIntervalEnd <= demoIntervalEnd - EPSILON){
				// refinedInterval gives end, so advance refined by one
				currEnd = refineByIntervalEnd;
				curStateIntervalIdx++;
			}
			else {
				// both give end, so advance both by one [although we take demo as offical end]
				currEnd = demoIntervalEnd;
				curDemoIntervalIdx++;
				curStateIntervalIdx++;
			}
			
			newIntervals.add(new Interval(currStart, currEnd));			
			currStart = currEnd;			
		}
		
		this.epochList = newIntervals.toArray(new Interval[]{});
		assert (refinedIntervalToDemoInterval.size() == this.epochList.length);
		
		// now copy the stuff from the original demography appropriately
		// the treePartionList
		ArrayList<ArrayList<TreeSet<Integer>>> newPartitionList = new ArrayList<ArrayList<TreeSet<Integer>>>();
		// then the migration matrices
		ArrayList<double[][]> newMigrationList = new ArrayList<double[][]>();
		// and the popSizes
		ArrayList<double[]> newPopSizes = new ArrayList<double[]>();
		// the growth rates
		ArrayList<double[]> newGrowthRates = null;
		if (demo.expGrowthRates != null) {
			newGrowthRates = new ArrayList<double[]>();
		}
		// and the pulse migration
		ArrayList<double[][]> newPulseList = new ArrayList<double[][]>();
		// build the new guys
		int prevEpoch = -1;
		for (Integer i : refinedIntervalToDemoInterval) {
			
			if (demo.isPulse(i)) {
				assert (i == prevEpoch + 1);
			} else {
				assert (i == prevEpoch + 1 || i == prevEpoch);
			}
			
			newPartitionList.add (demo.treePartitionList.get(i));
			newMigrationList.add (demo.migrationMatrixList.get(i));
			// what about them growth rates being special
			if (demo.expGrowthRates == null) {
				// only the sizes
				newPopSizes.add (demo.popSizes.get(i));
			}
			else {
				// growth rates and sizes
				// growth rates stay the same in each interval
				double[] currGrowthRates = demo.expGrowthRates.get(i);
				newGrowthRates.add (currGrowthRates);
				// the sizes might have to be modified
				// odd way to get the right one
				double demoEpochEnd = demo.epochList[i].endPoint;
				double refEpochEnd = this.epochList[newPopSizes.size()].endPoint;
				if ((demoEpochEnd == refEpochEnd) || (demo.popSizes.get(i) == null)) {
					// no time shifted, so don't modify size
					newPopSizes.add (demo.popSizes.get(i));
				}
				else {
					// time shifted, so modify sizes
					double elapsedTime = demoEpochEnd - refEpochEnd;
					assert (elapsedTime > 0);
					double[] finalSizes = demo.popSizes.get(i);
					double[] newEndSizes = new double[finalSizes.length];
					assert (finalSizes.length == currGrowthRates.length);
					for (int k=0; k<finalSizes.length; k++) {
						if (Math.abs(currGrowthRates[k]) < EPSILON) {
							newEndSizes[k] = finalSizes[k];
						}
						else {
							// this would not make sense
							assert (elapsedTime != Double.POSITIVE_INFINITY);
							newEndSizes[k] = Math.exp(elapsedTime * currGrowthRates[k]) * finalSizes[k];
						}
					}
					// add them
					newPopSizes.add (newEndSizes);
				}
			}
			newPulseList.add (demo.pulseMigrationMatrixList.get(i));
			
			prevEpoch = i;
		}
		// really copy it
		this.treePartitionList = deepCopyPartitionList (newPartitionList);
		this.migrationMatrixList = deepCopyMatrixList (newMigrationList);
		this.popSizes = deepCopyPopSizes (newPopSizes);
		if (newGrowthRates != null) {
			this.expGrowthRates = deepCopyPopSizes(newGrowthRates);
		}
		else {
			this.expGrowthRates = null;
		}
		this.pulseMigrationMatrixList = deepCopyMatrixList (newPulseList);
		// should be done
	}
	
	public boolean isPulse (Integer epoch) {
		if (this.epochList[epoch].endPoint == this.epochList[epoch].startPoint){
			// should be a pulse
			assert ((this.pulseMigrationMatrixList.get(epoch) != null) && (this.migrationMatrixList.get(epoch) == null) && (this.popSizes.get(epoch) == null));
			if(this.expGrowthRates != null){assert(this.expGrowthRates.get(epoch) == null);}
			return true;
		}
		else {
			// should not be a pulse
			assert ((this.pulseMigrationMatrixList.get(epoch) == null) && (this.migrationMatrixList.get(epoch) != null) && (this.popSizes.get(epoch) != null));
			if(this.expGrowthRates != null){assert(this.expGrowthRates.get(epoch) != null);}
			return false;
		}
	}

	protected boolean checkIntervalLengths (double EPSILON) {
		// got through intervals
		for (int i=0; i<this.epochList.length; i++) {
			double thisLength = this.epochList[i].endPoint - this.epochList[i].startPoint;
			if (thisLength < EPSILON) {
				if (!this.isPulse(i)) return false;
			}
			else {
				if (this.isPulse(i)) return false;
			}
		}
		// we went through, so everything fine
		return true;
	}

	public void dump (PrintStream outStream) {
		// go through the epochs
		for (int e=0; e<this.epochList.length; e++) {
			outStream.println ("==== Epoch " + e + " ====");
			// say the interval
			outStream.println (this.epochList[e].startPoint + "\t" + this.epochList[e].endPoint);
			// say the partition
			outStream.print("{");
			Iterator<TreeSet<Integer>> itBlock = this.treePartitionList.get(e).iterator();
			while (itBlock.hasNext()) {
				TreeSet<Integer> block = itBlock.next();
				outStream.print("{");
				Iterator<Integer> itElement = block.iterator();
				while (itElement.hasNext()) {
					Integer element = itElement.next();
					outStream.print(element);
					if (itElement.hasNext()) {
						outStream.print (", ");
					}
				}
				if (itBlock.hasNext()) {
					outStream.print ("}, ");
				}
				else {
					outStream.print ("}");
				}
			}
			outStream.println ("}");
			// say the sizes
			if (this.popSizes.get(e) == null) {
				outStream.print("null");
			} else {
				for (double currSize : this.popSizes.get(e)) {
					outStream.print (currSize + "\t");
				}
			}
			outStream.println();
			
			// maybe say the exp rates, if there are some
			outStream.print ("rates:\t");
			if (this.expGrowthRates == null || (this.expGrowthRates.get(e) == null)) {
				outStream.print("null");
			} else {
				for (double currRate : this.expGrowthRates.get(e)) {
					outStream.print (currRate + "\t");
				}
			}
			outStream.println();
			
			// say the pulse matrix
			if (this.pulseMigrationMatrixList.get(e) == null) {
				outStream.println("null");
			} else {
				for (double[] currRow : this.pulseMigrationMatrixList.get(e)) {
					for (double currVal : currRow) {
						outStream.print (currVal + "\t");
					}
					outStream.println();
				}
			}
			
			// say the mig matrix
			if (this.migrationMatrixList.get(e) == null) {
				outStream.println("null");
			} else {
				for (double[] currRow : this.migrationMatrixList.get(e)) {
					for (double currVal : currRow) {
						outStream.print (currVal + "\t");
					}
					outStream.println();
				}
			}

		}
	}
	
	public TreeSet<Integer> getMemberDemesIndices(int demeIdx, int epochIdx) {
		assert (epochIdx > 0);
		// return set
		TreeSet<Integer> returnSet = new TreeSet<Integer>();
		// loop over all previous guys
		for (int prevDeme=0; prevDeme<this.treePartitionList.get(epochIdx-1).size(); prevDeme++) {
			// remember index if its in there
			if (this.treePartitionList.get(epochIdx).get(demeIdx).containsAll(this.treePartitionList.get(epochIdx-1).get(prevDeme))) {
				returnSet.add (new Integer(prevDeme));
			}
		}
		// return the set
		return returnSet;
	}

	// just some wrappers to make it easier
	public int numIntervals() {
		return this.epochList.length;
	}

	public int numAncientDemes(int e) {
		return this.treePartitionList.get(e).size();
	}

	// leftLimit means that we keep a cut epoch of length 0
	public Demography cutDemography(double time, boolean discardZeroIntervals) {
		int startEpoch = 0;
		while (epochList[startEpoch].endPoint < time) startEpoch++;
		if (discardZeroIntervals && epochList[startEpoch].endPoint == time) startEpoch++;
		
		Interval[] cutEpochList = new Interval[this.epochList.length - startEpoch];
		for (int i = startEpoch; i < this.epochList.length; i++) {
			cutEpochList[i - startEpoch] = new Interval(Math.max(0, this.epochList[i].startPoint - time), this.epochList[i].endPoint - time);
		}
		
		ArrayList<double[][]> cutMigMatrixList = removeFirstEntries(deepCopyMatrixList (migrationMatrixList), startEpoch);
		ArrayList<double[][]> cutPulseMigList = removeFirstEntries(deepCopyMatrixList (pulseMigrationMatrixList), startEpoch);
		ArrayList<double[]> cutPopSizesList = removeFirstEntries(deepCopyPopSizes (popSizes), startEpoch);
		ArrayList<double[]> cutExpGrowthRates = expGrowthRates == null ? null : removeFirstEntries(deepCopyPopSizes(expGrowthRates), startEpoch);
		
		if (cutEpochList[0].endPoint == cutEpochList[0].startPoint && cutPulseMigList.get(0) == null) {
			assert !discardZeroIntervals;
			
			int numPops = numAncientDemes(startEpoch);
			double[][] tempMat= new double[numPops][numPops];
			for(int j= 0; j < numPops; j++){
				tempMat[j][j]= 1;
			}
			cutPulseMigList.set(0, tempMat);
			cutMigMatrixList.set(0, null);
			cutPopSizesList.set(0, null);
			if (cutExpGrowthRates != null) cutExpGrowthRates.set(0, null);
		}
		
		// make tree partition list
		ArrayList<ArrayList<TreeSet<Integer>>> cutPartitionList = new ArrayList<ArrayList<TreeSet<Integer>>>();
		// make the first list of sets
		ArrayList<TreeSet<Integer>> bottomLeafsPartition = new ArrayList<TreeSet<Integer>>();
		for (int pop = 0; pop < numAncientDemes(startEpoch); pop++) {
			bottomLeafsPartition.add(new TreeSet<Integer>());
			bottomLeafsPartition.get(pop).add(pop);
		}
		cutPartitionList.add(bottomLeafsPartition);
		
		// go thru and make the other partitions
		for (int epoch = startEpoch+1; epoch < this.epochList.length; epoch++) {
			ArrayList<TreeSet<Integer>> prevPartition = cutPartitionList.get(epoch - startEpoch - 1);
			ArrayList<TreeSet<Integer>> currPartition = new ArrayList<TreeSet<Integer>>();
			
			for (int currDeme = 0; currDeme < this.treePartitionList.get(epoch).size(); currDeme++) {
				TreeSet<Integer> currSet = new TreeSet<Integer>();
				for (int memberDeme : this.getMemberDemesIndices(currDeme, epoch)) {
					currSet.addAll(prevPartition.get(memberDeme));
				}
				currPartition.add(currSet);
			}
			
			cutPartitionList.add(currPartition);
		}
		
		return new Demography(cutEpochList, cutPartitionList, cutPopSizesList, cutExpGrowthRates, cutPulseMigList, cutMigMatrixList);
	}
	
	private static <E> ArrayList<E> removeFirstEntries(ArrayList<E> list, int numEntries) {
		for (int i = 0; i < numEntries; i++) {
			list.remove(0);
		}
		return list;
	}
	
	public double getPopSize(int demeIdx, int epochIdx, double currTime) {
		
		if (this.popSizes.get(epochIdx) == null) {
			throw new RuntimeException("Asked for undefined population size during an infinitessimal admixture epoch.");
		}
		if (currTime > this.epochList[epochIdx].endPoint || currTime < this.epochList[epochIdx].startPoint) {
			throw new RuntimeException("Asked for population size at invalid time.");
		}

		
		double ancestralPop = this.popSizes.get(epochIdx)[demeIdx];
		
		double growthRate = this.expGrowthRates == null ? 0 : this.expGrowthRates.get(epochIdx)[demeIdx];
		if(this.epochList[epochIdx].endPoint == Double.POSITIVE_INFINITY){
			assert growthRate == 0;
			return(ancestralPop);
		}
		
		return ancestralPop*Math.exp(growthRate*(this.epochList[epochIdx].endPoint - currTime));
	}

	
	public double getTwiceMigrationRate(int epochIndex, int childDeme, int parentDeme) {
		if(!this.isPulse(epochIndex)){
			return this.migrationMatrixList.get(epochIndex)[childDeme][parentDeme];
		}
		else{throw new RuntimeException("Don't ask for what you can't have!");}
	}

	public void addTrivialExponentialGrowthRates() {
		assert this.expGrowthRates == null;
		this.expGrowthRates = new ArrayList<double[]>();
		for (int epoch = 0; epoch < this.epochList.length; epoch++) {
			this.expGrowthRates.add(new double[this.popSizes.get(epoch).length]);
		}
	}

	public boolean hasPulseMigration() {
		for (int epoch=0; epoch<this.numIntervals(); epoch++) {
			if (this.isPulse(epoch)) {
				return true;
			}
		}
		return false;
	}
}
