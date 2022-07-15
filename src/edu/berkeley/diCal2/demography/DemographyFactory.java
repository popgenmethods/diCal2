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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.berkeley.diCal2.utility.RateMatrixTools;
import edu.berkeley.diCal2.utility.SimplePair;
import edu.berkeley.diCal2.utility.RealPartition.Interval;

public abstract class DemographyFactory {
	
	protected MutableDouble[] currPoint = null;
	
	// the epochs in this demography
	protected final ArrayList<SimplePair<MutableDouble, MutableDouble>> currEpochList;
	// partitions giving the tree structure of this demography
	protected final ArrayList<ArrayList<TreeSet<Integer>>> treePartitionList;
	// the migration rates in this demography (a list of them for each epoch)
	protected final ArrayList<MutableDouble[][]> currMigrationMatrixList;
	// and the sizes of the different demes in each epoch
	protected final ArrayList<MutableDouble[]> currPopSizesList;
	// and some pulse migration matrices
	protected final ArrayList<MutableDouble[][]> currPulseMigrationMatrixList;
	// possibly some exponential rates
	// no final here for now
	protected ArrayList<MutableDouble[]> currExpRatesList;
	
	// them bounds
	private final double[][] bounds;
	
	// and some epsilon
	public final double EPSILON;

	// half the migration rate? (e.g. for migrating trunk)
	protected final boolean halfMigrationRate;
	
	public int getNumPresentDemes () {
		return currPopSizesList.get(0).length;
	}
	
	// we want nothing too fancy for now
	// in the future a cemetary state should solve this issue
	public static boolean isReasonable (Demography demo, double epsilon) {
		if (demo == null) return false;
		
		// some partitions
		ArrayList<TreeSet<Integer>> lastPartition = demo.treePartitionList.get(demo.treePartitionList.size()-1);
		// is the last one the whole thing?
		if (lastPartition.size() == 1) {
			// the inclusion from somewhere else should make sure that this is the whole set
			return true;
		}
		else {
			// get the last migration matrix
			double[][] theMigration =  demo.migrationMatrixList.get(demo.migrationMatrixList.size()-1);
			// can we leave every deme?
			for (int i=0; i<theMigration.length; i++) {
				// can we leave this one?
				if (theMigration[i][i] > -epsilon) return false;
			}
			// we made it
			return true;
		}
	}
	
	public static class MutableDouble {
		private double value;

		public MutableDouble (double value) {
			this.value = value;
		}

		public double getValue() {
			return value;
		}

		public void setValue(double value) {
			this.value = value;
		}
		
		public String toString () {
			return Double.toString(this.value); 
		}
	}
	
	
	public DemographyFactory (double EPSILON, boolean halfMigrationRate, double[][] bounds) {
		this.currEpochList = new ArrayList<SimplePair<MutableDouble,MutableDouble>>();
		this.treePartitionList = new ArrayList<ArrayList<TreeSet<Integer>>>();
		this.currMigrationMatrixList = new ArrayList<DemographyFactory.MutableDouble[][]>();
		this.currPopSizesList = new ArrayList<DemographyFactory.MutableDouble[]>();
		this.currPulseMigrationMatrixList = new ArrayList<DemographyFactory.MutableDouble[][]>();
		this.EPSILON = EPSILON;
		this.halfMigrationRate = halfMigrationRate;
		this.bounds = bounds;
	}
	
	public int getDimension () {
		assert ((this.bounds == null) || (this.bounds.length == this.currPoint.length));
		return this.currPoint.length;
	}
		
	synchronized public Demography getDemography (double[] point) {
		setCurrPoint(point);
		Demography toReturn = getCurrDemography();
		return toReturn;
	}
	
	private void setCurrPoint(double[] point) {
		assert (point.length == getDimension());
		
		for (int i = 0; i < point.length; i++){
			this.currPoint[i].setValue(point[i]);
		}
	}
	
	private Demography getCurrDemography () {
		
		if (!currTimesValid()) return null;
		
		// check whether we are in the bounds
		if (this.bounds != null) {
			assert (this.currPoint.length == this.bounds.length);
			for (int i=0; i< this.bounds.length; i++) {
				if ((this.currPoint[i].value < this.bounds[i][0]) || (this.bounds[i][1] < this.currPoint[i].value)) {
					// out of bounds
					return null;
				}
			}
		}
				
		int numEpochs = currEpochList.size();
		assert(treePartitionList.size() == numEpochs);
		assert(currMigrationMatrixList.size() == numEpochs);
		assert(currPopSizesList.size() == numEpochs);
		assert(this.currPulseMigrationMatrixList.size() == numEpochs);
		assert (this.currExpRatesList == null || this.currExpRatesList.size() == numEpochs);
		
		// some temps
		Interval[] epochList = new Interval[numEpochs];
		ArrayList<double[]> popSizesList = new ArrayList<double[]>();
		// dem rates are special (for now)
		ArrayList<double[]> expRatesList = null;
		if (this.currExpRatesList != null) {
			expRatesList = new ArrayList<double[]>();
		}
		ArrayList<double[][]> migrationMatrixList = new ArrayList<double[][]>();
		ArrayList<double[][]> pulseMigrationList = new ArrayList<double[][]>();
		
		// make it happen
		for (int epoch = 0; epoch < numEpochs; epoch++) {
			// add to epochList
			SimplePair<MutableDouble, MutableDouble> currEpoch = currEpochList.get(epoch);			
			epochList[epoch] = new Interval(currEpoch.first().getValue(), currEpoch.second().getValue());
			
			int numPops = this.treePartitionList.get(epoch).size();
			
			// add to popsizes (only if it is not a pulse interval)
			MutableDouble[] mutablePopSizes = currPopSizesList.get(epoch);
			if (mutablePopSizes != null) {
				double[] popSizes = new double[numPops];
				popSizesList.add(popSizes);
				for (int pop=0; pop < numPops; pop++){
					popSizes[pop] = mutablePopSizes[pop].getValue();
					if (popSizes[pop] <= EPSILON) return null;
				}
			}
			else {
				popSizesList.add(null);
			}
			
			// some exponential rates, if we want
			if (this.currExpRatesList != null) {
				MutableDouble[] mutableRates = this.currExpRatesList.get(epoch);
				if (mutableRates != null) {
					assert (mutableRates.length == numPops);
					double[] expRates = new double[numPops];
					expRatesList.add(expRates);
					for (int pop=0; pop < numPops; pop++){
						expRates[pop] = mutableRates[pop].getValue();
						// rates can be positive and negative, and 0 is actually very ok
					}
				}
				else {
					// null allowed cause of pulse epochs
					expRatesList.add(null);
				}
				// done
			}
			
			// add to migration matrix list
			MutableDouble[][] mutableMigMatrix = currMigrationMatrixList.get(epoch);
			if (mutableMigMatrix == null) {
				migrationMatrixList.add(null);
			}
			else {
				assert (mutableMigMatrix.length == numPops);
				
				double[][] migMatrix = new double[numPops][numPops];
				migrationMatrixList.add(migMatrix);
				
				for (int row = 0; row < numPops; row++){
					MutableDouble[] mutableRow = mutableMigMatrix[row];
					assert (mutableRow.length == numPops);
					
					double totalMigRate = 0d;
					for (int column = 0; column < numPops; column++) {
						if (row == column) continue;
						totalMigRate += (migMatrix[row][column] = mutableRow[column].getValue());
					}
					migMatrix[row][row] = - totalMigRate;
				}
				if (!RateMatrixTools.isProperRateMatrix (migMatrix, EPSILON)) return null;
			}
				
			// add to pulse migration matrix list
			MutableDouble[][] mutablePulseMatrix = this.currPulseMigrationMatrixList.get(epoch);
			if (mutablePulseMatrix == null) {
				pulseMigrationList.add(null);
			}
			else {
				assert (mutablePulseMatrix.length == numPops);
				
				double[][] pulseMatrix = new double[numPops][numPops];
				pulseMigrationList.add(pulseMatrix);
				
				for (int row = 0; row < numPops; row++){
					MutableDouble[] mutableRow = mutablePulseMatrix[row];
					assert (mutableRow.length == numPops);
					
					double totalPulseProb = 0d;
					for (int column = 0; column < numPops; column++) {
						if (row == column) continue;
						totalPulseProb += (pulseMatrix[row][column] = mutableRow[column].getValue());
					}
					pulseMatrix[row][row] = 1 - totalPulseProb;
				}
				if (!RateMatrixTools.isProperStochaticMatrix (pulseMatrix, EPSILON)) return null;
			}

		}
		
		if (this.halfMigrationRate) {
			for (double[][] migMatrix : migrationMatrixList) {
				
				if (migMatrix == null) continue;
				
				for (int i = 0; i < migMatrix.length; i++) {
					for (int j = 0; j < migMatrix[i].length; j++) {
						migMatrix[i][j] /= 2d;
					}
				}
			}
		}
		
		// get demography
		return  new Demography (epochList, treePartitionList, popSizesList, expRatesList, pulseMigrationList, migrationMatrixList);
	}
	
	public boolean currTimesValid() {
		assert (this.currEpochList != null);
		assert (this.currEpochList.size() > 0);
		
		if (this.currEpochList.get(0).first().getValue() != 0d) return false;
		if (this.currEpochList.get(this.currEpochList.size() - 1).second().getValue() != Double.POSITIVE_INFINITY) return false;
		
		double lastPoint = 0d;
		for (int i=0; i<this.currEpochList.size(); i++){
			// get pair
			SimplePair<MutableDouble, MutableDouble> pair = this.currEpochList.get(i);
			
			assert (pair.first().getValue() == lastPoint);
			
			if (this.isPulse  (i)) {
				// they have to be equal
				if (pair.second().getValue() != pair.first().getValue()) return false;
			}
			else {
				// they should be at least EPSILON apart
				if (pair.second().getValue() <= pair.first().getValue() + EPSILON) return false;
			}
			
			// prepare next round
			lastPoint = pair.second().getValue();
		}
		
		return true;
	}
		



	private boolean isPulse (int i) {
		// some more assertion should come lat0r
		return (this.currMigrationMatrixList.get(i) == null);
	}


	public static class MsDemoFactory extends DemographyFactory {
		/*TODO
		 * - keep numbering of pop equal after merging. Done, not certain if correct
		 * - check time tol, working
		 * - specify times out of order, check splitting case
		 * - Check cases where a new mutable double should be created, and where the old one should be referenced
		 * */
		
		private int currentEpoch = 0;
		private double currentTime = 0.0f;
		private float timeTol = 1e-4f;
		
		TreeMap<Integer, MutableDouble> paramsToVary = new TreeMap<Integer, MutableDouble>();
		
		/*
		 * Parses the MS string and returns an MSParser object that has the 
		 * same structure as a Demography object
		 */
		public MsDemoFactory(String msString, double EPSILON, boolean halfMigrationRate, double[][] bounds) throws IOException {
			// super duper
			super(EPSILON, halfMigrationRate, bounds);
			
			String[] argString = msString.split("\\s+"); //split string input by spaces.
			
			ArrayList<String> sortedArg = sortByTime(argString);
			
			int i = 0;
			String curFlag = null;
			while(i<sortedArg.size()){
				try{
					curFlag = sortedArg.get(i);
				}catch(IndexOutOfBoundsException outBound){
					System.err.println("Index out of bounds when reading command line input.\n");
					System.exit(1);
				}
								
				try{
					if (curFlag.equals("ms")){
						//number of populations
						int sampleSize = Integer.parseInt(sortedArg.get(i+1));
						int numSamples = Integer.parseInt(sortedArg.get(i+2));
						i += 3;
					}else if(curFlag.equals("-I")){
						//Initialize epoch
						MutableDouble start = new MutableDouble(0.0);
						MutableDouble end = new MutableDouble(Double.POSITIVE_INFINITY);
						
						SimplePair<MutableDouble,MutableDouble> initPair = new SimplePair<MutableDouble,MutableDouble>(start,end);
						currEpochList.add(0, initPair);
						
						//initialize population
						int numPop = Integer.parseInt(sortedArg.get(i+1));

						MutableDouble[] initialPop = new MutableDouble[numPop];
						for (int j = 0; j<numPop; j++){
							initialPop[j] = getMutableDouble("1.0", paramsToVary);
						}
						currPopSizesList.add(0, initialPop);
						
						//Initialize treeset
						ArrayList<TreeSet<Integer>> initialTree = new ArrayList<TreeSet<Integer>>(numPop);
						for (int j = 0; j<numPop; j++){
							TreeSet<Integer> initialTreeSet = new TreeSet<Integer>(); 
							initialTreeSet.add(j+1); //each pop in its own set
							initialTree.add(initialTreeSet);
						}
						treePartitionList.add(initialTree);
						
						//Initialize Matrix
						MutableDouble[][] initialMatrix = new MutableDouble[numPop][numPop];
						for (int j=0; j<numPop; j++){
							for (int k=0; k<numPop; k++){
								initialMatrix[j][k] = getMutableDouble("0.0",paramsToVary);
							}
						}
						currMigrationMatrixList.add(0, initialMatrix);
										
						i += numPop + 2;
					}else if(curFlag.equals("-n")){
						//set size of population i to x*N, assumed at time 0
						int popToChange = Integer.parseInt(sortedArg.get(i+1)) - 1;//indexed from 0
						MutableDouble fraction = getMutableDouble(sortedArg.get(i+2), paramsToVary);
						
						MutableDouble[] curPop = currPopSizesList.get(0); //always updates initial epoch
						curPop[popToChange] = fraction;
						
						i += 3;
					}else if(curFlag.equals("-en")){
						
						float time = Float.parseFloat(sortedArg.get(i+1)); //set equal to right epoch
						boolean timeUpdated = updateCurrentState(time);
						int popToChange = Integer.parseInt(sortedArg.get(i+2)) - 1;
						MutableDouble fraction = getMutableDouble(sortedArg.get(i+3), paramsToVary);
						//fraction.setValue(-1.0 * fraction.getValue());
						
						MutableDouble[] curPop = currPopSizesList.get(currentEpoch); //update current epoch
						curPop[popToChange] = fraction;
												
																	
						i += 4;
					}else if(curFlag.equals("-m")){
						//set i,j element in migration matrix at time 0
						MutableDouble[][] curMat = currMigrationMatrixList.get(0);
						
						int x,y;
						x = Integer.parseInt(sortedArg.get(i+1)) - 1;
						y = Integer.parseInt(sortedArg.get(i+2)) - 1;
												
						curMat[x][y] = getMutableDouble(sortedArg.get(i+3), paramsToVary);
						
						i += 4;
					}else if(curFlag.equals("-em")){
						//set i,j element in migration matrix at current time
						float time = Float.parseFloat(sortedArg.get(i+1));
						boolean timeUpdated = updateCurrentState(time);
						
						MutableDouble[][] curMat = currMigrationMatrixList.get(currentEpoch);
						
						int x,y;
						x = Integer.parseInt(sortedArg.get(i+2)) - 1;
						y = Integer.parseInt(sortedArg.get(i+3)) - 1;
																	
						curMat[x][y] = getMutableDouble(sortedArg.get(i+4), paramsToVary); //new MutableDouble(val);
						
						
						i += 5;
					}else if(curFlag.equals("-ma")){
						
						//set all elements of migration matrix (row major order) at epoch 0. 
						MutableDouble[][] curMat = currMigrationMatrixList.get(0);
						int m = curMat.length; int n = curMat[0].length; //match initial matrix dimensions
						
						int j = i + 1;
						for (int x = 0; x<m; x++){
							for(int y = 0; y<n; y++){
								curMat[x][y] = getMutableDouble(sortedArg.get(j), paramsToVary);
								j++;
							}
						}
										
						i += m*n + 1;
					}else if(curFlag.equals("-ema")){
							//set all elements of migration matrix (row major order) at this epoch
							float time = Float.parseFloat(sortedArg.get(i+1));
							boolean timeUpdated = updateCurrentState(time);
							
							int numPop = Integer.parseInt(sortedArg.get(i+2));
							
							MutableDouble[][] curMat = new MutableDouble[numPop][numPop];
							
							int j = i + 3;
							for (int x = 0; x<numPop; x++){
								for(int y = 0; y<numPop; y++){
									curMat[x][y] = getMutableDouble(sortedArg.get(j), paramsToVary);
									j++;
								}
							}
							
							currMigrationMatrixList.remove(currentEpoch);
							currMigrationMatrixList.add(currentEpoch, curMat);
						
											
							i += numPop*numPop + 3;
					}else if(curFlag.equals("-ej")){
						//combine populations
						float time = Float.parseFloat(sortedArg.get(i+1));
						boolean timeUpdated = updateCurrentState(time);
						
						//move pop1 to pop2, leave pop2 size unchanged
						int pop1 = Integer.parseInt(sortedArg.get(i+2)) - 1;//indexed from 0
						int pop2 = Integer.parseInt(sortedArg.get(i+3)) - 1;
						
						ArrayList<TreeSet<Integer>> curTreeSet = treePartitionList.get(currentEpoch);
						curTreeSet.get(pop2).addAll(curTreeSet.get(pop1));
						curTreeSet.get(pop1).clear();
						
						
						MutableDouble[] curPopSize = currPopSizesList.get(currentEpoch);
						//curPopSize[pop2] += curPopSize[pop1]; //do or do not change pop size?
						MutableDouble[] newPopSize = new MutableDouble[curPopSize.length];
						for (int x =0; x<curPopSize.length; x++){
							if (x!=pop1){
								newPopSize[x] = curPopSize[x];
							}else{
								newPopSize[x] = null;
							}
						}
						currPopSizesList.remove(currentEpoch);
						currPopSizesList.add(currentEpoch, newPopSize);
						
						//update currMigrationMatrixList by removing row and column for pop1.
						MutableDouble[][] curMat = currMigrationMatrixList.get(currentEpoch);
						MutableDouble[][] newMat = new MutableDouble[newPopSize.length][newPopSize.length];
						for (int x = 0; x<curMat.length; x++){
							for (int y = 0; y<curMat.length; y++){
								if (x==pop1){
									newMat[x][y] = null; //getMutableDouble("0.0",paramsToVary);
								}else{
									newMat[x][y] = curMat[x][y];
								}
							}
						}
						currMigrationMatrixList.remove(currentEpoch);
						currMigrationMatrixList.add(currentEpoch, newMat);
						
						i += 4;
					}else if(curFlag.equals("-eg")){
						//duplicate current state but don't change anything
						float time = Float.parseFloat(sortedArg.get(i+1));
						boolean timeUpdated = updateCurrentState(time);
						
						i += 4;
					}else{
						if (!curFlag.substring(0,1).equals("-")){
							System.out.println("Flag \""+curFlag+"\" missing hyphen. Input prior might be bad. Continuing...\n");
						}else{
							System.out.println("Unrecognized command \""+curFlag+"\". Continuing...\n");
						}
						i++;
					}
				}catch(IndexOutOfBoundsException outBound){
					System.err.println("Error when reading flag \""+curFlag+"\" in constructor. Check numbers of arguments and population index.\n");
					System.exit(1);
				}catch(NumberFormatException nfe){
					System.err.println("Error when reading flag \""+curFlag+"\". Unexpected format after flag.\n");
					System.exit(1);
				}
			}
			
			// set currPoint
			int numParamsToVary = this.paramsToVary.size();
			this.currPoint = new MutableDouble[numParamsToVary];
			
			for (i = 0; i < numParamsToVary; i++) {
				if (!this.paramsToVary.containsKey(i)) throw new IOException("[ERROR] The wildcard numbers for the variable demographic parameters are not consecutive.");
				this.currPoint[i] = this.paramsToVary.get(i);
			}
			
			for (int j=0; j<currentEpoch+1; j++){
				this.currPulseMigrationMatrixList.add(null);
				
			}
			
			//clean lists
			this.clearNullElements();

		}
		
		/*
		 * Sort demography events by the time of event. Change to faster sort later
		 */
		private ArrayList<String> sortByTime(String[] argString){
			
			ArrayList<String> sortedArray = new ArrayList<String>();
			
			//find all time 0 flags
			int i = 0;
			while (i<argString.length){
				String x = argString[i];
				if (x.startsWith("-") && !x.substring(1,2).equals("e")){
					sortedArray.add(x);
					int j = i+1;
					while (j<argString.length && !argString[j].startsWith("-")){
						sortedArray.add(argString[j]);
						j++;
					}
					i = j;
				}else{
					i++;
				}
			}
			
			//now find other times
			i = 0;
			double minTime = 0d;
			double bestTime = Double.POSITIVE_INFINITY;
			String bestTimeString = "";
			
			while (sortedArray.size()<argString.length){
				//find the best time right now
				bestTime = Double.POSITIVE_INFINITY;
				i = 0;
				while (i<argString.length){
					String x = argString[i];
					if (x.startsWith("-e")){
						double time = Double.parseDouble(argString[i+1]);
						
						if (time>minTime && time<bestTime){
							bestTime = time;
							bestTimeString = argString[i+1];
						}
					}
					i++;
				}

				minTime = bestTime;
				
				//add all flags to sorted array that match bestTime
				i = 0;
				while (i<argString.length){
					String x = argString[i];
					//System.out.println(sortedArray);
					if (x.startsWith("-e") && argString[i+1].equals(bestTimeString)){
						sortedArray.add(x);
						int j = i+1;
						while(j<argString.length && !argString[j].startsWith("-")){
							sortedArray.add(argString[j]);
							j++;
						}
						i = j;
					}else{
						i++;
					}
				}
			}
			return sortedArray;
		}
		
		/*
		 * Switch currentEpoch corresponding to time. Or, if time is new, create a
		 * new epoch by duplicating the most immediate preceding epoch
		 */
		private boolean updateCurrentState(double time){
			if (time - currentTime> timeTol){
				//first make copy of previous state
				int previousEpoch = currentEpoch;
				currentEpoch++;
				currentTime = time;
				
				
				MutableDouble dtime = new MutableDouble(time);
								
				SimplePair<MutableDouble,MutableDouble> prevEpoch = currEpochList.get(previousEpoch);
				
				SimplePair<MutableDouble,MutableDouble> inv1 = new SimplePair<MutableDouble,MutableDouble>(prevEpoch.first(),dtime);
				SimplePair<MutableDouble,MutableDouble> inv2 = new SimplePair<MutableDouble,MutableDouble>(dtime,prevEpoch.second());
				
				currEpochList.remove(previousEpoch);
				currEpochList.add(previousEpoch, inv1);
				currEpochList.add(currentEpoch, inv2);
				
				MutableDouble[] curPop = this.currPopSizesList.get(previousEpoch);
				int numPop = curPop.length;
				
				//duplicate population sizes
				MutableDouble[] newPop = new MutableDouble[numPop];
				for (int i=0;i<numPop;i++){
					newPop[i] = curPop[i];
				}
				
				this.currPopSizesList.add(currentEpoch, newPop);
								
				//duplicate tree partition list
				ArrayList<TreeSet<Integer>> curTreeSet = treePartitionList.get(previousEpoch);
				ArrayList<TreeSet<Integer>> newTreeSet = new ArrayList<TreeSet<Integer>>();
				
				for (int j = 0; j<curTreeSet.size(); j++){
					newTreeSet.add (new TreeSet<Integer>(curTreeSet.get(j)));
				}
				
				this.treePartitionList.add(currentEpoch, newTreeSet);

				//duplicate migration matrix
				MutableDouble[][] curMatrix = this.currMigrationMatrixList.get(previousEpoch);
				MutableDouble[][] newMatrix = new MutableDouble[curMatrix.length][curMatrix.length];
				for (int j = 0; j<curMatrix.length; j++){
					for (int k =0; k<curMatrix.length; k++){
						newMatrix[j][k] = curMatrix[j][k]; //new MutableDouble(curMatrix[j][k].getValue());
					}
				}
				this.currMigrationMatrixList.add(currentEpoch, newMatrix);

				
				return true;
			}
			return false;
		}
		
		/*
		 * Resize all elements so that they match numPop at each epoch
		 */
		private void clearNullElements(){
			
			for (int i=0; i<currentEpoch+1; i++){
				//clear pop
				MutableDouble[] curPopSizes = currPopSizesList.get(i);
				
				ArrayList<MutableDouble> newPops = new ArrayList<MutableDouble>();
				int numNull = 0;
				for(int j=0; j<curPopSizes.length; j++){
					if (curPopSizes[j]==null){
						numNull++;
					}else{
						newPops.add(curPopSizes[j]);
					}
				}
				
				int newSize = newPops.size();
				
				//toArray doesn't work
				MutableDouble[] popArray = new MutableDouble[newPops.size()];
				for (int j=0; j<popArray.length; j++){
					popArray[j] = newPops.get(j);
				}
				
				currPopSizesList.remove(i);
				currPopSizesList.add(i, popArray);
				
				//clear treeset
				ArrayList<TreeSet<Integer>> curTreeSet = treePartitionList.get(i);
				Iterator<TreeSet<Integer>> iter = curTreeSet.iterator();
				TreeSet<Integer> curSet;
				ArrayList<TreeSet<Integer>> newTreeSet = new ArrayList<TreeSet<Integer>>();
				while(iter.hasNext()){
					curSet = iter.next();
					if ( !curSet.isEmpty() ){
						newTreeSet.add(curSet);
					}
				}
				
				treePartitionList.remove(i);
				treePartitionList.add(i, newTreeSet);
				
				//clear pop
				MutableDouble[][] curMatrix = currMigrationMatrixList.get(i);
				MutableDouble[][] newMatrix = new MutableDouble[newSize][newSize];
				
				int x = 0;
				for(int j=0; j<curMatrix.length; j++){
					int y = 0;
					if (curMatrix[j][0]==null){
						continue;
					}else{
						for (int k=0; k<curMatrix.length; k++){
							if (curMatrix[k][j]==null){
								continue;
							}else{
								newMatrix[x][y] = curMatrix[k][j];
								y++;
							}

						}
						x++;
					}
					
				}
				
				currMigrationMatrixList.remove(i);
				currMigrationMatrixList.add(i, newMatrix);

			}
		}
		
		private static MutableDouble getMutableDouble(String s, TreeMap<Integer, MutableDouble> paramsToVary){
			MutableDouble toReturn = null;
			if (s.startsWith("?")) {
				int currParamIndex = Integer.parseInt(s.substring(1));
				
				if (!paramsToVary.containsKey(currParamIndex)) {
					paramsToVary.put(currParamIndex, new MutableDouble(Double.NaN));
				}
				
				toReturn = paramsToVary.get(currParamIndex);
			} else {
				toReturn = new MutableDouble(Double.parseDouble(s));
			}
			
			return toReturn;
		}
		
		/*
		 * Print out the MSParser object in the Demography format.
		 */
		public String toString(){
			//print out this MSParser object
			String outString = "[";
			
			//first print the epochs
			int i = 0;
			for (i = 0; i<currEpochList.size()-1; i++){
				outString += currEpochList.get(i).second();
				if (i<currEpochList.size()-2){
					outString += ", ";
				}
			}
			outString += "]\n";
			
			//now print out current state for each epoch 
			for (i = 0; i<currEpochList.size(); i++){
				MutableDouble[] population = currPopSizesList.get(i);
				ArrayList<TreeSet<Integer>> tree = treePartitionList.get(i);
				MutableDouble[][] matrix = currMigrationMatrixList.get(i);
				
				//print population sets
				outString += "{";
				for (int j = 0; j<tree.size(); j++){
					outString += "{";
					
					TreeSet<Integer> curSet = tree.get(j);
					Iterator<Integer> itr = curSet.iterator();
					while (itr.hasNext()){
						int pop = itr.next();
						outString += pop;
						
						
						if (itr.hasNext()){
							outString += ",";
						}
					}
					outString += "}";
					if (j<tree.size()-1){
						outString += ",";
					}
				}
				outString += "}\n";
				
				//print population sizes
				for (int j = 0; j<population.length; j++){
					if (population[j]!=null){
						outString += population[j].getValue() + " ";
					}
				}
				outString += "\n";
				
				//Print pulse migration matrix
				outString += "null\n";
				
				//print migration matrix
				for (int j = 0; j<matrix.length; j++){
					if (population[j]==null){
						continue;
					}
					
					for (int k = 0; k<matrix.length; k++){
						if (matrix[k][j]!=null){
							outString += matrix[j][k].getValue() + " ";
						}
					}
					outString += "\n";
					
				}
			}
			
			return outString;
		}
	}
	
	
	
	public static class ParamFileDemoFactory extends DemographyFactory {
		
		public ParamFileDemoFactory (Reader demoReader, Reader expRatesReader, double EPSILON, boolean halfMigrationRate, double[][] bounds, boolean useEigenCore) throws IOException {
			super (EPSILON, halfMigrationRate, bounds);
						
			// read in the configuration first for some meta data
			BufferedReader bufferedReader = new BufferedReader(demoReader);
			
			TreeMap<Integer, MutableDouble> paramsToVary = new TreeMap<Integer, MutableDouble>();
			
			// go through lines
			boolean expectTimes = true;
			boolean expectPartition = false;
			boolean expectPopSizes = false;
			boolean expectPulseMatrix = false;
			boolean expectMigMatrix = false; 
			String line = null;
			int currMatrixLines = -1;
			int currNumPops = -1;
			MutableDouble[][] matrix = null;
			int numPresentDemes = -1;
			while ((line = bufferedReader.readLine()) != null) {
				// always ignore comment
				if (line.trim().startsWith("#") || line.trim().equals("")) continue;
				
				// what do we have
				if (expectTimes) {
					// get the times from the line
					String[] fields = line.trim().split("\\[");
					if (fields.length != 2 || !fields[0].equals("")) {
						throw new IOException("[ERROR] In demography file: Times are formatted incorrectly.");
					}
					fields = fields[1].split("\\]");
					if (fields.length != 1 && fields.length != 0) {
						throw new IOException("[ERROR] In demography file: Times are formatted incorrectly.");
					}
					
					// now take the middle
					if (fields.length > 0) fields = fields[0].split(",");
					
					ArrayList<MutableDouble> breakpoints = new ArrayList<MutableDouble>();
					for (String field : fields) {
												
						field = field.trim();
						
						if (field.equals("")) continue;
						
						breakpoints.add(getMutableDouble(field, paramsToVary));
					}
					
					// now make the epoch list
					MutableDouble lastEnd = new MutableDouble(0d);
					for (int i = 0; i < breakpoints.size(); i++) {
						MutableDouble currBreakpoint = breakpoints.get(i);
						
						if (currBreakpoint.getValue() < lastEnd.getValue()) throw new IOException("[ERROR] In demography file: Times are not in order.");
						
						this.currEpochList.add(new SimplePair<MutableDouble, MutableDouble>(lastEnd, currBreakpoint));
						lastEnd = currBreakpoint;
					}
					this.currEpochList.add(new SimplePair<MutableDouble, MutableDouble>(lastEnd, new MutableDouble(Double.POSITIVE_INFINITY)));
					
					
					// what to expect next
					expectTimes = false;
					expectPartition = true;
					expectPopSizes = false;
					expectPulseMatrix = false;
					expectMigMatrix = false;
				}
				else if (expectPartition) {

					line = line.trim();
					if (line.charAt(0) != '{' || line.charAt(line.length()-1) != '}'){
						throw new IOException("[ERROR] In demography file: Partition formatted incorrectly.");
					}
					
					String[] fields = line.substring(1, line.length()-1).split("(\t| )*\\{(\t| )*|(\t| )*\\}(\t| )*,(\t| )*\\{(\t| )*|(\t| )*\\}(\t| )*");
					if (!fields[0].equals("")) throw new IOException("[ERROR] In demography file: Partition formatted incorrectly.");
					
					ArrayList<TreeSet<Integer>> popList = new ArrayList<TreeSet<Integer>>();
					int numPops = 0;
					
					int j = 0;
					for (String field : fields) {
					
						if (j++ == 0) continue;
						
						field = field.trim();
						String[] subfields = field.split(",");

						TreeSet<Integer> ancientDeme = new TreeSet<Integer>();
						popList.add(ancientDeme);
						for (String subfield : subfields) {
							ancientDeme.add(Integer.parseInt(subfield.trim()));
							numPops++;
							
						}
					}

					
					// check numPops equals numPresentDemes
					if (this.treePartitionList.size() == 0) {
						numPresentDemes = numPops;
					} else if (numPresentDemes != numPops) {
						throw new IOException("[ERROR] In demography file: Number of present demes doesn't match number of populations at present.");
					} else {
						// check that popList is a refinement of the previous popList
						ArrayList<TreeSet<Integer>> prevPopList = this.treePartitionList.get(treePartitionList.size() - 1);
						
						// check some pop contains it
						for (TreeSet<Integer> childDeme : prevPopList) {
							int timesContained = 0;
							for (TreeSet<Integer> parentDeme : popList) {
								if (parentDeme.containsAll(childDeme)) timesContained++;
							}
							if (timesContained != 1) throw new IOException("[ERROR] In demography file: A partition in the list decribing the population tree is not a proper refinement.");
						}
						
					}
					
					// check popList contains each present deme exactly once
					TreeSet<Integer> combinedPop = new TreeSet<Integer>();
					for (TreeSet<Integer> ancientDeme : popList){
						combinedPop.addAll(ancientDeme);
					}
					for (int i = 0; i < numPops; i++) {
						if (!combinedPop.contains(i)) throw new IOException ("[ERROR] In demography file: A partition in the list decribing the population tree does not contain all present populations.");
					}
					
					// add it to the list
					this.treePartitionList.add(popList);
					
					currNumPops = popList.size();
					if (currNumPops <= 0) throw new IOException ("[ERROR] In demography file: A partition in the list decribing the population tree is empty.");
					
					// what to expect next
					expectTimes = false;
					expectPartition = false;
					expectPopSizes = true;
					expectPulseMatrix = false;
					expectMigMatrix = false;
				}
				else if (expectPopSizes) {
					
					
					String[] fields = line.trim().split("(\t| )+");
					if (fields.length != currNumPops) throw new IOException("[ERROR] In demography file: Number of population sizes given doesn't match with partition.");
					
					MutableDouble[] popSizes = new MutableDouble[currNumPops];
					this.currPopSizesList.add(popSizes);
					
					for (int i = 0; i < currNumPops; i++) {
						popSizes[i] = getMutableDouble(fields[i], paramsToVary);
						if (!Double.isNaN(popSizes[i].getValue()) && (popSizes[i].getValue() <= 0)) throw new IOException("[ERROR] In demography file: Negative population size given.");
					}
					
					// what to expect next
					expectTimes = false;
					expectPartition = false;
					expectPopSizes = false;
					expectPulseMatrix = true;
					expectMigMatrix = false;
					currMatrixLines = 0;
				}
				else if (expectPulseMatrix) {
					
					if (line.trim().equals("null")) {
						// no pulse matrix this time
						// what to expect next
						expectTimes = false;
						expectPartition = false;
						expectPopSizes = false;
						expectPulseMatrix = false;
						expectMigMatrix = true;
						currMatrixLines = 0;
						// we need it still
						matrix = null;
						// on to the next line
						continue;
					}
					
					assert (currMatrixLines >= 0);
					assert (currNumPops > 0);
					
					if (currMatrixLines == 0) {						
						matrix = new MutableDouble[currNumPops][currNumPops];
						this.currPulseMigrationMatrixList.add(matrix);
					}
					
					
					String[] fields = line.trim().split("(\t| )+");
					if (fields.length != currNumPops) throw new IOException ("[ERROR] In demography file: Pulse migration matrix is not a square matrix.");
					
					for (int i = 0; i < currNumPops; i++) {
						matrix[currMatrixLines][i] = getMutableDouble(fields[i], paramsToVary);
						if (currMatrixLines == i && (matrix[i][i].getValue() != 0 || Double.isNaN(matrix[i][i].getValue())) ) throw new IOException("[ERROR] In demography file: Diagonal values in pulse migration matrix are wrong.");
						else if (!Double.isNaN(matrix[currMatrixLines][i].getValue()) && (matrix[currMatrixLines][i].getValue() < 0)) throw new IOException("[ERROR] In demography file: Pulse migration matrix has negative entries.");
					}

					
					assert (currMatrixLines < currNumPops);
					if (currMatrixLines == currNumPops - 1) {
						// add infinitessimal interval before the last epoch
						int currEpochIdx = this.currMigrationMatrixList.size();
						// the epoch list is already full
						SimplePair<MutableDouble, MutableDouble> nextEpoch = this.currEpochList.get(currEpochIdx);
						SimplePair<MutableDouble, MutableDouble> currEpoch = new SimplePair<MutableDouble, MutableDouble>(nextEpoch.first(), nextEpoch.first());
						
						// ATTENTION, this shifts things
						this.currEpochList.add(currEpochIdx, currEpoch);
						
						// add the null for the migration matrix (we are at the right position, i think)
						this.currMigrationMatrixList.add(null);
						
						// add the pulse migration matrix
						// THIS IS ALREADY DONE
						
						// and the pop sizes and tree partitions have to be copied (shallow, cause should be the same mutable doubles)
						assert (this.currPopSizesList.size() - 1 == currEpochIdx);
						assert (this.treePartitionList.size()-1 == currEpochIdx);
						this.currPopSizesList.add(currEpochIdx, null);
						this.treePartitionList.add(this.treePartitionList.get(currEpochIdx));
						
						
						// what to expect next
						expectTimes = false;
						expectPartition = false;
						expectPopSizes = false;
						expectPulseMatrix = false;
						expectMigMatrix = true;
						currMatrixLines = 0;
						matrix = null;
					} else {
						currMatrixLines++;
					}
				}				
				else if (expectMigMatrix) {
					
					assert (currMatrixLines >= 0);
					assert (currNumPops > 0);
					
					if (currMatrixLines == 0) {						
						matrix = new MutableDouble[currNumPops][currNumPops];
						this.currMigrationMatrixList.add(matrix);
					}
					
					
					String[] fields = line.trim().split("(\t| )+");
					if (fields.length != currNumPops) throw new IOException("[ERROR] In demography file: Migration matrix not square.");
					
					for (int i = 0; i < currNumPops; i++) {
						matrix[currMatrixLines][i] = getMutableDouble (fields[i], paramsToVary);
						if (currMatrixLines == i && (matrix[i][i].getValue() > 0 || Double.isNaN(matrix[i][i].getValue())) ) throw new IOException("[ERROR] In demography file: Entry on diagonal of migration matrix greater than zero.");
						else if (!Double.isNaN(matrix[currMatrixLines][i].getValue()) && (matrix[currMatrixLines][i].getValue() < 0)) throw new IOException("[ERROR] In demography file: Off-diagonal entry of migration matrix smaller than zero.");
					}

					
					assert (currMatrixLines < currNumPops);
					if (currMatrixLines == currNumPops - 1) {
						// add an empty pulse matrix
						this.currPulseMigrationMatrixList.add(null);
						// what to expect next
						expectTimes = false;
						expectPartition = true;
						expectPopSizes = false;
						expectMigMatrix = false;
						expectPulseMatrix = false;
						currMatrixLines = -1;
						currNumPops = -1;
						matrix = null;
					} else {
						currMatrixLines++;
					}
				}				
			}
			// done reading demography
			
			// check it			
			int numIntervals = this.currEpochList.size();
			if (this.currMigrationMatrixList.size() != numIntervals) throw new IOException("[ERROR] In demography file: Not enough migration matrices.");
			if (this.treePartitionList.size() != numIntervals) throw new IOException("[ERROR] In demography file: Not enough partitions.");
			if (this.currPopSizesList.size() != numIntervals) throw new IOException("[ERROR] In demography file: Not enough population size sets.");
			if (this.currPulseMigrationMatrixList.size() != numIntervals) throw new IOException("[ERROR] In demography file: Not enough pulse migration matrices.");
			if (this.treePartitionList.get(this.treePartitionList.size()-1).size() > 1) {
				System.out.println("# [WARNING] In demography file: The last epoch has more than one deme. This could lead to long runtimes.");
			}
			
			// now maybe read exponential rates
			if (expRatesReader != null) {
				// read some rates
				this.currExpRatesList = new ArrayList<MutableDouble[]>();
				
				// if there are '?'-wildcards here, they have to be bigger then the ones in the demography file for now (if this is actually ever checked)
				// open reader (again)
				BufferedReader bufferedRatesReader = new BufferedReader (expRatesReader);
				
				// we need an index because of pulses
				int intervalIdx = 0;
				// and go through lines
				while ((line = bufferedRatesReader.readLine()) != null) {
					// always ignore comment
					if (line.trim().startsWith("#") || line.trim().equals("")) continue;

					// now it should be a comma separated lines of doubles (we check the size later)
					String[] fields = line.trim().split("(\t| )+");
					
					// container for the rates
					MutableDouble[] rates = new MutableDouble[fields.length];
					// pulse gets nothing
					if (this.currPopSizesList.get(intervalIdx) == null) {
						this.currExpRatesList.add (null);
						// and we jump by one
						intervalIdx++;
					}
					this.currExpRatesList.add(rates);
					
					// now read 'em
					for (int i = 0; i < rates.length; i++) {
						rates[i] = getMutableDouble(fields[i], paramsToVary);
						// can be arbitrary double, no restrictions here
					}
					
					//up with the index
					intervalIdx++;
				}
				
				// and make sure that they match with the already read demography
				if (this.currExpRatesList.size() != numIntervals) throw new IOException();
				for (int i=0; i<numIntervals; i++) {
					// account for pulse epochs
					if ((this.currPopSizesList.get(i) == null) != (this.currExpRatesList.get(i) == null)) throw new IOException("[ERROR] In demography file: Each epoch needs population sizes or exponential rates (at least).");
					if (this.currExpRatesList.get(i) != null) {
						if (this.currExpRatesList.get(i).length != this.currPopSizesList.get(i).length) throw new IOException("[ERROR] In demography file: Number of exponential rates doesn't match number of populations.");
					}
				}
			}
			else {
				// no rates wanted
				
				if (useEigenCore) {
					// really none
					this.currExpRatesList = null;
				}
				else {
					// make up dummy rates
					this.currExpRatesList = new ArrayList<MutableDouble[]>();
					for (int i=0; i<numIntervals; i++) {
						// account for pulse epochs
						if (this.currPopSizesList.get(i) == null) {
							this.currExpRatesList.add(null);
						}
						else {
							MutableDouble[] themRates = new MutableDouble[this.currPopSizesList.get(i).length];
							
							for (int l = 0; l < themRates.length; l++) {
								// zero rates
								themRates[l] = new MutableDouble(0.0);
							}
					
							this.currExpRatesList.add(themRates);
						}
					}
				}

			}
			
			
			// set currPoint
			int numParamsToVary = paramsToVary.size();
			if ((bounds != null) && (numParamsToVary != bounds.length)) {
				throw new IOException ("Number of Parameters to vary doesn't match number of bounds given.");
			}
			this.currPoint = new MutableDouble[numParamsToVary];
			
			for (int i = 0; i < numParamsToVary; i++) {
				if (!paramsToVary.containsKey(i)) throw new IOException ("Parameter numbering in demography-file not contiguous.");
				this.currPoint[i] = paramsToVary.get(i);
			}			
			
		}
	}
	
	private static MutableDouble getMutableDouble (String s, TreeMap<Integer, MutableDouble> paramsToVary){
		MutableDouble toReturn = null;
		if (s.startsWith("?")) {
			int currParamIndex = Integer.parseInt(s.substring(1));
			
			if (!paramsToVary.containsKey(currParamIndex)) {
				paramsToVary.put(currParamIndex, new MutableDouble(Double.NaN));
			}
			
			toReturn = paramsToVary.get(currParamIndex);
		} else {
			toReturn = new MutableDouble(Double.parseDouble(s));
		}
		
		return toReturn;
	}

	public double getMinNonZeroMigrationRate() {
		double minRate = Double.POSITIVE_INFINITY;
		// go through matrices and get min
		// if we estimate parameters, they should be NaN
		for (MutableDouble[][] currMatrix : this.currMigrationMatrixList) {
			// deal with zero length epochs
			if (currMatrix == null) continue;
			for (int i=0; i<currMatrix.length; i++) {
				for (int j=0; j<currMatrix[i].length; j++) {
					// no diagonal
					if (i != j) {
						double rate = currMatrix[i][j].getValue();
						if ((rate > this.EPSILON) && (rate != Double.NaN)) {
							minRate = Math.min (minRate, rate);
						}
					}
				}
			}
		}
		// return it
		return minRate;
	}

	public void addMinimalMigration(double minimalRate) {
		// go through matrices and add minimal rates
		// the estimation values should be NaN
		for (MutableDouble[][] currMatrix : this.currMigrationMatrixList) {
			// again, ignore epochs of length zero
			if (currMatrix == null) continue;
			for (int i=0; i<currMatrix.length; i++) {
				for (int j=0; j<currMatrix[i].length; j++) {
					// no diagonal
					if (i != j) {
						double rate = currMatrix[i][j].getValue();
						assert (rate >= 0d || Double.isNaN(rate));
						if ((rate < minimalRate) && (rate != Double.NaN)) {
							currMatrix[i][j].setValue(minimalRate);
						}
					}
				}
			}
		}
	}

	public double[][] getBounds() {
		return this.bounds;
	}

	private String realValue (MutableDouble d) {
		for (int i=0; i<this.currPoint.length; i++) {
			// this should work
			if (this.currPoint[i] == d) {
				return "?" + i;
			}
		}
		return Double.toString (d.getValue());
	}
	
	public void dump (PrintStream outStream, boolean addCommentChar) {
		assert (addCommentChar == true);
		// go through the epochs
		for (int e=0; e<this.currEpochList.size(); e++) {
			outStream.println ("# ==== Epoch " + e + " ====");
			// say the interval
			outStream.println ("# " + this.realValue(this.currEpochList.get(e).first()) + "\t" + this.realValue(this.currEpochList.get(e).second()));
			// say the partition
			outStream.print("# {");
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
			if (this.currPopSizesList.get(e) == null) {
				outStream.print("# null");
			} else {
				outStream.print("# ");
				for (MutableDouble currSize : this.currPopSizesList.get(e)) {
					outStream.print (this.realValue(currSize) + "\t");
				}
			}
			outStream.println();
			
			// maybe say the exp rates, if there are some
			outStream.print ("# rates:\t");
			if (this.currExpRatesList == null || (this.currExpRatesList.get(e) == null)) {
				outStream.print("null");
			} else {
				for (MutableDouble currRate : this.currExpRatesList.get(e)) {
					outStream.print (this.realValue(currRate) + "\t");
				}
			}
			outStream.println();
			
			// say the pulse matrix
			if (this.currPulseMigrationMatrixList.get(e) == null) {
				outStream.println("# null");
			} else {
				for (MutableDouble[] currRow : this.currPulseMigrationMatrixList.get(e)) {
					outStream.print("# ");
					for (MutableDouble currVal : currRow) {
						outStream.print (this.realValue(currVal) + "\t");
					}
					outStream.println();
				}
			}
			
			// say the mig matrix
			if (this.currMigrationMatrixList.get(e) == null) {
				outStream.println("# null");
			} else {
				for (MutableDouble[] currRow : this.currMigrationMatrixList.get(e)) {
					outStream.print("# ");
					for (MutableDouble currVal : currRow) {
						outStream.print (this.realValue(currVal) + "\t");
					}
					outStream.println();
				}
			}

		}
	}
}
