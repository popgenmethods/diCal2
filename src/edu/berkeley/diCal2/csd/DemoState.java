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
import java.util.List;

import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.haplotype.FSAHaplotype;
import edu.berkeley.diCal2.utility.SimplePair;
import edu.berkeley.diCal2.utility.RealPartition.Interval;

// a DemoState is a collection of (epoch, ancientDeme), along with a single presentDeme
// the hidden states of the diCal HMM consist of pairs of (DemoState, haplotype), for haplotypes residing in the DemoState's presentDeme
public class DemoState {
	
	public int presentDeme;
	List<SimplePair<Integer,Integer>> EpochAncientDemePairs = new ArrayList<SimplePair<Integer,Integer>>(); // a list of (epoch, ancientDeme) residing in this demoState
	
	// a class to make DemoStateCollections
	public static class DemoStateFactory {
		final IntervalFactory intervalFactory;
		final int additionalCakeIntervals;
		final IntervalFactory additionalCakeIntervalFactory;
		final double EPSILON;
		
		// use the old cake intervals?		
		public static boolean OLD_CAKE_INTERVALS = false;	
		
		public DemoStateFactory (IntervalFactory intervalFactory, int additionalCakeIntervals, double EPSILON) {
			super();
			this.intervalFactory = intervalFactory;
			this.additionalCakeIntervals = additionalCakeIntervals;
			this.EPSILON = EPSILON;
			if (!DemoStateFactory.OLD_CAKE_INTERVALS){
				this.additionalCakeIntervalFactory = additionalCakeIntervals == 0 ? null : new IntervalFactory.MeanRateIntervalFactory (additionalCakeIntervals+1, false);
			} else {
				this.additionalCakeIntervalFactory = additionalCakeIntervals == 0 ? null : new IntervalFactory.OldCakeIntervalFactory (additionalCakeIntervals, false);
			}
		}
		
		public synchronized <H extends FSAHaplotype> DemoStateCollection getDemoStates(Demography demo, DemoConfiguration<H> config, int additionalDeme) {
			if (config.getTotalGeneticTypes() == 0) return null;
			
			Interval[] additionalCakeIntervals = this.additionalCakeIntervalFactory == null ? null : this.additionalCakeIntervalFactory.getIntervals (demo, config, additionalDeme);
			if (intervalFactory == null) {
				return new DemoStateCollection(demo, additionalCakeIntervals, this.EPSILON);
			} else {
				Interval[] intervals = this.intervalFactory.getIntervals (demo, config, additionalDeme);
				return new DemoStateCollection(demo, intervals, additionalCakeIntervals, this.EPSILON);
			}
		}
	}
	
	
	// an object to keep track of all DemoStates
	// WARNING: when creating DemoStateCollections several times over an EM run, you MUST make sure the indices of the DemoState's match up!
	public static class DemoStateCollection {
		
		public enum DemoStateType {
			IntervalPresentDeme, IntervalPresentAncientDeme
		}

		private final int[] numIntervalsInHiddenTime;
		
		public final DemoStateType type;
 
		public final Demography refinedDemography;
		
		// [refinedInterval][ancientDeme][presentDeme]
		// the demoState for each (refinedInterval, ancientDeme, presentDeme) triplet
		private final List<List<int[]>> demoStateMap;
		
		// a list of the demoStates
		private final List<DemoState> demoStateList;
		
		// start and end time for each demoState
		private final double[] demoStateStartTimes;
		private final double[] demoStateEndTimes;
		
	
		public int getNumIntervalsInHiddenTime(int hiddenTimeIdx){
			return this.numIntervalsInHiddenTime[hiddenTimeIdx];
		}
		
		public int getNumHiddenTimes(){
			return this.numIntervalsInHiddenTime.length;
		}
		
		
		public int getDemoState(int epoch, int ancientDeme, int presentDeme) {
			return this.demoStateMap.get(epoch).get(ancientDeme)[presentDeme];
		}
		
		public int getPresentDeme (int demoStateIdx) {
			assert (demoStateIdx >= 0 && demoStateIdx < this.numDemoStates());
			return this.demoStateList.get(demoStateIdx).presentDeme;
		}
		
		public <H extends FSAHaplotype> HaplotypeConfiguration<H> getPresentHapConfig (DemoConfiguration<H> config, int intervalDemeIdx) {
			assert (intervalDemeIdx >= 0 && intervalDemeIdx < this.numDemoStates());
			
			return config.getPopulation(this.getPresentDeme(intervalDemeIdx));
		}
		
		public int numDemoStates () {
			return this.demoStateList.size();
		}
		
		public double startTime(int demoState) {
			return demoStateStartTimes[demoState];
		}
		
		public double endTime(int demoState) {
			return demoStateEndTimes[demoState];
		}
		
		
		
		// CONSTRUCTORS
		
		public static void initializeDemoStateMap(Demography demo, List<List<int[]>> demoStateMap) {
			int numEpochs = demo.epochList.length;
			int numPresentDemes = demo.numAncientDemes(0);
			for (int epoch = 0; epoch < numEpochs; epoch++) {
				if (demo.isPulse(epoch)) {
					demoStateMap.add(null);
					continue;
				}
				
				demoStateMap.add(new ArrayList<int[]>());
				int numAncientDemes = demo.popSizes.get(epoch).length;
				for (int ancientDeme = 0; ancientDeme < numAncientDemes; ancientDeme++) {
					demoStateMap.get(epoch).add(new int[numPresentDemes]);
				}
			}
		}
		
		private static void initializeDemoStateTimes(Demography refinedDemo, List<DemoState> demoStateList, double[] startTimes, double[] endTimes) {
			assert(startTimes.length == endTimes.length);
			assert(startTimes.length == demoStateList.size());
			
			for (int demoStateIdx = 0; demoStateIdx < demoStateList.size(); demoStateIdx++) {
				double start = Double.POSITIVE_INFINITY;
				double end = Double.NEGATIVE_INFINITY;
				
				DemoState demoState = demoStateList.get(demoStateIdx);
				for (SimplePair<Integer,Integer> epochAncientDeme : demoState.EpochAncientDemePairs) {
					Interval epoch = refinedDemo.epochList[epochAncientDeme.first()];
					start = Math.min(start, epoch.startPoint);
					end = Math.max(end, epoch.endPoint);
				}
				
				assert (start >= 0d);
				assert (start < Double.POSITIVE_INFINITY);
				assert (end > 0d);

				startTimes[demoStateIdx] = start;
				endTimes[demoStateIdx] = end;
			}
		}
		
		
		// this constructs a demoState for each epoch-ancientDeme-presentDeme triplet
		public DemoStateCollection(Demography demo, Interval[] additionalCakeIntervals, double EPSILON) {
			//WARNING: This Constructor breaks with the whole multiple theta thing
			// this is not used at the moment
			assert(false);
			
			
			this.numIntervalsInHiddenTime= new int[demo.epochList.length];
			
			// set the type
			this.type = DemoStateType.IntervalPresentAncientDeme;
			
			// get the original epochs
			Interval[] origEpochs = demo.epochList;
			
			// now map the refined intervals back to the original epochs
			ArrayList<Integer> refinedIntervalsToEpochs = new ArrayList<Integer>();
			this.refinedDemography = new Demography(demo, additionalCakeIntervals, null, refinedIntervalsToEpochs, EPSILON);
			
			// initialize the demoState variables
			this.demoStateList = new ArrayList<DemoState>();
			this.demoStateMap = new ArrayList<List<int[]>>();
			
			// get the refined intervals by epoch
			List<List<Integer>> epochToRefinedIntervals = new ArrayList<List<Integer>>();
			for (Interval epoch : origEpochs) {
				epochToRefinedIntervals.add(new ArrayList<Integer>());
			}
			for (int refinedInterval = 0; refinedInterval < this.refinedDemography.epochList.length; refinedInterval++) {
				int epoch = refinedIntervalsToEpochs.get(refinedInterval);
				epochToRefinedIntervals.get(epoch).add(refinedInterval);
				this.numIntervalsInHiddenTime[epoch]++;
			}
			
			// create storage for demoStateMap
			initializeDemoStateMap(this.refinedDemography, this.demoStateMap);
			
			// create demoStates by going through each (epoch,ancientDeme,presentDeme) triplet
			for (int epoch = 0; epoch < origEpochs.length; epoch++) {
				
				if (demo.isPulse(epoch)) {
					for (int refinedInterval : epochToRefinedIntervals.get(epoch)) {
						assert (this.refinedDemography.isPulse(refinedInterval));
					}
					continue;
				}
				
				for (int ancientDeme = 0; ancientDeme < demo.popSizes.get(epoch).length; ancientDeme++) {
					for (int presentDeme = 0; presentDeme < demo.popSizes.get(0).length; presentDeme++) {
						// make a DemoState for this (epoch,ancientDeme,presentDeme)
						DemoState currDemoState = new DemoState();
						
						// store it in the demoStateList and get its index
						int currDemoStateIdx = demoStateList.size();
						this.demoStateList.add(currDemoState);
						
						// set the presentDeme of the DemoState
						currDemoState.presentDeme = presentDeme;

						// go thru the refined intervals in the epoch, and add them
						for (int refinedInterval : epochToRefinedIntervals.get(epoch)) {
							
							assert (!this.refinedDemography.isPulse(refinedInterval));
							
							// add the epoch-ancientDeme pair to the DemoState
							currDemoState.EpochAncientDemePairs.add(new SimplePair<Integer,Integer>(refinedInterval, ancientDeme));
							
							// add the DemoState to the map
							this.demoStateMap.get(refinedInterval).get(ancientDeme)[presentDeme] = currDemoStateIdx;
						}
						
					}
				}
			}
			
			this.demoStateStartTimes = new double[this.numDemoStates()];
			this.demoStateEndTimes = new double[this.numDemoStates()];
			initializeDemoStateTimes(refinedDemography, demoStateList, demoStateStartTimes, demoStateEndTimes);
			
		}
		
		// this constructs a demoState for each epoch-presentDeme pair
		// so all ancientDemes within an epoch lie in the same DemoState
		public DemoStateCollection (Demography origDemo, Interval[] hiddenStateIntervals, Interval[] additionalCakeIntervals, double EPSILON) {
			
			this.numIntervalsInHiddenTime= new int[hiddenStateIntervals.length];
			
			// set the type
			this.type = DemoStateType.IntervalPresentDeme;
			

			assert (hiddenStateIntervals != null);
			
			int numStateIntervals = hiddenStateIntervals.length;
			ArrayList<Integer> epochToStateInterval = new ArrayList<Integer>();
			
			// initialize private variables
			origDemo = additionalCakeIntervals == null ? origDemo : new Demography(origDemo, additionalCakeIntervals, new ArrayList<Integer>(), EPSILON);
			this.refinedDemography = new Demography (origDemo, hiddenStateIntervals, epochToStateInterval, EPSILON);
			this.demoStateList = new ArrayList<DemoState>();
			this.demoStateMap = new ArrayList<List<int[]>>();
			
			// create storage for demoStateMap
			initializeDemoStateMap(this.refinedDemography, this.demoStateMap);
			
			// for each state interval, create a list of epochs that belong to it
			List<List<Integer>> epochsByStateInterval = new ArrayList<List<Integer>>();
			for (int stateInterval = 0; stateInterval < numStateIntervals; stateInterval++) {
				epochsByStateInterval.add(new ArrayList<Integer>());
			}
			// fill epochsByStateInterval
			for (int epoch = 0; epoch < this.refinedDemography.epochList.length; epoch++) {
				int stateInterval = epochToStateInterval.get(epoch);
				epochsByStateInterval.get(stateInterval).add(epoch);
				this.numIntervalsInHiddenTime[stateInterval]++;
			}
			
			// now create the demoStates
			for (int stateInterval = 0; stateInterval < numStateIntervals; stateInterval++) {			
				for (int presentDeme = 0; presentDeme < this.refinedDemography.numAncientDemes(0); presentDeme++) {
					// make a DemoState for this stateInterval,presentDeme
					DemoState currDemoState = new DemoState();
					
					int currDemoStateIdx = demoStateList.size();
					this.demoStateList.add(currDemoState);
					
					// set the presentDeme of the DemoState
					currDemoState.presentDeme = presentDeme;
					
					// add all epochs that belong in the state interval
					for (int epoch : epochsByStateInterval.get(stateInterval)) {
						
						if (this.refinedDemography.isPulse(epoch)) {
							continue;
						}
						
						for (int ancientDeme = 0; ancientDeme < this.refinedDemography.popSizes.get(epoch).length; ancientDeme++) {
							currDemoState.EpochAncientDemePairs.add(new SimplePair<Integer,Integer>(epoch, ancientDeme));
							this.demoStateMap.get(epoch).get(ancientDeme)[presentDeme] = currDemoStateIdx;
						}
					}
				}
			}
			// done with constructor
			
			this.demoStateStartTimes = new double[this.numDemoStates()];
			this.demoStateEndTimes = new double[this.numDemoStates()];
			initializeDemoStateTimes(refinedDemography, demoStateList, demoStateStartTimes, demoStateEndTimes);
		}
	}
}
	
