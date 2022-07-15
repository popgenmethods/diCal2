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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.berkeley.diCal2.csd.auxiliary.LaGuerre;
import edu.berkeley.diCal2.csd.auxiliary.LaGuerre.QuadraturePoint;
import edu.berkeley.diCal2.demography.ConstDemography;
import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.haplotype.GeneticType;
import edu.berkeley.diCal2.maximum_likelihood.MetaOptimization;
import edu.berkeley.diCal2.utility.RealPartition;
import edu.berkeley.diCal2.utility.RealPartition.Interval;

public abstract class IntervalFactory {
	
	private final boolean printIntervals;

	public IntervalFactory (boolean printIntervals) {
		this.printIntervals = printIntervals;
	}

	protected abstract <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract (Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme);

	// wrapper
	public <H extends GeneticType> Interval[] getNontrivialIntervals (Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
		
		// get intervals from the internal function
		Interval[] thePartition = this.getNontrivialIntervalsAbstract (demo, presentConfig, additionalDeme);
		
		// debug, if wanted
		if (this.printIntervals) {
			String debugString = "# theIntervals:\t" + Arrays.toString (thePartition);
			MetaOptimization.synchronizedPrintln (debugString);
		}
		
		// and then return them
		return thePartition;
	}
	
	
	public abstract String factoryType();
	public abstract String factoryParams();
	
	
	
	public static IntervalFactory getIntervalFactory (String factoryName, String factoryParams, boolean printIntervals) throws IOException {
		String[] intervalFactoryParams = factoryParams.split(",");
		
		IntervalFactory intervalFactory = null;
		
		if (factoryName.equals("single")) {

			intervalFactory = new IntervalFactory.SingleIntervalFactory (printIntervals);
			
			if (intervalFactoryParams.length != 1 || !intervalFactoryParams[0].equals("")){
				System.err.println("# Warning: single interval method does not require parameters. Ignoring inputted params.");
			}
		} else if (factoryName.equals("simple")) {
			
			intervalFactory = new IntervalFactory.SimpleIntervalFactory (printIntervals);
			
			if (intervalFactoryParams.length != 1 || !intervalFactoryParams[0].equals("")){
				System.err.println("# Warning: simple interval method does not require parameters. Ignoring inputted params.");
			}
		} else if (factoryName.equals("meanrate")) {
			
			
			boolean validParams = intervalFactoryParams.length == 1;
			if (validParams) {
				try {
					int numIntervals = Integer.parseInt(intervalFactoryParams[0]);
					if (numIntervals <= 0) validParams = false;
					intervalFactory = new IntervalFactory.MeanRateIntervalFactory (numIntervals, printIntervals);
				} catch (NumberFormatException ex) {
					validParams = false;
				}
			}
			
			if (!validParams) {
				System.err.println("# Invalid parameters for mean rate interval method.");
				return null;
			}
			
		} else if (factoryName.equals("mixtureuniform")) {
			
			boolean validParams = intervalFactoryParams.length == 1;
			if (validParams) {
				try {
					int numIntervals = Integer.parseInt(intervalFactoryParams[0]);
					if (numIntervals <= 0) validParams = false;
					intervalFactory = new IntervalFactory.MixtureUniformIntervalFactory (numIntervals, UberDemographyCore.ONELOCUS_EPSILON, printIntervals);
				} catch (NumberFormatException ex) {
					validParams = false;
				}
			}
			
			if (!validParams) {
				System.err.println("# Invalid parameters for uniform mixture interval method.");
				return null;
			}
		} else if (factoryName.equals("old")) {
			
			boolean validParams = intervalFactoryParams.length == 1;
			if (validParams) {
				try {
					int numIntervals = Integer.parseInt(intervalFactoryParams[0]);
					if (numIntervals <= 0) validParams = false;
					intervalFactory = new IntervalFactory.OldIntervalFactory (numIntervals, printIntervals);
				} catch (NumberFormatException ex) {
					validParams = false;
				}
			}
			
			if (!validParams) {
				System.err.println("# Invalid parameters for old interval method.");
				return null;
			}
		} else if (factoryName.equals("fixeduniform")) {
			
			boolean validParams = intervalFactoryParams.length == 2;
			if (validParams) {
				try {
					int numIntervals = Integer.parseInt(intervalFactoryParams[0]);
					double effPopSize = Double.parseDouble(intervalFactoryParams[1]);
					
					if (numIntervals <= 0 || effPopSize <= 0) validParams = false;
					intervalFactory = new IntervalFactory.FixedUniformIntervalFactory (numIntervals, effPopSize, printIntervals);
				} catch (NumberFormatException ex) {
					validParams = false;
				}
			}
			
			if (!validParams) {
				System.err.println("# Invalid parameters for fixedUniform interval method.");
				return null;
			}
			
		} else if (factoryName.equals("quadrature")) {
			
			boolean validParams = intervalFactoryParams.length == 2;
			if (validParams) {
				try {
					int numIntervals = Integer.parseInt(intervalFactoryParams[0]);
					double rescaling = Double.parseDouble(intervalFactoryParams[1]);
					
					if (numIntervals <= 0 || rescaling <= 0d) validParams = false;
					intervalFactory = new IntervalFactory.QuadratureFactory (numIntervals, rescaling, printIntervals);
				} catch (NumberFormatException ex) {
					validParams = false;
				}
			}
			
			if (!validParams) {
				System.err.println("# Invalid parameters for quadrature interval method.");
				return null;
			}
		} else if (factoryName.equals("quadrature2")) {
			
			boolean validParams = intervalFactoryParams.length == 2;
			if (validParams) {
				try {
					int numIntervals = Integer.parseInt(intervalFactoryParams[0]);
					double rescaling = Double.parseDouble(intervalFactoryParams[1]);
					
					if (numIntervals <= 0 || rescaling <= 0d) validParams = false;
					intervalFactory = new IntervalFactory.QuadratureFactory2 (numIntervals, rescaling, printIntervals);
				} catch (NumberFormatException ex) {
					validParams = false;
				}
			}
			
			if (!validParams) {
				System.err.println("# Invalid parameters for quadrature interval method.");
				return null;
			}
		} else if (factoryName.equals("custom") || factoryName.equals("customfixed")) {
			
			boolean validParams = intervalFactoryParams.length >= 1;
			
			if (validParams) {
				try {
					
					double[] points = new double[intervalFactoryParams.length];
					
					double lastPoint = 0d;
					for (int i = 0; i < intervalFactoryParams.length; i++) {
						double nextPoint = Double.parseDouble(intervalFactoryParams[i]);
						if (nextPoint <= lastPoint) validParams = false;
						points[i] = nextPoint;
						lastPoint = nextPoint;
					}
					if (lastPoint == Double.POSITIVE_INFINITY) validParams = false;
							
					if (factoryName.equals("custom")) {
						intervalFactory = new IntervalFactory.CustomIntervalFactory (points, printIntervals);
					} else if (factoryName.equals("customfixed")) {
						intervalFactory = new IntervalFactory.CustomFixedIntervalFactory (points, printIntervals);
					} else {
						assert (false);
						return null;
					}
				} catch (NumberFormatException ex) {
					validParams = false;
				}
			}
			
			if (!validParams) {
				System.err.println("# Invalid parameters for custom interval method.");
				return null;
			}

		} else if (factoryName.equals("customfixedFile")) {
			FileReader reader = new FileReader(factoryParams);
			BufferedReader bufferedReader = new BufferedReader(reader);
			String line;
			List<String> lines = new ArrayList<String>();
			while ((line = bufferedReader.readLine()) != null) {
				// always ignore comment
				if (line.trim().startsWith("#") || line.trim().equals("")) continue;
				lines.add(line);
			}
			bufferedReader.close();
			
			if (lines.size() != 1) {
				System.err.println("# Invalid file format for customfixed File");
				return null;
			}
			line = lines.get(0);
			String[] valStrings = line.split(",");
			double[] points = new double[valStrings.length];
			for (int i = 0; i < points.length; i++) points[i] = Double.parseDouble(valStrings[i]);
			
			intervalFactory = new IntervalFactory.CustomFixedIntervalFactory (points, printIntervals);
			
		} else if (factoryName.equals("combined")) {
			String[] subfactoryParams = factoryParams.split("/");
			boolean validParams = subfactoryParams.length > 0;
			
			List<IntervalFactory> memberFactories = new ArrayList<IntervalFactory>();
			
			for (String subfacPrm : subfactoryParams) {
				String[] nameAndParams = subfacPrm.split(":");
				if (nameAndParams.length != 2) {
					validParams = false;
					break;
				}
				
				String name = nameAndParams[0];
				String params = nameAndParams[1];
				
				memberFactories.add(getIntervalFactory (name, params, printIntervals));
			}
			intervalFactory = new CombinedIntervalFactory (memberFactories, UberDemographyCore.ONELOCUS_EPSILON, printIntervals);
			
			if (!validParams) {
				System.err.println("# Invalid parameters for combined interval method.");
				return null;
			}

		} else if (factoryName.equals("loguniform")) {

			boolean validParams = (intervalFactoryParams.length % 2  == 1);
			
			if (validParams) {
				try {
					// get the parameters (this should maybe go into the constructor of the interval factory, but then again so should splitting up the params by comma)
					List<Integer> intervalNumbers = new ArrayList<Integer>();
					List<Double> intervalBorders = new ArrayList<Double>();
					for (int i=0; i<(intervalFactoryParams.length/2); i++) {
						intervalNumbers.add(Integer.parseInt(intervalFactoryParams[i]));
					}
					for (int i=(intervalFactoryParams.length/2); i<intervalFactoryParams.length; i++) {
						intervalBorders.add(Double.parseDouble(intervalFactoryParams[i]));
					}
					assert ((intervalNumbers.size()+1) == intervalBorders.size());
					
					// see whether all is good
					boolean good = true;
					for (int i=0; i<intervalNumbers.size(); i++) {
						good &= (intervalNumbers.get(i) > 0);
						good &= (intervalBorders.get(i) < intervalBorders.get(i+1));
					}
					
					if (good) {
						intervalFactory = new IntervalFactory.LogUniformFactory (intervalNumbers, intervalBorders, printIntervals);
					}
					else {
						validParams = false;
					}
				} catch (NumberFormatException ex) {
					validParams = false;
				}
			}
			
			// everything ok? 
			if (!validParams) {
				throw new IOException ("# Invalid parameters for loguniform interval method. Have to give n integers and n+1 boundaries");
			}
		} else if (factoryName.equals("exponentialGrowth")) {
			intervalFactory = new ExponentialGrowthIntervalFactory(Double.parseDouble(factoryParams), printIntervals);
		} else if (factoryName.equals("fixedGrid")) {
			int gridPoints = Integer.parseInt(intervalFactoryParams[0]);
			double finalTime = Double.parseDouble(intervalFactoryParams[1]);
			intervalFactory = new FixedGridIntervalFactory(gridPoints, finalTime, printIntervals);
		}
		else {
			System.err.println("# Invalid interval method");
			return null;
		}
		
		return intervalFactory;
	}
	
	
	// returns a single interval [0, \infty)
	public static class SingleIntervalFactory extends IntervalFactory {

		public SingleIntervalFactory (boolean printIntervals) {
			super (printIntervals);
		}
		
		@Override
		public <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			Interval[] toReturn = new Interval[1];
			toReturn[0] = new Interval(0, Double.POSITIVE_INFINITY);
			return toReturn;
		}

		@Override
		public String factoryType() {
			return "single";
		}

		@Override
		public String factoryParams() {
			return null;
		}

	}
	
	public static class LogUniformFactory extends IntervalFactory {

		final private List<Integer> intervalNumbers;
		final private List<Double> intervalBorders;
		
		final private Interval[] intervals;

		public LogUniformFactory (List<Integer> intervalNumbers, List<Double> intervalBorders, boolean printIntervals) throws IOException {
			super (printIntervals);
			
			// remember parameters
			this.intervalNumbers = intervalNumbers;
			this.intervalBorders = intervalBorders;
			assert (intervalNumbers.size()+1 == intervalBorders.size());
			assert (intervalBorders.get(0) > 0);
			
			// build them intervals
			List<Interval> themIntervals = new ArrayList<Interval>();
			
			// add the first interval
			themIntervals.add (new Interval(0d, intervalBorders.get(0)));
			
			// now go through each of the ranges
			for (int rangeIdx=0; rangeIdx<this.intervalNumbers.size(); rangeIdx++) {
				// stuff for this range
				int rangeIntervals = this.intervalNumbers.get(rangeIdx);
				assert (rangeIntervals > 0);
				double rangeLower = this.intervalBorders.get(rangeIdx);
				double rangeUpper = this.intervalBorders.get(rangeIdx+1);
				assert (rangeLower < rangeUpper);
				
				double prevLast = rangeLower;
				for (int i=1; i<=rangeIntervals; i++) {
					double thisLast = 0d;
					if (i == rangeIntervals) {
						thisLast = rangeUpper;
					}
					else {
						 thisLast = MetaOptimization.logify (i/((float)rangeIntervals), rangeLower, rangeUpper);
					}
					themIntervals.add (new Interval(prevLast, thisLast));
					prevLast = thisLast;
				}
			}
			
			// add the last interval
			themIntervals.add (new Interval(intervalBorders.get(intervalBorders.size()-1), Double.POSITIVE_INFINITY));

			// remember
			this.intervals = themIntervals.toArray (new Interval[] {});
		}
		
		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract (Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			return this.intervals;
		}

		@Override
		public String factoryType() {
			return "logUniform";
		}

		@Override
		public String factoryParams() {
			String returnString = "";
			for (Integer currInt : this.intervalNumbers) returnString += "" + currInt + ",";
			for (int i=0; i<(this.intervalBorders.size()-1); i++) returnString += "" + this.intervalBorders.get(i) + ",";
			returnString += "" + this.intervalBorders.get(this.intervalBorders.size()-1);
			return returnString;
		}

	}

	// makes the intervals the same as the demographic epochs (but not the ones of length zero)
	public static class SimpleIntervalFactory extends IntervalFactory {
		
		public SimpleIntervalFactory (boolean printIntervals) {
			super (printIntervals);
		}
		
		@Override
		public <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract (Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme){
			return getNonPulseIntervals(demo);
		}

		
		public static Interval[] getNonPulseIntervals(Demography demo) {
			ArrayList<Interval> toReturn = new ArrayList<Interval>();
			// weed out the empty ones
			for (Interval thisInt : demo.epochList) {
				if (thisInt.endPoint > thisInt.startPoint) toReturn.add(thisInt);
			}
			return toReturn.toArray(new Interval[] {});
		}
		
		@Override
		public String factoryType() {
			return "simple";
		}

		@Override
		public String factoryParams() {
			return null;
		}
	}
	
	static Interval[] getIntervalsFromEndpoints(double[] points) {
		// make intervals from the points
		List<Interval> intervals = new ArrayList<Interval>();
		
		double lastPoint = 0d;
		for (double nextPoint : points) {
			assert nextPoint > lastPoint;
			intervals.add(new Interval(lastPoint, nextPoint));
			lastPoint = nextPoint;
		}
		assert lastPoint < Double.POSITIVE_INFINITY;
		intervals.add(new Interval(lastPoint, Double.POSITIVE_INFINITY));

		// and save them
		return intervals.toArray(new Interval[]{});
	}
	
	public static class CustomFixedIntervalFactory extends IntervalFactory {

		private final Interval[] intervals;
		private final String paramString;
		
		public CustomFixedIntervalFactory (double[] points, boolean printIntervals) {
			super (printIntervals);

			// save a string for the params
			this.paramString = Arrays.toString(points);
			
			this.intervals = getIntervalsFromEndpoints(points);
		}
		
		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			// give them away
			return this.intervals;
		}
		
		@Override
		public String factoryType() {
			return "customFixed";
		}

		@Override
		public String factoryParams() {
			return this.paramString; 
		}
	}

	
	public static class CustomIntervalFactory extends IntervalFactory {

		private final double[] points;
		
		public CustomIntervalFactory (double[] points, boolean printIntervals) {
			super (printIntervals);
			this.points = points;
		}
		
		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			List<Interval> intervals = new ArrayList<Interval>();
			
			double lastPoint = 0d;
			for (double nextPoint : this.points) {
				nextPoint /= Math.max(1d, presentConfig.getPopulation(additionalDeme).totalGeneticTypes());
				assert nextPoint > lastPoint;
				intervals.add(new Interval(lastPoint, nextPoint));
				lastPoint = nextPoint;
			}
			assert lastPoint < Double.POSITIVE_INFINITY;
			intervals.add(new Interval(lastPoint, Double.POSITIVE_INFINITY));
			
			return intervals.toArray(new Interval[]{});
		}

		@Override
		public String factoryType() {
			return "custom";
		}

		@Override
		public String factoryParams() {
			return Arrays.toString(this.points); 
		}
	}
	
	public static class FixedGridIntervalFactory extends IntervalFactory {

		private final int gridPoints;
		private final double finalTime;
		private final Interval[] grid;
		
		public FixedGridIntervalFactory(int gridPoints, double finalTime, boolean printIntervals) {
			super (printIntervals);
			this.gridPoints = gridPoints;
			this.finalTime = finalTime;
			this.grid = new Interval[gridPoints + 1];
			
			double intervalLength = finalTime / (double) gridPoints;
			for (int i = 0; i < gridPoints; i++) {
				this.grid[i] = new Interval(i * intervalLength, (i+1) * intervalLength);
			}
			this.grid[gridPoints] = new Interval(gridPoints * intervalLength, Double.POSITIVE_INFINITY);
		}

		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			return grid;
		}

		@Override
		public String factoryType() {
			return "fixedGrid";
		}

		@Override
		public String factoryParams() {
			return "" + gridPoints + "," + finalTime;
		}

	}

	
	public static class FixedUniformIntervalFactory extends IntervalFactory {

		private final int numIntervals;
		private final double effectivePopSize;
		
		public FixedUniformIntervalFactory(int numIntervals, double effectivePopSize, boolean printIntervals) {
			super (printIntervals);
			this.numIntervals = numIntervals;
			this.effectivePopSize = effectivePopSize;
		}

		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			// get the intervals
			Interval[] thePartition = RealPartition.getExponentialPartition(numIntervals, Math.max(1d, presentConfig.getPopulation(additionalDeme).totalGeneticTypes()) / this.effectivePopSize);

			// give it away now
			return thePartition;
		}

		@Override
		public String factoryType() {
			return "fixedUniform";
		}

		@Override
		public String factoryParams() {
			return "" + numIntervals + "," + effectivePopSize;
		}

	}

	public <H extends GeneticType> Interval[] getIntervals (Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
		if (presentConfig.getTotalGeneticTypes() == 0) {
			return null;
		}
		else {
			return getNontrivialIntervals (demo, presentConfig, additionalDeme);
		}
	}
	
	
	// This computes a mean rate for each epoch in a demography, then constructs a balanced partition based on those rates
	public static class MeanRateIntervalFactory extends IntervalFactory {		
		
		public MeanRateIntervalFactory (int numIntervals, boolean printIntervals) {
			super (printIntervals);
			this.numIntervals = numIntervals;
		}
		
		@Override
		public <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract (Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			
			assert (presentConfig.getTotalGeneticTypes() > 0);
			
			ArrayList<Interval> intervalList = new ArrayList<Interval>();
			
			// get balanced intervals
			Interval[] rateOneIntervals = RealPartition.getBalancedPartition(numIntervals);

			double currRateOneTime = rateOneIntervals[0].startPoint;
			int currRateOneIdx = 0;
			assert (currRateOneTime == 0d);
			
			// start off
			double currTime = 0d;
			double lastTime = 0d;
			int currEpochIdx = 0;
			
			// initial rate
			double currRate = getMeanRate(currEpochIdx, demo, presentConfig);
			
			// loop and step through rate-one intervals and real intervals simultaneously
			while (true) {
				// whats the stepsize in the rate-one world
				double currRateOneIncrement = rateOneIntervals[currRateOneIdx].endPoint - currRateOneTime;
				double currEpochEndPoint = demo.epochList[currEpochIdx].endPoint;
				
				// eat of a piece with the given rate (we might overstep, and thats ok)
				double nextTime = currTime + currRateOneIncrement / currRate;
				
				// step in real world
				// are we still in the currEpoch, after eating this bit
				if (nextTime <= currEpochEndPoint){
					intervalList.add(new Interval(lastTime, nextTime));
					lastTime = nextTime;
					
					currRateOneIdx++;
				}
				// did we go over the epoch end and have to change the epoch
				if (nextTime >= currEpochEndPoint){					
					nextTime = currEpochEndPoint;					
					currEpochIdx++;
				}

				// step in the rate one world
				currRateOneTime += (nextTime - currTime) * currRate;
				currTime = nextTime;

				
				// make sure that the ending is proper
				if (currEpochIdx >= demo.epochList.length){
					assert (currRateOneIdx == rateOneIntervals.length);
					assert (currTime == Double.POSITIVE_INFINITY);
					assert (currRateOneTime == Double.POSITIVE_INFINITY);
					break;
				}
				
				
				// update rate
				// but only if it is a non pulse deme (that we stepped over in the real world in the last round)
				if (!demo.isPulse(currEpochIdx)) {
					currRate = getMeanRate (currEpochIdx, demo, presentConfig);
				}
			}
			
			// return it
			return intervalList.toArray(new Interval[] {});						
		}

		
		public <H extends GeneticType> double getMeanRate(int epochIdx, Demography demo, DemoConfiguration<H> presentConfig) {
			
			if (demo.isPulse(epochIdx)) {
				return 1d;
			}
			
			double rateSum = 0d;
			for (int ancientDeme = 0; ancientDeme < demo.treePartitionList.get(epochIdx).size(); ancientDeme++){
				for (int presentDeme : demo.treePartitionList.get(epochIdx).get(ancientDeme) ){
					rateSum += presentConfig.getPopulation(presentDeme).totalGeneticTypes() / demo.popSizes.get(epochIdx)[ancientDeme];
				}
			}
			return rateSum / demo.treePartitionList.get(epochIdx).size();
		}

		
		
		@Override
		public String factoryType() {
			return "meanrate";
		}

		@Override
		public String factoryParams() {
			return "" + numIntervals;
		}

		private final int numIntervals;
	}
	
	// for each epoch we empirically invert the cumulative distribution function of the competing absorbing exponentials from the different demes
	// and use this to get the appropriate number of intervals
	public static class MixtureUniformIntervalFactory extends IntervalFactory {

		public MixtureUniformIntervalFactory (int numIntervals, double EPSILON, boolean printIntervals){
			super (printIntervals);
			assert (numIntervals > 0);
			this.numIntervals = numIntervals;
			this.EPSILON = EPSILON;
		}
		
		
		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			assert (presentConfig.getTotalGeneticTypes() > 0);
			
			ArrayList<Interval> intervalList = new ArrayList<Interval>();

			// initialize stuff
			double currTime = 0d;
			double lastTime = 0d;
			int currEpochIdx = 0;
			
			double remainingMass = 1d / this.numIntervals;
			
			while (true) {				
				assert (intervalList.size() <= this.numIntervals - 1);
				
				// should we add last interval, and break the loop?
				if (intervalList.size() >= this.numIntervals - 1) {
					// remaining mass should be close to end
					assert (Math.abs(remainingMass - 1d) < this.EPSILON);
					// then add the interval
					intervalList.add(new Interval(lastTime, Double.POSITIVE_INFINITY));
					break;
				}

				if (demo.isPulse(currEpochIdx)) {
					currEpochIdx++;
					continue;
				}
				
				// get the absorption rates in those demes, where we can get absorbed
				double[] rates = getNonzeroRates(currEpochIdx, demo, presentConfig, this.EPSILON);
				
				// we need the exponential sum w such that 1 - (w / rates.length) = remainingMass
				double w = rates.length * (1d - remainingMass);
				// and add as much time as the inverse of the exponential sum suggests
				double nextTime = currTime + inverseExponentialSum(rates, w, this.EPSILON);
				
				// see where the endpoint of the current epoch is
				double currEpochEndPoint = demo.epochList[currEpochIdx].endPoint;
				
				// either we ate up remaining mass before hitting the end of epoch or the end of epoch is nearer
				if (nextTime <= currEpochEndPoint + this.EPSILON){
					// we ate up the mass before the endpoint
					
					// so add an interval
					intervalList.add(new Interval(lastTime, nextTime));
					lastTime = nextTime;
					
					// reset remaining mass
					remainingMass = 1d / (this.numIntervals - intervalList.size());
					
					// should we go to next epoch?
					if (nextTime >= currEpochEndPoint - this.EPSILON) {
						currEpochIdx++;
					}
					
				}
				// or the end of epoch comes first
				else {
					// just increase the epoch and adjust the remaining mass
					nextTime = currEpochEndPoint;					
					currEpochIdx++;
					// eat of as much as you get until the end of the epoch
					remainingMass -= (1 - exponentialSum(rates, (nextTime - currTime), false) / rates.length);
										
				}
				
				currTime = nextTime;				
			}
						
			return intervalList.toArray(new Interval[] {});									
		}
		
		// gets all the nonzero rates
		private static final <H extends GeneticType> double[] getNonzeroRates (int epochIdx, Demography demo, DemoConfiguration<H> presentConfig, double EPSILON){
			ArrayList<Double> rateList = new ArrayList<Double>();
			
			for (int ancientDeme = 0; ancientDeme < demo.treePartitionList.get(epochIdx).size(); ancientDeme++){
				double currRate = 0d;
				for (int presentDeme : demo.treePartitionList.get(epochIdx).get(ancientDeme) ){
					currRate += presentConfig.getPopulation(presentDeme).totalGeneticTypes();
				}
				
				currRate /= demo.popSizes.get(epochIdx)[ancientDeme];
				
				if (currRate > EPSILON) {
					rateList.add(currRate);
				}
			}
			
			double[] toReturn = new double[rateList.size()];
			for (int i = 0; i < toReturn.length; i++){
				toReturn[i] = rateList.get(i);
			}
			
			return toReturn;
		}
		
		
		// this finds the x such that w = \sum_i^n \exp( - rates[i] * x), assuming all rates are < 0, and n >= w >= 0
		private static final double inverseExponentialSum(double[] rates, double w, double EPSILON){
			
			assert (rates.length > 0);
			for (double rate : rates){
				assert (rate > 0d);
			}
			
			assert (w <= rates.length && w >= 0);
			if (w == rates.length) return 0d;
			if (w == 0d) return Double.POSITIVE_INFINITY;
						
			double currPoint = 0d;
			double currValue = exponentialSum(rates, currPoint, false);
			
			while (true) {
				
				if (Math.abs(currValue - w) < EPSILON) return currPoint;
				
				// update currPoint
				double derivative = exponentialSum(rates, currPoint, true);
				
				double nextPoint = currPoint + (w - currValue) / derivative;
				
				currPoint = nextPoint;
				currValue = exponentialSum(rates, currPoint, false);
			}
			
		}
		
		// returns \sum_i \exp( - rates[i] * x), or the derivative of this function
		private static double exponentialSum(double[] rates, double x, boolean isDerivative) {
			double sum = 0d;
			for (double rate : rates){
				assert (rate > 0d);
				sum += Math.exp(- rate * x) * (isDerivative ? - rate : 1d);
			}
			return sum;
		}
		
		
		@Override
		public String factoryType() {
			return "mixtureuniform";
		}


		@Override
		public String factoryParams() {
			return "" + numIntervals;
		}


		private final int numIntervals;
		private final double EPSILON;
		
	}
	
	public static class CombinedIntervalFactory extends IntervalFactory {

		private final List<IntervalFactory> memberFactories;
		private final double EPSILON;
		
		public CombinedIntervalFactory(List<IntervalFactory> memberFactories, double EPSILON, boolean printIntervals) {
			super (printIntervals);
			this.memberFactories = memberFactories;
			this.EPSILON = EPSILON;			
		}
		
		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(
				Demography demo, DemoConfiguration<H> presentConfig,
				int additionalDeme) {
			
			// create a dummy demography with only one epoch
			Demography tmpDemo = new ConstDemography(1, 0);
			
			// go thru the interval factories
			for (IntervalFactory factory : memberFactories) {
				// get the current intervals
				Interval[] currIntervals = factory.getIntervals (demo, presentConfig, additionalDeme);
				
				// refine the dummy demography
				tmpDemo = new Demography(tmpDemo, currIntervals, null, EPSILON);
			}
			
			return tmpDemo.epochList;
		}

		@Override
		public String factoryType() {
			return "combined";
		}

		@Override
		public String factoryParams() {
			String toReturn = "";
			for (IntervalFactory factory : memberFactories) {
				toReturn += factory.factoryType() + ":" + factory.factoryParams() + ";";
			}
			return toReturn;
		}
		
	}
		
	
	// old heuristic used for the structured guy
	// just mean absorption rate (w/o any cake) of present demes
	public static class OldIntervalFactory extends IntervalFactory {
		
		public OldIntervalFactory(int numAdditionalIntervals, boolean printIntervals){
			super (printIntervals);
			this.numAdditionalIntervals = numAdditionalIntervals;
		}	
		
		@Override
		public <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			Interval[] refineByIntervals = getRefineByIntervals(demo, presentConfig);
			
			// forget about the mapping later
			Demography tmpDemo = new Demography(demo, refineByIntervals, new ArrayList<Integer>(), UberDemographyCore.ONELOCUS_EPSILON);
			
			// filter the intervals, drop out empty
			ArrayList<Interval> toReturn = new ArrayList<RealPartition.Interval>();
			for (int e=0; e<tmpDemo.epochList.length; e++) {
				if (!tmpDemo.isPulse(e)) toReturn.add(tmpDemo.epochList[e]);
			}
			
			// and return it
			return toReturn.toArray (new Interval[] {});			
		}
		
		private <H extends GeneticType> Interval[] getRefineByIntervals(Demography demo, DemoConfiguration<H> presentConfig) {
			
			// get the rate for the refinement
			double expRate = 0d;
			// for now take harmonic mean
			// WARNING be aware of the ghost demes, they can mess with this, we should include something, though
			for (int i=0; i<presentConfig.getNumberDemes(); i++) {
				if (presentConfig.getPopulation(i).totalGeneticTypes() == 0) {
					// we have a ghost deme, so lets put some migration rates out of this deme
					// minus cause diagonal of migration rate is negative
					// the two is just so it works better
					if (demo.migrationMatrixList.get(0)[i][i] == 0) {
						// for now do this 
						continue;
					}
					else {
						expRate -= 1/(2*demo.migrationMatrixList.get(0)[i][i]);
					}
				}
				else {
					// normal deme
					expRate += demo.popSizes.get(0)[i] / presentConfig.getPopulation(i).totalGeneticTypes();
				}
			}
			// and normalize it
			expRate = presentConfig.getNumberDemes() / expRate;
			// should have a nice expRate
			
			// get the refinement
			// additional intervals is the number of breakpoints for the refinement
			Interval[] refinement = RealPartition.getExponentialPartition (numAdditionalIntervals+1, expRate);
			
			return refinement;
		}
		
		
		
		@Override
		public String factoryType() {
			return "old";
		}

		@Override
		public String factoryParams() {
			return "" + numAdditionalIntervals;
		}



		private final int numAdditionalIntervals;
	}
	
	
	
	
	
	// this was the old way we got cake intervals. keep it around so we can compare with previous output
	public static class OldCakeIntervalFactory extends IntervalFactory {
		
		final int additionalIntervals;
		
		public OldCakeIntervalFactory(int additionalIntervals, boolean printIntervals) {
			super (printIntervals);
			this.additionalIntervals = additionalIntervals;
		}

		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			return RealPartition.getBalancedPartition(additionalIntervals+1);
		}

		@Override
		public String factoryType() {
			return null;
		}

		@Override
		public String factoryParams() {
			return null;
		}
	}
	
	
	public static class QuadratureFactory extends IntervalFactory {

		final int numIntervals;
		final double effPopSize;
		
		public QuadratureFactory(int numIntervals, double effPopSize, boolean printIntervals) {
			super (printIntervals);
			this.numIntervals = numIntervals;
			this.effPopSize = effPopSize;
		}

		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(Demography demo, DemoConfiguration<H> presentConfig, int additionalDeme) {
			
			int nTrunk = presentConfig.getTotalGeneticTypes();
			
			Interval[] qPoints = new Interval[this.numIntervals];
			
			double currStartPoint = 0d;
			double currEndPoint = 0d;
			
			int tIdx = 0;
			for (QuadraturePoint q : LaGuerre.getQuadraturePoints(this.numIntervals)){
				
				double nbar = getNbar(nTrunk, currStartPoint, demo);
				double increment = (q.endPoint - q.startPoint) / (nbar / this.effPopSize);

				currEndPoint = currStartPoint + increment;
				qPoints[tIdx++] = new Interval(currStartPoint, currEndPoint);
						
				currStartPoint = currEndPoint;
			}
		
			return qPoints;
		}

		@Override
		public String factoryType() {
			return "quadrature";
		}

		@Override
		public String factoryParams() {
			return "" + numIntervals + "," + effPopSize;
		}

		
		protected double getNbar (int nTrunk, double time, Demography demo) {
			return nTrunk * TrunkProcess.getCakeFactors(nTrunk, new double[] {time})[0];
		}
	}
	
	public static class QuadratureFactory2 extends QuadratureFactory {

		public QuadratureFactory2(int numIntervals, double effPopSize, boolean printIntervals) {
			super(numIntervals, effPopSize, printIntervals);
		}

		@Override
		public String factoryType() {
			return "quadrature2";
		}

		@Override
		protected double getNbar(int nTrunk, double time, Demography demo) {
			double[] popSizes = new double[demo.popSizes.size()];
			
			for (int i = 0; i < popSizes.length; i++) {
				// assert one population only
				assert demo.popSizes.get(i).length == 1;
				popSizes[i] = demo.popSizes.get(i)[0];
			}
			
			return nTrunk * TrunkProcess.getCakeFactors(nTrunk, new double[] {time}, demo.epochList, popSizes)[0];
		}
	
		
		
	}
	
	// makes epochs small enough so no deme grows by more than a factor of k within an epoch
	public static class ExponentialGrowthIntervalFactory extends IntervalFactory{

		final double growthFactor;
		
		public ExponentialGrowthIntervalFactory(double growthFactor, boolean printIntervals) {
			super(printIntervals);
			this.growthFactor = growthFactor;
		}

		@Override
		public String factoryType() {
			return "exponentialGrowth";
		}

		@Override
		public String factoryParams() {
			return "" + growthFactor; 
		}
		
		private static double[] getExpGrowthPoints(Demography demo, double growthFactor) {
			
			ArrayList<Double> pointsList = new ArrayList<Double>();
			
			for (int epochIdx = 0; epochIdx < demo.numIntervals(); epochIdx++) {
				if (demo.isPulse(epochIdx)) continue;
				
				double endPoint = demo.epochList[epochIdx].endPoint;
				if (endPoint == Double.POSITIVE_INFINITY) continue;
				
				double maxGrowthRate = 0d;
				for (double growthRate : demo.expGrowthRates.get(epochIdx)) {
					growthRate = Math.abs(growthRate);
					if (growthRate > maxGrowthRate) maxGrowthRate = growthRate;
				}
				
				if (maxGrowthRate == 0d) continue;
				
				// how long does it take population to change by growthFactor?
				// popSize = ancestralPop*Math.exp(growthRate*(this.epochList[epochIdx].endPoint - currTime))
				// so need Math.exp(growthRate * t) = growthFactor
				// so growthRate * t = log(growthFactor)
				double halfLife = Math.abs(Math.log(growthFactor)) / maxGrowthRate;
				
				double startPoint = demo.epochList[epochIdx].startPoint;
				for (double currTime = startPoint + halfLife; currTime < endPoint; currTime += halfLife) {
					pointsList.add(currTime);
				}
				pointsList.add(endPoint);
			}
			
			double[] toReturn = new double[pointsList.size()];
			for (int i = 0; i < pointsList.size(); i++) {
				toReturn[i] = pointsList.get(i);
			}
			return toReturn;
		}

		@Override
		protected <H extends GeneticType> Interval[] getNontrivialIntervalsAbstract(
				Demography demo, DemoConfiguration<H> presentConfig,
				int additionalDeme) {
			return getIntervalsFromEndpoints(getExpGrowthPoints(demo, growthFactor));
		}
		
	}
}