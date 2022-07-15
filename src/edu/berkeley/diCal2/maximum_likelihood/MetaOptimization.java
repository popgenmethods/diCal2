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

package edu.berkeley.diCal2.maximum_likelihood;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.Switch;
import com.martiansoftware.jsap.stringparsers.FileStringParser;

import edu.berkeley.diCal2.csd.UberDemographyCore;
import edu.berkeley.diCal2.csd.auxiliary.DelayedRandom;
import edu.berkeley.diCal2.demography.DemographyFactory;
import edu.berkeley.diCal2.maximum_likelihood.StructureEstimationEM.DiCalObjectiveFunction;
import edu.berkeley.diCal2.maximum_likelihood.StructureEstimationEM.EmResult;

public class MetaOptimization {

	public static final int MAX_NEW_POINT_TRIES = 50;
	public static final double GAUSSIAN_FACTOR = 1d;
	public static final double GAUSSIAN_DISPLACE_FACTOR = 1;
	public static final double GAUSSIAN_STRETCH_FACTOR = 5 * GAUSSIAN_DISPLACE_FACTOR;
	static double EPSILON = UberDemographyCore.DEFAULT_RENORMALIZATION_EPSILON;
	
	public static class MetaOptimizationArguments {
		public boolean gridStart;
		public Integer numStartPoints = null;
		public double[][] startingRanges;
		public int numMetaIterations;
		public int metaParallelEmSteps;
		public Integer metaKeepBest = null;
		public Integer metaNumPoints = null;
		public double stretchProportion;
		public double disperseFactor;
		public File metaStartFile = null;
		public double sdPercentage;
		
		public Parameter[] register (Parameter[] additionalParams) {

			// we need some file parsing things
			FileStringParser myFileParser = FileStringParser.getParser();
			myFileParser.setMustBeFile(true);
			myFileParser.setMustExist(true);

			// here go the new parameters
			Parameter[] metaParams = new Parameter[] {
					new FlaggedOption ("metaStartFile", myFileParser, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaStartFile", "A file containing the start points to be used in tab-separated format." ),
					new FlaggedOption ("metaNumStartPoints", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaNumStartPoints", "The number of starting points (in each dimension if we have a grid)."),
					new FlaggedOption ("metaNumIterations", JSAP.INTEGER_PARSER, "1", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaNumIterations", "The number of meta optimization steps for the gentic algorithm."),
					new FlaggedOption ("metaParallelEmSteps", JSAP.INTEGER_PARSER, "1", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaParallelEmSteps", "Number of EM steps to be exectued in parallel during gentic algorithm."),
					new FlaggedOption ("metaKeepBest", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaKeepBest", "Number of best points to keep after a meta step."),
					new FlaggedOption ("metaNumPoints", JSAP.INTEGER_PARSER, JSAP.NO_DEFAULT, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaNumPoints", "Number of points to be used in each meta iteration."),
					new Switch ("metaGridStart", JSAP.NO_SHORTFLAG, "metaGridStart", "If flag is set we start with a grid of points, otherwise the start points are chosen randomly.")
			};
			
			// add them to the list
			List<Parameter> paramList = new ArrayList<Parameter>();
			for (Parameter param : metaParams) paramList.add(param);
			for (Parameter param : additionalParams) paramList.add(param);
			
			// give them away
			return paramList.toArray (new Parameter[] {});
		}
		
		public Parameter[] registerHidden(Parameter[] additionalHiddenParams) {
			Parameter[] metaHiddenParams = new Parameter[] {
					new FlaggedOption ("metaStretchProportion", JSAP.DOUBLE_PARSER, "0.25", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaStretchProportion", "Probability that a new point is chosen in a ridge-like manner, as opposed to each dimension independently."),
					new FlaggedOption ("metaDisperseFactor", JSAP.DOUBLE_PARSER, "1", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaDisperseFactor", "SD for new Gaussian points is multiplied by this factor."),
					new FlaggedOption ("metaSDPercentageIfZero", JSAP.DOUBLE_PARSER, "0.2", JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG, "metaSDPercentageIfZero", "If SD is zero (maybe only one point), then it is replaced by this percentage of the total value of the best point.")
			};
			
			// add them to the list
			List<Parameter> paramHiddenList = new ArrayList<Parameter>();
			for (Parameter param : metaHiddenParams) paramHiddenList.add(param);
			for (Parameter param : additionalHiddenParams) paramHiddenList.add(param);
			
			// give them away
			return paramHiddenList.toArray (new Parameter[] {});
		}

		
		void get (JSAPResult jsapParams, double[][] bounds) throws IOException {
			// copy them bounds, for now
			if (bounds != null) {
				this.startingRanges = new double[bounds.length][bounds[0].length];
				for (int i=0; i<bounds.length; i++) {
					this.startingRanges[i] = bounds[i];
				}
			}
			
			// get your parameters
			if (jsapParams.contains ("metaStartFile")) this.metaStartFile = jsapParams.getFile("metaStartFile");
			this.numMetaIterations = jsapParams.getInt("metaNumIterations");
			if (jsapParams.contains("metaNumStartPoints")) this.numStartPoints = jsapParams.getInt("metaNumStartPoints");
			this.gridStart = jsapParams.getBoolean("metaGridStart");
			this.metaParallelEmSteps = jsapParams.getInt("metaParallelEmSteps");
			if (jsapParams.contains("metaKeepBest")) this.metaKeepBest = jsapParams.getInt("metaKeepBest");
			if (jsapParams.contains("metaNumPoints")) this.metaNumPoints = jsapParams.getInt("metaNumPoints");
			this.stretchProportion = jsapParams.getDouble("metaStretchProportion");
			this.disperseFactor = jsapParams.getDouble("metaDisperseFactor");
			this.sdPercentage = jsapParams.getDouble("metaSDPercentageIfZero");
			
			boolean moreThanOneMeta = this.numMetaIterations > 1;
			boolean moreThanOneStart = ((this.numStartPoints != null) && (this.numStartPoints > 1)) || (this.metaStartFile != null);
			if (moreThanOneMeta) {
				if (!moreThanOneStart) {
					throw new IOException ("Meta iteration only allowed with more than one start point.");
				}
				if (this.metaKeepBest == null) {
					throw new IOException ("Need to specify metaKeepBest for meta iteration.");
				}
				else if (this.metaNumPoints == null) {
					throw new IOException ("Need to specify metaNumPoints for meta iteration.");
				}
				else {
					if (this.metaKeepBest > this.metaNumPoints) {
						throw new IOException ("metaNumPoints has to be bigger then metaKeepBest.");
					}
				}
			}
			else {
				if (this.metaKeepBest != null) {
					throw new IOException ("No metaKeepBest allowed without meta iteration.");
				}
				if (this.metaNumPoints != null) {
					throw new IOException ("No metaNumPoints allowed without meta iteration.");
				}
			}
			if (this.metaStartFile != null) {
				if (this.numStartPoints != null || this.gridStart) {
					throw new IOException ("If a file for the start points is specified, no numStartPoints or gridStart allowed.");
				}
			}
		}

		public void print (PrintStream outStream) {
			outStream.println("# numMetaIterations = " + this.numMetaIterations);
			outStream.println("# numStartPoints = " + this.numStartPoints);
			outStream.println("# metaParallelEmSteps = " + this.metaParallelEmSteps);
			outStream.println("# gridStart = " + this.gridStart);
		}

	}

	public static double[] addRecoMutEstimation(double[] startingPoint, boolean estimateRecomScaling, int thetaDim, double[] thetaInitial) {
		double[] currPoint= new double[startingPoint.length + (estimateRecomScaling ? 1 : 0) + thetaDim];
		// the demography parameter part
		for (int i = 0; i < startingPoint.length; i++) {
			currPoint[i] = startingPoint[i];
		}
		// recom?
		if (estimateRecomScaling) {
			currPoint[currPoint.length - thetaDim - 1] = 1d;
		}
		//Initialize theta to be what was read in from param file
		for(int i = 0; i < thetaDim; i++){
			currPoint[currPoint.length - thetaDim + i] = thetaInitial[i];
		}
		return currPoint;
	}

	public static List<double[]> buildGrid (double[][] startingRanges, int numStartPoints) {

		assert (numStartPoints > 1);
		
		// TODO we have to check for time ordered parameters, since they might not be allowed
		// then again, it's of ok to just return -Infty; just have to be aware of it
		
		int dimension = startingRanges.length;
		
		// to go here
		List<double[]> points = new ArrayList<double[]>();
		
		// gets initialized with zero anyways
		int[] currPoint = new int[dimension];
		// fancy pants odometer algorithm
		boolean finished = false;
		while (!finished) {
			// add a points for this
			double[] unitIntervalPoint = new double[currPoint.length];
			for (int i=0; i<unitIntervalPoint.length; i++) unitIntervalPoint[i] = currPoint[i]/(float)(numStartPoints-1);
			
			points.add (unitIntervalToGrid (unitIntervalPoint, startingRanges, numStartPoints));
			
			// and go to next
			for (int i=0; i<dimension; i++) {
				currPoint[i]++;
				if (currPoint[i] >= numStartPoints) {
					currPoint[i] = 0;
					if (i == dimension-1) finished = true;
				}
				else {
					break;
				}
			}
		}
		
		// give it away now
		return points;
	}

	private static double[] unitIntervalToGrid (double[] unitPoint, double[][] startingRanges, int numStartPoints) {
		// transform it
		double[] transformedPoint = new double[unitPoint.length];
		assert (unitPoint.length == startingRanges.length);
		for (int i=0; i<unitPoint.length; i++) {
			// logify
			assert ((startingRanges[i].length == 2) && (startingRanges[i][0] < startingRanges[i][1]));
			transformedPoint[i] = logify (unitPoint[i], startingRanges[i][0], startingRanges[i][1]);
		}
		return transformedPoint;
	}

	public static double logify (double fraction, double min, double max) {
		// Round[max^Range[Log[max, min], 1, (1 - Log[max, min])/(numPoints - 1)]*10^digits]/10^digits
		// upperBound can't be one
		if (Math.abs(max - 1d) < EPSILON) {
			if (max - 1d < 0) max = 1d + EPSILON;
			else max = 1d - EPSILON;
		}
		// no calculate
		double oneTo = 1d;
		double zeroTo = myLog  (max, min);
		double slope = oneTo - zeroTo;
		double gridValue = Math.pow (max, slope * fraction + zeroTo);
		
		// sometimes it is inaccurate
		if (gridValue < min) gridValue = min;
		if (gridValue > max) gridValue = max;
		
		return gridValue;
	}

	private static double myLog (double base, double x) {
		return Math.log(x)/Math.log(base);
	}

	public static List<double[]> randomPoints (double[][] startingRanges, int numStartPoints, DemographyFactory demoFactory, boolean estimateRecomScaling, boolean estimateTheta, int thetaDim, DelayedRandom rGen) {

		int dim = startingRanges.length;
		
		List<double[]> pointList = new ArrayList<double[]>();
		
		// make a number of points
		for (int i=0; i<numStartPoints; i++) {
			
			// make a new point
			double[] thePoint = null;
			// try a few times
			int tries = 0;
			while ((thePoint == null) || (demoFactory.getDemography (DiCalObjectiveFunction.getDemoPointStatic (thePoint, estimateRecomScaling, estimateTheta, thetaDim)) == null) || (tries >= MetaOptimization.MAX_NEW_POINT_TRIES)) {
				// make this try
				thePoint = new double[dim];
				for (int j=0; j<thePoint.length; j++) {
					double rand = rGen.nextDouble();
					thePoint[j] = logify (rand, startingRanges[j][0], startingRanges[j][1]);
				}
				// count up
				tries++;
			}
			// if we hit the max, issue warning
			if (tries >= MetaOptimization.MAX_NEW_POINT_TRIES)  {
				System.out.println ("# [WARNING] Could not find valid start point within max iterations.");
			}
			// add the point
			pointList.add(thePoint);
		}
		
		return pointList;
	}

	public static synchronized void synchronizedPrintln (String debugString) {
		System.out.println (debugString);
	}

	public static List<double[]> nextGeneration (List<EmResult> prevGeneration, int numKeepBest, int metaNumPoints, double stretchProportion, double disperseFactor, double sdPercentage, DemographyFactory demoFactory, boolean estimateRecomScaling, boolean estimateTheta, int thetaDim, DelayedRandom rGen) {

		assert (numKeepBest >= 1);
		
		// sort the results by likelihood, smallest first
		Collections.sort (prevGeneration, new CompareResults());
		// want highest first
		Collections.reverse (prevGeneration);
						
		
		System.out.println ("# ==== Sorted generation:");
		for (EmResult test : prevGeneration) {
			System.out.println("# " + test.logLikelihood + "\t" + Arrays.toString(test.mle));
		}
		
		
		// get the best out
		System.out.println ("# ==== we take:");
		List<double[]> theBest = new ArrayList<double[]>();
		for (int i=0; i<Math.min(numKeepBest,prevGeneration.size()); i++) {
			theBest.add (prevGeneration.get(i).mle);
			System.out.println ("# " + prevGeneration.get(i).logLikelihood + "\t" + Arrays.toString(prevGeneration.get(i).mle));
		}
		
		// we need the standard deviation of the best points
		double[] sd = marginalStandardDeviation  (theBest);
		
		System.out.println("# ==== sd: " + Arrays.toString(sd));
		
		// see whether sds are zero and we need to get a new effective sd
		double[]	 theBestPoint = 	prevGeneration.get(0).mle;
		// go through all individual sds
		assert (sd.length == theBestPoint.length);
		for (int j=0; j<sd.length; j++) {
			// is it zero?
			if (sd[j] < EPSILON) {
				// yes, so put something in there
				// percentage of coordinate of best point 
				sd[j] = sdPercentage * theBestPoint[j];
			}
		}
				
		System.out.println("# ==== effective sd: " + Arrays.toString(sd));
		
		// here goes the next generation
		List<double[]> nextGeneration = new ArrayList<double[]>();
		
		// take the best from the previous run (we want to remember the best)
		for (double[] p : theBest) nextGeneration.add(p);
		
		int dim = theBest.get(0).length;
		
		// deal with them bounds
		double[][] bounds = new double[dim][2];
		double[][] demoBounds = demoFactory.getBounds();
		// copy the demoBounds in the beginning
		for (int i=0; i<demoBounds.length; i++) {
			bounds[i] = demoBounds[i];
		}
		// for now, all the reco and mut stuff, if any, gets no bounds
		for (int i=demoBounds.length; i<dim; i++) {
			// this should work
			bounds[i] = null;
		}
		
		// and populate it until full
		while (nextGeneration.size() < metaNumPoints) {
			
			// first choose a random parent
			int parentIdx = (int)(rGen.nextDouble() * theBest.size());
			double[] parent = theBest.get(parentIdx);
			
			// find a new point
			double[] newPoint = null;
			int tries = 0;
			while ((newPoint == null) || (demoFactory.getDemography (DiCalObjectiveFunction.getDemoPointStatic (newPoint, estimateRecomScaling, estimateTheta, thetaDim)) == null) || (tries >= MetaOptimization.MAX_NEW_POINT_TRIES)) {
				newPoint = getNewPoint (parent, sd, bounds, stretchProportion, disperseFactor, rGen);
				tries++;
			}

			// put him in the generation
			nextGeneration.add (newPoint);
		}
		
		
		System.out.println ("# ==== next generation: ");
		for (double[] point : nextGeneration) {
			System.out.println("# " + Arrays.toString(point));
		}
		
		
		
		// return it
		return nextGeneration;
	}
	
	public static double[] getNewPoint (double[] parent, double[] sd, double[][] bounds, double stretchProportion, double disperseFactor, DelayedRandom rGen) {
		assert (parent.length == sd.length);
		assert (parent.length == bounds.length);
		assert ((0d <= stretchProportion) && (stretchProportion <= 1d));
		int dim = parent.length;
		double[] newPoint = new double[dim];
		// then choose whether we do a ridge step, or a random step
		if (rGen.nextDouble() > stretchProportion) {
			// no stretch, one gaussian for everybody
			// add something everywhere
			for (int d=0; d<newPoint.length; d++) {
				double rand = rGen.nextGaussian();
				double theGaussian = Math.sqrt(parent[d] * sd[d]) * disperseFactor * rand * MetaOptimization.GAUSSIAN_DISPLACE_FACTOR;
				newPoint[d] = reflectAtBounds (parent[d] + theGaussian, bounds[d]);
			}
		}
		else {
			// yes stretch, get one gaussian
			double theGaussian = rGen.nextGaussian();
			double avgSd = 1;
			for (double v : sd) avgSd *= v;
			avgSd = Math.pow (avgSd, 1d/sd.length);
			double stretchFactor = reflectAtBounds (avgSd * disperseFactor * theGaussian * MetaOptimization.GAUSSIAN_STRETCH_FACTOR + 1, new double[] {0d, Double.POSITIVE_INFINITY});
			// add something everywhere
			for (int d=0; d<newPoint.length; d++) {
				newPoint[d] = reflectAtBounds (parent[d] * stretchFactor, bounds[d]);
			}
		}
		return newPoint;
	}
	
	public static double reflectAtBounds (double x, double[] bounds) {
		double returnValue = x;
		// no problem here
		if (bounds == null) {
			 return returnValue;
		}
		// normal cases
		if (returnValue < bounds[0]) {
			returnValue = bounds[0] + (bounds[0] - returnValue);
		}
		else if (bounds[1] < returnValue) {
			returnValue = bounds[1] - (returnValue - bounds[1]);
		}
		// now it should be ok or reflected, unless we chose the bounds stupidly
		if ((returnValue < bounds[0]) || (bounds[1] < returnValue)) {
			// whatever
			returnValue = (bounds[0] + bounds[1])/2d;
		}
		// now we sould definitely have something in bounds
		assert ((bounds[0] <= returnValue) && (returnValue <= bounds[1]));
		return returnValue;
	}
	
	private static double[] marginalStandardDeviation (List<double[]> theBest) {
		assert (theBest.size() > 0);
		int n = theBest.size();
		if (n == 1) {
			// should be initialized to zero
			return new double[theBest.get(0).length];
		}
		else {
			int dim = theBest.get(0).length;
			// this is the real case
			double[] mean = new double[dim];
			double[] sd = new double[dim];
			
			// now go through
			for (int d=0; d<dim; d++) {
				// first get the mean
				mean[d] = 0;
				for (int i=0; i<n; i++) {
					mean[d] += theBest.get(i)[d];
				}
				mean[d] /= n;
				
				// then the variance
				sd[d] = 0;
				for (int i=0; i<n; i++) {
					sd[d] += Math.pow (theBest.get(i)[d] - mean[d], 2);
				}
				// lets do n instead of n-1
				sd[d] /= n;
				// don't forget root
				sd[d] = Math.sqrt(sd[d]);
			}
			// done
			return sd;
		}
	}

	public static class CompareResults implements Comparator<EmResult> {
		@Override
		public int compare (EmResult arg0, EmResult arg1) {
			assert ((arg0 != null) && (arg1 != null));
			return Double.compare (arg0.logLikelihood, arg1.logLikelihood);
		}		
	}

	
	public static List<List<String>> TSVReader (File fileName) throws IOException {
		
		Set<Character> commentCharacterSet = new TreeSet<Character>();
		for (char c : new char[] {'#'}) commentCharacterSet.add(c);

		// get a reader for that file
		BufferedReader bufferedCsdReader = new BufferedReader(new FileReader(fileName));
				
		List<List<String>> toReturn = new ArrayList<List<String>>();
		// go through the file
		String line = null;
		while ((line = bufferedCsdReader.readLine()) != null) {
			// always ignore comment and empty lines
			if (line.trim().equals("") || commentCharacterSet.contains(line.charAt(0))) continue;
			
			// now we should have a list of tab seperated strings
			String[] fields = line.split("\\s+");
			List<String> thisList = new ArrayList<String>();
			toReturn.add(thisList);
			for (String value : fields) thisList.add (value);
		}
		// now we should have our list		
		
		// be nice
		bufferedCsdReader.close();
		
		return toReturn;
	}
	
	public static List<double[]> readPointsFromFile (File fileName, int dim) throws IOException {
		
		List<double[]> pointList = new ArrayList<double[]>();
		
		// read it from file
		List<List<String>> tsvStrings = TSVReader (fileName);
		
		// and convert the strings to doubles
		for (List<String> line : tsvStrings) {
			double[] point = new double[dim];
			pointList.add (point);
			if (line.size() != dim) throw new IOException ("Dimension of point in starting file doesn't match with given demogrpahy.");
			for (int i=0; i<line.size(); i++) {
				point[i] = Double.parseDouble (line.get(i));
			}
		}
		
		// check size
		if (pointList.size() <= 1) {
			throw new IOException ("Specify at least two points in the start file, or use command line option to specify one point.");
		}
		
		// return it
		return pointList;
	}
}
