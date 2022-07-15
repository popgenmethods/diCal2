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

package edu.berkeley.diCal2.utility;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class RealPartition {

	/// get some balanced partition, where the exponential weight (with rate one) of each interval is equal [this makes it some how logarithmic]
	public static Interval[] getBalancedPartition (int nrIntervals) {
		assert (nrIntervals > 0);
		// get some place to store it
		List<Interval> returnList = new ArrayList<Interval>();
		
		
		if (nrIntervals > 1){
			// helpful
			double weightPerInterval = 1.0/nrIntervals;
			
			// and fill it (reverse order)
			double endPoint = Double.POSITIVE_INFINITY;
			double startPoint;
			
			// first one has to be done like this
			startPoint = - Math.log (weightPerInterval);
			returnList.add (new Interval (startPoint, endPoint));
			
			// then the rest
			for (int i=2; i<nrIntervals; i++) {
				// update endPoint to startPoint of previous interval [actually the interval to the right]
				endPoint = startPoint;
				// calculate new startPoint
				startPoint = - Math.log (i * weightPerInterval);
				// and add it to the list
				returnList.add (new Interval (startPoint, endPoint));
			}
			
			// also make the last one special
			endPoint = startPoint;
			startPoint = 0.0;
			returnList.add (new Interval (startPoint, endPoint));
	
			// give it away now
			// but be careful, its is in reverse order
			Collections.reverse (returnList);
		} else {
			returnList.add(new Interval(0d, Double.POSITIVE_INFINITY));
		}
		
		// now should be in right order
		return returnList.toArray (new Interval[0]);
	}
	
	/// boundary points of an interval [subclass]
	public static class Interval {
		
		// constructors
		public Interval(double startPoint, double endPoint) {
			this.startPoint = startPoint;
			this.endPoint = endPoint;
		}
		
		
		// just for debug
		public String toString () {
			return "[" + this.startPoint + ", " + this.endPoint + "]";
		}
		
		// the values
		public final double startPoint;
		public final double endPoint;
	}

	public static Interval[] getExponentialPartition(int nrIntervals, double expRate) {
		// get one for rate 1
		Interval[] preIntervals = getBalancedPartition (nrIntervals);
		
		// now transform the interval points
		Interval[] postIntervals = new Interval[preIntervals.length];
		for (int i=0; i<preIntervals.length; i++) {
			postIntervals[i] = new Interval (preIntervals[i].startPoint / expRate, preIntervals[i].endPoint / expRate);
		}

		// return it
		return postIntervals;
	}

	public static boolean allBiggerThanEpsilon (Interval[] intervals, double epsilon) {
		
		if (intervals == null) return true;
		
		for (Interval now : intervals) {
			if (now.endPoint - now.startPoint < epsilon) return false;
		}
		return true;
	}
	
	public static boolean isConsecutive (Interval[] intervals) {
		
		if (intervals == null) return true;
		
		double lastEnd = 0d;

		for (Interval now : intervals) {
			if (now.startPoint != lastEnd) return false;
			lastEnd = now.endPoint;
		}
		return (intervals[intervals.length-1].endPoint == Double.POSITIVE_INFINITY);
	}
}
