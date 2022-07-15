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

public final class LogSum {
	private final double[] logSummandArray;
	private int currSize;
	
	private double maxLogSummand;
	
	public LogSum(int maxEntries) {
		this.logSummandArray = new double[maxEntries];
		
		reset();
	}
	
	public final void reset() {
		this.currSize = 0;
		this.maxLogSummand = Double.NEGATIVE_INFINITY;
	}
	
	public synchronized final void addLogSummand(double logSummand) {
		this.logSummandArray[currSize++] = logSummand;
		this.maxLogSummand = Math.max(this.maxLogSummand, logSummand);
	}
	
	public final double retrieveLogSum() {
		if (maxLogSummand == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
		
		double factorSum = 0;
		for (int i = 0; i < currSize; i++) {
			factorSum += Math.exp(logSummandArray[i] - maxLogSummand);
		}
		
		return Math.log(factorSum) + maxLogSummand;
	}
	
	public final static double computePairLogSum(double ls1, double ls2) {
		double maxLogSummand = Math.max(ls1, ls2);
		if (maxLogSummand == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
		
		double factorSum = Math.exp(ls1 - maxLogSummand) + Math.exp(ls2 - maxLogSummand);
		return Math.log(factorSum) + maxLogSummand;
	}
	
	public final static double computePairLogDifference(double ls1, double ls2) {
		double maxLogSummand = Math.max(ls1, ls2);
		if (maxLogSummand == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
		
		double factorSum = Math.exp(ls1 - maxLogSummand) - Math.exp(ls2 - maxLogSummand);
		return Math.log(factorSum) + maxLogSummand;
	}
}