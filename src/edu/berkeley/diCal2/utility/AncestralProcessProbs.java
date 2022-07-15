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

import java.math.BigDecimal;
import java.math.MathContext;

public class AncestralProcessProbs {
	
	private static BigDecimal factorial(int n) {
		return fallingFactorial(n, n);
	}
	
	private static BigDecimal risingFactorial(int n, int k) {
		BigDecimal r = BigDecimal.ONE;
		for (; k >= 1; k--) {
			r = r.multiply(BigDecimal.valueOf(n++));
		}
		return r;
	}
	
	private static BigDecimal fallingFactorial(int n, int k) {
		BigDecimal r = BigDecimal.ONE;
		for (; k >= 1; k--) {
			r = r.multiply(BigDecimal.valueOf(n--));
		}
		return r;
	}
	
	
	/**
	 * Returns the probability of there being j ancestors at time t given a 
	 * sample of size n taken at time zero.
	 * @param n Number of lineages at time 0.
	 * @param j Number of lineages at time t.
	 * @param t
	 * @return Probability.
	 */
	public static double getProb(int n, int j, double t) {
		// Exact probability computation
		MathContext mc = new MathContext(100);
		BigDecimal r = new BigDecimal(0, mc);
		BigDecimal s;
		for (int k=j; k <= n; k++) {
			// BigDecimal won't do fractional powers. Instead, write pow as 
			// pow = int(pow) + dec(pow) and multiply.
			double pow = t * k * (k-1) / 2;
			int ipow = (int)Math.floor(pow);
			double dpow = pow - ipow;
			s = BigDecimal.valueOf(Math.pow(Math.E, -dpow));
			s = s.multiply(BigDecimal.valueOf(-1).pow(k-j, mc));
			s = s.multiply(BigDecimal.valueOf(2*k - 1), mc);
			s = s.multiply(risingFactorial(j, k-1), mc);
			s = s.multiply(fallingFactorial(n, k), mc);
			s = s.divide(BigDecimal.valueOf(Math.E).pow(ipow), mc);
			s = s.divide(factorial(j), mc);
			s = s.divide(factorial(k - j), mc);
			s = s.divide(risingFactorial(n, k), mc);
			r = r.add(s, mc);
		}
		
		return r.doubleValue();
	}
}