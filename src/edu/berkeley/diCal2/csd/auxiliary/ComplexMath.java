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

package edu.berkeley.diCal2.csd.auxiliary;

import Jampack.Z;

public class ComplexMath {
	/// do the exponential
	public static Z exponential (Z z) {
		double expRe = Math.exp(z.re);
		return new Z (expRe*Math.cos (z.im), expRe*Math.sin (z.im));
	}
	
	// new (better) times [convention: 0 \times \infty = 0]
	public static Z times (double a, Z b) {
		// any of the factors zero?
		if ((Math.abs (a) < EPSILON) || (Z.abs(b) < EPSILON)) {
			// then the product is zero
			return new Z(0d, 0d);
		}
		else {
			// normal product
			return (new Z()).Times(a, b);
		}
	}
	
	// new (better) times [convention: 0 \times \infty = 0]
	public static Z ourTimes (Z a, Z b) {
		// any of the factors zero?
		if ((Z.abs (a) < EPSILON) || (Z.abs(b) < EPSILON)) {
			// then the product is zero
			return new Z(0d, 0d);
		}
		else {
			// normal product
			return (new Z()).Times(a, b);
		}
	}

	// the EPSILON
	public static double EPSILON = 1e-15;
}
