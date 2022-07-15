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

import org.apache.commons.math3.ode.nonstiff.HighamHall54Integrator;

public class RobustHighamHall extends HighamHall54Integrator {

	public RobustHighamHall(double minStep, double maxStep,
			double scalAbsoluteTolerance, double scalRelativeTolerance) {
		super(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
	}

	public RobustHighamHall(double minStep, double maxStep,
			double[] vecAbsoluteTolerance, double[] vecRelativeTolerance) {
		super(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
	}

	@Override
	// override this to make it more robust
	public double initializeStep(boolean forward, int order, double[] scale,
			double t0, double[] y0, double[] yDot0, double[] y1, double[] yDot1) {
		double toReturn = super.initializeStep(forward, order, scale, t0, y0, yDot0, y1, yDot1);
		if (Double.isNaN(toReturn)) {
			toReturn = getMaxStep();
		}
		return toReturn;
	}
}
