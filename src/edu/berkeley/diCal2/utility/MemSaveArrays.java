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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MemSaveArrays {
	
	public static abstract class BigArrayShort {
		
		// the max size
		static long SINGLE_ARRAY_MAX_SIZE = Integer.MAX_VALUE - 10;

		// get the right array
		public static BigArrayShort getArray (long size) {
			if (size <= SINGLE_ARRAY_MAX_SIZE) {
				return new BigArrayShort.BigArrayShortSingle (size);
			}
			else {
				return new BigArrayShort.BigArrayShortMulti (size);
			}
		}

		// value access
		public abstract short get (long index);
		public abstract void set (long index, short value);
		public abstract void adjust (long index, short adjustment);

		// implementation with single array
		static class BigArrayShortSingle extends BigArrayShort {

			private short[] data;
			
			// constructor just makes an array
			public BigArrayShortSingle (long size) {
				assert (size <= SINGLE_ARRAY_MAX_SIZE);
				System.out.println ("# single array of size " + size);
				
				this.data = new short[(int) size];
			}
			
			@Override
			public short get (long index) {
				return this.data[(int) index];
			}

			@Override
			public void set (long index, short value) {
				this.data[(int) index] = value;
			}

			@Override
			public void adjust(long index, short adjustment) {
				this.data[(int) index] += adjustment;
			}
		}
		
		// implementation with multiple arrays
		static class BigArrayShortMulti extends BigArrayShort {

			private List<short[]> data;
			
			public BigArrayShortMulti (long size) {
				// prepare the space
				this.data = new ArrayList<short[]>();
				// running variable
				long runningSize = size;
				// loop
				while (runningSize > SINGLE_ARRAY_MAX_SIZE) {
					// add container
					this.data.add (new short[(int) SINGLE_ARRAY_MAX_SIZE]);
					// and decrease 
					runningSize -= SINGLE_ARRAY_MAX_SIZE;
				}
				assert (runningSize <= SINGLE_ARRAY_MAX_SIZE);
				// rest
				if (runningSize > 0) {
					// last one
					this.data.add (new short[(int) runningSize]);
				}
				System.out.println ("# multi array: " + this.data.size() + " arrays of size " + SINGLE_ARRAY_MAX_SIZE + " (or less)");				
			}
			
			public static int metaArrayIdx (long index) {
				return (int) (index / SINGLE_ARRAY_MAX_SIZE);
			}
			
			public static int withinIdx (long index) {
				return (int) (index % SINGLE_ARRAY_MAX_SIZE);
			}
			
			@Override
			public short get (long index) {
				return this.data.get(metaArrayIdx (index))[withinIdx (index)];
			}

			@Override
			public void set (long index, short value) {
				this.data.get(metaArrayIdx (index))[withinIdx (index)] = value;
			}

			@Override
			public void adjust (long index, short adjustment) {
				this.data.get(metaArrayIdx (index))[withinIdx (index)] += adjustment;
			}
		}
	}
	
	public static class FiveDimShortArray {
		
		// storage for data
		private final BigArrayShort data;
		
		// also store the dimensions
		private final long firstSize;
		private final long secondSize;
		private final long thirdSize;
		private final long fourthSize;
		private final long fifthSize;
		
		// get them dimensions
		public long getFirstSize() {
			return firstSize;
		}

		public long getSecondSize() {
			return secondSize;
		}

		public long getThirdSize() {
			return thirdSize;
		}

		public long getFourthSize() {
			return fourthSize;
		}

		public long getFifthSize() {
			return fifthSize;
		}

		/// constructor takes five dimensions
		public FiveDimShortArray (long firstSize, long secondSize, long thirdSize, long fourthSize, long fifthSize) throws IOException {
			// remember sizes
			this.firstSize = firstSize;
			this.secondSize = secondSize;
			this.thirdSize = thirdSize;
			this.fourthSize = fourthSize;
			this.fifthSize = fifthSize;
			
			// java does not allow to have arrays that are than Integer.MAX_VALUE
			// (with some buffer, 4 according to stackoverflow)
			long number = this.firstSize * this.secondSize * this.thirdSize * this.fourthSize * this.fifthSize;
			assert (number >= 0);
			System.out.println ("# 5-D array: (" + this.firstSize + "x" + this.secondSize + "x" + this.thirdSize + "x" + this.fourthSize + "x" + this.fifthSize + " = " + number + ")");
			
			// and create and array that's large enough
			this.data = BigArrayShort.getArray (number);
		}
		
		// getter and setter
		public short get (long firstIdx, long secondIdx, long thirdIdx, long fourthIdx, long fifthIdx) {
			return this.data.get (this.computeIndex (firstIdx, secondIdx, thirdIdx, fourthIdx, fifthIdx));
		}
		
		public void set (long firstIdx, long secondIdx, long thirdIdx, long fourthIdx, long fifthIdx, short value) {
			this.data.set (this.computeIndex (firstIdx, secondIdx, thirdIdx, fourthIdx, fifthIdx), value);
		}

		public void adjust (long firstIdx, long secondIdx, long thirdIdx, long fourthIdx, long fifthIdx, short adjustment) {
			this.data.adjust (this.computeIndex (firstIdx, secondIdx, thirdIdx, fourthIdx, fifthIdx), adjustment);
		}

		// the index conversion
		private long computeIndex (long firstIdx, long secondIdx, long thirdIdx, long fourthIdx, long fifthIdx) {
			long returnIdx = (((firstIdx * this.secondSize + secondIdx) * this.thirdSize + thirdIdx) * this.fourthSize + fourthIdx) * this.fifthSize + fifthIdx;
			assert (returnIdx >= 0);
			return returnIdx;
		}
	}
	
	public static class ThreeDimDoubleArray {
		
		// storage for data
		private final double[] data;
		
		// also store the dimensions
		private final int firstSize;
		private final int secondSize;
		private final int thirdSize;
		
		/// constructor takes five dimensions
		public ThreeDimDoubleArray (int firstSize, int secondSize, int thirdSize) {
			// remember sizes
			this.firstSize = firstSize;
			this.secondSize = secondSize;
			this.thirdSize = thirdSize;
			
			// and create and array that's large enough
			this.data = new double [this.firstSize * this.secondSize * this.thirdSize];
		}
		
		// getter and setter
		public double get (int firstIdx, int secondIdx, int thirdIdx) {
			return this.data[this.computeIndex (firstIdx, secondIdx, thirdIdx)];
		}
		
		public void set (int firstIdx, int secondIdx, int thirdIdx, double value) {
			this.data[this.computeIndex (firstIdx, secondIdx, thirdIdx)] = value;
		}

		public void adjust (int firstIdx, int secondIdx, int thirdIdx, double adjustment) {
			this.data[this.computeIndex (firstIdx, secondIdx, thirdIdx)] += adjustment;
		}

		// the index conversion
		private int computeIndex (int firstIdx, int secondIdx, int thirdIdx) {
			return (firstIdx * this.secondSize + secondIdx) * this.thirdSize + thirdIdx;
		}
	}

	public static class TwoDimDoubleArray {
		
		// storage for data
		private final double[] data;
		
		// also store the dimensions
		private final int firstSize;
		private final int secondSize;
		
		/// constructor takes five dimensions
		public TwoDimDoubleArray (int firstSize, int secondSize) {
			// remember sizes
			this.firstSize = firstSize;
			this.secondSize = secondSize;
			
			// and create and array that's large enough
			this.data = new double [this.firstSize * this.secondSize];
		}
		
		public int getFirstSize () {
			return this.firstSize;
		}
		
		public int getSecondSize () {
			return this.secondSize;
		}
		
		// getter and setter
		public double get (int firstIdx, int secondIdx) {
			return this.data[this.computeIndex (firstIdx, secondIdx)];
		}
		
		public void set (int firstIdx, int secondIdx, double value) {
			this.data[this.computeIndex (firstIdx, secondIdx)] = value;
		}

		public void adjust (int firstIdx, int secondIdx, double adjustment) {
			this.data[this.computeIndex (firstIdx, secondIdx)] += adjustment;
		}

		// the index conversion
		private int computeIndex (int firstIdx, int secondIdx) {
			return firstIdx * this.secondSize + secondIdx;
		}

		public void setRow (int firstIdx, double[] values) {
			assert (values.length == secondSize);

			// now go through and store
			for (int i=0; i<values.length; i++) {
				this.data[this.computeIndex (firstIdx, i)] = values[i];
			}
		}
	}

	public static class TwoDimFloatArray {
		
		// storage for data
		private final float[] data;
		
		// also store the dimensions
		private final int firstSize;
		private final int secondSize;
		
		/// constructor takes five dimensions
		public TwoDimFloatArray (int firstSize, int secondSize) {
			// remember sizes
			this.firstSize = firstSize;
			this.secondSize = secondSize;
			
			// and create and array that's large enough
			this.data = new float [this.firstSize * this.secondSize];
		}
		
		public int getFirstSize () {
			return this.firstSize;
		}
		
		public int getSecondSize () {
			return this.secondSize;
		}
		
		// getter and setter
		public float get (int firstIdx, int secondIdx) {
			return this.data[this.computeIndex (firstIdx, secondIdx)];
		}
		
		public void set (int firstIdx, int secondIdx, float value) {
			this.data[this.computeIndex (firstIdx, secondIdx)] = value;
		}

		public void adjust (int firstIdx, int secondIdx, float adjustment) {
			this.data[this.computeIndex (firstIdx, secondIdx)] += adjustment;
		}

		// the index conversion
		private int computeIndex (int firstIdx, int secondIdx) {
			return firstIdx * this.secondSize + secondIdx;
		}

		public void setRow (int firstIdx, float[] values) {
			assert (values.length == secondSize);

			// now go through and store
			for (int i=0; i<values.length; i++) {
				this.data[this.computeIndex (firstIdx, i)] = values[i];
			}
		}
	}
}
