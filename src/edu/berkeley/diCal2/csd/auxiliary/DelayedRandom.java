package edu.berkeley.diCal2.csd.auxiliary;

import java.util.Random;

public class DelayedRandom {
	
	private final Random myRandom;

	public DelayedRandom (Long seed) {
		if (seed == null) {
			// we do not have a seed, so we cannot use this RNG
			this.myRandom = null;
		}
		else {
			// we have a seed, so initialize the RNG properly
			this.myRandom = new Random (seed.longValue());
		}
	}
	
	public boolean isProper () {
		return (this.myRandom != null);
	}
	
	public Random getInternalRandom () {
		// check if proper
		if (this.isProper() ) {
			return this.myRandom;
		}
		else {
			this.fail();
			assert (false);
			return null;
		}
	}

	private void fail() {
		System.err.println("The specified analysis involves randomness. Please specifiy a seed with --seed.");
		System.exit(1);
	}

	public DelayedRandom spawnOffspring() {
		// spawn a delayed random offspring
		// that guy can also be improper
		if (this.isProper()) {
			return new DelayedRandom (new Long (this.myRandom.nextLong()));
		}
		else {
			return new DelayedRandom (null);
		}
	}

	public double nextDouble() {
		if (this.isProper()) {
			return this.myRandom.nextDouble();
		}
		else {
			this.fail();
			assert (false);
			return Double.NaN;
		}
	}

	public boolean nextBoolean() {
		if (this.isProper()) {
			return this.myRandom.nextBoolean();
		}
		else {
			this.fail();
			assert (false);
			return false;
		}
	}

	public double nextGaussian() {
		if (this.isProper()) {
			return this.myRandom.nextGaussian();
		}
		else {
			this.fail();
			assert (false);
			return Double.NaN;
		}
	}

	public int nextInt(int upper) {
		if (this.isProper()) {
			return this.myRandom.nextInt(upper);
		}
		else {
			this.fail();
			assert (false);
			return -1;
		}
	}
}
