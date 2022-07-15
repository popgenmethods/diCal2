import msprime
import random

random.seed (4711)

def pieceWiseDemography():
	# here are some parameters
	# pop sizes
	# msprime convention is size at the more recent end counts
	N = [9771, 8000, 10000]
	# this corresponds to exp rate of 2
	growthRate = 0.0001
	# times
	# generation_time = 25
	t = [2000, 8000]

	# population configurations
	population_configurations = [
		msprime.PopulationConfiguration(
			sample_size=10, initial_size=N[0], growth_rate=growthRate)
	]
	# some size changes
	# AND IMPORTANTLY, THE GROWTH RATES ARE SET TO ZERO
	demographic_events = [
		msprime.PopulationParametersChange(
			time=t[i], initial_size=N[i+1], growth_rate=0, population_id=0)
		for i in range(len(t))
	]

	return (population_configurations, demographic_events)


def simulate (pop, demo, seqLen, mutRate, recoRate, seeed):
	tree_sequence = msprime.simulate(
		random_seed=seeed,
		Ne=10000,
		length=seqLen,
		mutation_rate=mutRate,
		recombination_rate=recoRate,
		population_configurations=pop,
		demographic_events=demo)

	return tree_sequence


def main():

	# sequence length
	contigLen = 100e6

	# make it
	(pop, demo) = pieceWiseDemography()

	# # print it
	# dd = msprime.DemographyDebugger(
	# 	population_configurations=pop,
	# 	# migration_matrix=mig,
	# 	demographic_events=demo)
	# dd.print_history()

	# simulate a sample
	mu = 0.0000000125
	r = mu

	# write the sample
	# 16 contigs
	for d in range(2):
		print ("[CONTIG_%d]" % d)
		tree_sequence = simulate (pop, demo, contigLen, mu, r, random.randint(0,99999999))
		ofs = open ("contig.%d.vcf" % d, "w")
		tree_sequence.write_vcf (ofs, 2)
		ofs.close()

	print ("[REFERENCE]")
	# need a reference sequences
	refSeq = ["A"]*int(contigLen)

	ofs = open ("output.ref", "w")
	ofs.write ("".join(refSeq) + "\n")
	ofs.close()


if __name__ == "__main__":
	main()