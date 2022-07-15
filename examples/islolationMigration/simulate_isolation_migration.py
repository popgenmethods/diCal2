import msprime
import random

random.seed (4725)

def migStopDemography():
	# here are some parameters
	# pop sizes
	N_1 = 5000
	N_2 = 5000
	N_A = 10000
	# times
	tDivInGen = 2000
	# Migration rate.
	m = 1e-4

	# population configurations
	population_configurations = [
		msprime.PopulationConfiguration(
			sample_size=4, initial_size=N_1),
		msprime.PopulationConfiguration(
			sample_size=4, initial_size=N_2)
	]
	migration_matrix = [
		[0,    m],
		[m,    0]
	]
	demographic_events = [
		# merge both populations
		msprime.MassMigration(
			time=tDivInGen, source=1, destination=0, proportion=1.0),
		msprime.PopulationParametersChange(
			time=tDivInGen, initial_size=N_A, growth_rate=0, population_id=0),
		# stop migration
        msprime.MigrationRateChange(time=tDivInGen, rate=0, matrix_index=(0, 1)),
        msprime.MigrationRateChange(time=tDivInGen, rate=0, matrix_index=(1, 0))
	]

	return (population_configurations, migration_matrix, demographic_events)


def simulate (pop, mig, demo, seqLen, mutRate, recoRate, randomSeed):
	tree_sequence = msprime.simulate(
		random_seed=randomSeed,
		Ne=10000,
		length=seqLen,
		mutation_rate=mutRate,
		recombination_rate=recoRate,
		population_configurations=pop,
		migration_matrix=mig,
		demographic_events=demo)

	return tree_sequence



def main():

	# sequence length
	contigLen = 100e6

	# make it
	(pop, mig, demo) = migStopDemography()

	# # print it
	# dd = msprime.DemographyDebugger(
	# 	population_configurations=pop,
	# 	migration_matrix=mig,
	# 	demographic_events=demo)
	# dd.print_history()

	# # simulate a sample
	mu = 0.0000000125
	r = mu

	# write the sample
	# 16 contigs
	for d in range(2):
		print ("[CONTIG_%d]" % d)
		tree_sequence = simulate (pop, mig, demo, contigLen, mu, r, random.randint(0,99999999))
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
