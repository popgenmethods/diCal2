import msprime
import random

random.seed (4724)

def migStopDemography():
	# here are some parameters
	# pop sizes
	N_1 = 10000
	N_2 = 10000
	N_3 = 10000
	N_12 = 10000
	N_A = 10000
	# times
	tDiv12 = 4000
	tDiv3 = 8000
	# introgression/admixture	
	tAdm = 1000
	admFraction = 0.03

	# population configurations
	population_configurations = [
		msprime.PopulationConfiguration(
			sample_size=4, initial_size=N_1),
		msprime.PopulationConfiguration(
			sample_size=4, initial_size=N_2),
		msprime.PopulationConfiguration(
			sample_size=4, initial_size=N_3)
	]
	migration_matrix = [
		[0,    0,    0],
		[0,    0,    0],
		[0,    0,    0]
	]
	demographic_events = [
		# admixture
		msprime.MassMigration(
			time=tAdm, source=0, destination=2, proportion=admFraction),
		# first merger
		msprime.MassMigration(
			time=tDiv12, source=1, destination=0, proportion=1.0),
		msprime.PopulationParametersChange(
			time=tDiv12, initial_size=N_12, growth_rate=0, population_id=0),
		msprime.PopulationParametersChange(
			time=tDiv12, initial_size=N_3, growth_rate=0, population_id=2),
		# second merger
		msprime.MassMigration(
			time=tDiv3, source=2, destination=0, proportion=1.0),
		msprime.PopulationParametersChange(
			time=tDiv3, initial_size=N_A, growth_rate=0, population_id=0)
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