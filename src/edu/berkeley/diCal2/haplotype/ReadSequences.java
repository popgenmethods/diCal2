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

package edu.berkeley.diCal2.haplotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.berkeley.diCal2.haplotype.FSAXTypeConfig.FSAReference;
import edu.berkeley.diCal2.haplotype.HFSAFullHaplotype.HFSAFullHapFactory;
import edu.berkeley.diCal2.utility.RealPartition.Interval;
import gnu.trove.list.TByteList;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.hash.TByteHashSet;
import gnu.trove.set.hash.TIntHashSet;

public class ReadSequences {

	public static final int VCF_NUMBER_DEFAULT_COLUMNS = 9;
	public static final int MISSING_ALLELE = -1;
//	public static final int SEG_SITE_IN_REFERENCE = Integer.MIN_VALUE;
	public static final int SEG_SITE_IN_REFERENCE = Byte.MIN_VALUE;
	public static final int REF_COLUMN = 3;
	public static final int ALT_COLUMN = 4;
	public static final int SNP_COLUMN = 1;
	public static final int FILTER_COLUMN = 6;
	private static final int FORMAT_COLUMN = 8;
	
	public enum ReferenceMode {
		FIXED, RANDOM, MISSING, NONE
	}

	public static class MutableInt {
		private int value;
		public MutableInt (int value) {
			this.value = value;
		}
		public void increase () {
			this.value++;
		}
		public int getValue () {
			return this.value;
		}
	}
	
	// reads a file where each line is a full sequence, then compresses it into a list of HFSAXFullHaplotype (which only store the haplotypes at the SNPs)
	public static List<HFSAXFullHaplotype> xReadFullSequenceFormat (Reader reader, List<Boolean> sequencesToRead, int numAlleles, char[] alleles, char[] commentCharacters, char[] missingCharacters) throws IOException {
		
		List<HFSAFullHaplotype> nonXhapsList = readSequences(new HFSAFullHapFactory(), new BufferedReader (reader), sequencesToRead, numAlleles, alleles, commentCharacters, missingCharacters);
		
		// make the raw reference
		int[] ref = nonXhapsList.get(0).getRawAlleleConfig();
		
		// find the snps
		for (HFSAFullHaplotype hap : nonXhapsList) {
			for (int l =0; l< hap.getNumLoci(); l++){
				if (hap.getAllele(l) != ref[l]) {
					ref[l] = ReadSequences.SEG_SITE_IN_REFERENCE;
				}
			}	
		}
		
		// store and order the snps
		TIntIntHashMap referenceIdxMap = new TIntIntHashMap();
		for (int l=0, i=0; l<ref.length; l++){
			if (ref[l] == ReadSequences.SEG_SITE_IN_REFERENCE) {
				referenceIdxMap.put(l, i++);
			}
		}
		
		// make the reference
		SimpleFSARef fsaReference = new SimpleFSARef(ref, referenceIdxMap);
		
		List<HFSAXFullHaplotype> xList = new ArrayList<HFSAXFullHaplotype>();
		for (HFSAFullHaplotype nonXhap : nonXhapsList) {
			
			// maybe some empty ones
			if (nonXhap == null) {
				xList.add (null);
				continue;
			}
			
			// make the x haplotype
			int[] alleleConfig = new int[fsaReference.getNumSegLoci()];
			for (int snpIdx = 0; snpIdx < fsaReference.getNumSegLoci(); snpIdx++) {
				alleleConfig[snpIdx] = nonXhap.getAllele(fsaReference.getSegLocus(snpIdx));				
			}
			HFSAXFullHaplotype xHap = new HFSAXFullHaplotype(fsaReference, alleleConfig);
			
			// add it
			xList.add(xHap);
		}
		
		return xList;
	}
	
	public static List<HFSAXFullHaplotype> xReadCompressedSequenceFormat (Reader reader, List<Boolean> sequencesToRead, int numAlleles, char[] alleles, char[] commentCharacters, char[] missingCharacters) throws IOException {
		// make comment characters
		TreeSet<Character> commentCharacterSet = new TreeSet<Character>();
		for (char c : commentCharacters) {
			commentCharacterSet.add(c);
		}
		
		// open the sequence file
		BufferedReader bSequenceReader = new BufferedReader(reader);
		
		int numLoci = -1;
		SimpleFSARef fsaRef = null;
		
		boolean readRef = false, readNumLoci = false, readSegSites = false, readHaps = false;

		String line = null;
		while ((line = bSequenceReader.readLine()) != null) {
			
			line = line.trim();
			
			if (line.equals("REFERENCE")) {
				readRef = readNumLoci = readSegSites = readHaps = false;
				readRef = true;
			} else if (line.equals("NUM LOCI")) {
				readRef = readNumLoci = readSegSites = readHaps = false;
				readNumLoci = true;
			} else if (line.equals("SEGREGATING SITES")) {
				readRef = readNumLoci = readSegSites = readHaps = false;
				readSegSites = true;
			} else if (line.equals("HAPLOTYPES")) {
				readRef = readNumLoci = readSegSites = readHaps = false;
				readHaps = true;
				break;
			} else if (readRef) {
				if (fsaRef != null) {
					throw new IOException ("Only one of REFERENCE or SEGREGATING SITES should be specified, and exactly once");
				}
				
				fsaRef = makeFSARef(line, numAlleles, alleles, missingCharacters);
				
				readRef = false;

			} else if (readNumLoci) {
				numLoci = Integer.parseInt(line);
				
				readNumLoci = false;
			} else if (readSegSites) {
				if (fsaRef != null) {
					throw new IOException ("Only one of REFERENCE or SEGREGATING SITES should be specified, and exactly once.");
				}
				
				if (numLoci <= 0) {
					throw new IOException ("NUM LOCI must be specified and positive before SEGREGATING SITES are read.");
				}
				
				fsaRef = makeFSARef(numLoci, line);
				
				readSegSites = false;
			}
		}
		
		if (readHaps == false || fsaRef == null) {
			throw new IOException ("Incorrectly formatted compressed sequence file (missing reference and haplotypes). Try --oldCompressedDataFormat.");
		}
		
		// create the haplotype factory
		FixedReferenceFSAXFactory<HFSAXFullHaplotype> hapFactory = new FixedReferenceFSAXFactory<HFSAXFullHaplotype>(new HFSAXFullHaplotype.HFSAXFullHapFactory(), fsaRef);
		
		return readSequences(hapFactory, bSequenceReader, sequencesToRead, numAlleles, alleles, commentCharacters, missingCharacters);
	}
	
	
	

	// function to read old format for compressed file sequence
	// reads a file where the first line is a list of the segregating sites, and each subsequent line is a haplotype at the SNPs
	// note that the non-segregating sites are always assumed to have allele 0
	public static List<HFSAXFullHaplotype> xReadOldCompressedSequenceFormat (Reader reader, List<Boolean> sequencesToRead, int numAlleles, char[] alleles, int numLoci, char[] commentCharacters, char[] missingCharacters) throws IOException {
		// make comment characters
		TreeSet<Character> commentCharacterSet = new TreeSet<Character>();
		for (char c : commentCharacters) {
			commentCharacterSet.add(c);
		}
		
		// open the sequence file
		BufferedReader bSequenceReader = new BufferedReader(reader);
		
		// read in the segregating sites on first line of sequence file
		String segSitesString = "";
		// skip blank and commented lines
		while (segSitesString.equals("") || commentCharacterSet.contains(segSitesString.charAt(0))) {			
			segSitesString = bSequenceReader.readLine().trim();		
		}
		
		// make the reference
		SimpleFSARef fsaReference = makeFSARef(numLoci, segSitesString);
		
		// make the haplotype factory
		FixedReferenceFSAXFactory<HFSAXFullHaplotype> hapFactory = new FixedReferenceFSAXFactory<HFSAXFullHaplotype>(new HFSAXFullHaplotype.HFSAXFullHapFactory(), fsaReference);
		
		return readSequences(hapFactory, bSequenceReader, sequencesToRead, numAlleles, alleles, commentCharacters, missingCharacters);
	}
	
	// segSitesString should be a tab-separated line of integers denoting the positions of segregating sites
	public static SimpleFSARef makeFSARef (int numLoci, String segSitesString) throws IOException {
		String[] fields = segSitesString.split("(\t| )+");
		
		// convert seg sites to integers
		List<Integer> segSites = new ArrayList<Integer>();
		for (String entry : fields) {
			int site = Integer.parseInt(entry);
			if (site >= numLoci) {
				throw new IOException("Index of segregating site in sequence-file is larger than sequence length (number of loci) provided.");
			}
			segSites.add(site);
		}
		
		// store and order the snps
		int[] ref = new int[numLoci];
		TIntIntHashMap referenceIdxMap = new TIntIntHashMap();
		int curSnp = 0;
		for (int site : segSites) {
			ref[site] = ReadSequences.SEG_SITE_IN_REFERENCE;
			referenceIdxMap.put(site, curSnp++);
		}
		
		// make the reference
		return new SimpleFSARef(ref, referenceIdxMap);
	}
	
	// somehow this reference might not be allowed to have missing characters in it. however, we need that for vcfs
	private static SimpleFSARef makeFSARef(String referenceString, int numAlleles, char[] alleles, char[] missingCharacters) {
		TreeSet<Character> missingCharSet = new TreeSet<Character>();
		for (char c : missingCharacters) missingCharSet.add(c);
		assert !missingCharSet.contains('x');
		
		// get rid of some things
		referenceString = referenceString.trim();
		
		int numLoci = referenceString.length();

		// and copy into int array
		int[] ref = new int[numLoci];
		TIntIntHashMap referenceIdxMap = new TIntIntHashMap();
		int curSnp = 0;		
		for (int l = 0; l < numLoci; l++) {
			
			if (referenceString.charAt(l) == 'x') {
				ref[l] = ReadSequences.SEG_SITE_IN_REFERENCE;
				referenceIdxMap.put(l, curSnp++);
			} else if (missingCharSet.contains(referenceString.charAt(l))) {
				ref[l] = ReadSequences.MISSING_ALLELE;
			} else {
				if (alleles != null) {
					for (int a = 0; a < alleles.length; a++) {
						if (referenceString.charAt(l) == alleles[a]) {
							ref[l] = a;
						}
					}
				} else {
					ref[l] = referenceString.charAt(l) - '0';
				}
				
				assert (ref[l] >= 0);
				assert (ref[l] < numAlleles);				
			}
			
		}
		return new SimpleFSARef(ref, referenceIdxMap);
	}
	
	public static List<int[]> byteArrayListToIntArrayList (List<byte[]> byteArrayList) {
		List<int[]> intArrayList = new ArrayList<int[]>();
		
		for (byte[] byteArray : byteArrayList) {
			intArrayList.add(SimpleFSARef.byteArrayToIntArray(byteArray));
		}
		
		return intArrayList;
	}
	
	public static <H extends FSAHaplotype> List<H> readSequences (FSAHapFactory<H> hapFactory, BufferedReader bufferedSequenceReader, List<Boolean> sequencesToRead, int numAlleles, char[] alleles, char[] commentCharacters, char[] missingCharacters) throws IOException {
		return makeSequences(byteArrayListToIntArrayList(readRawSequences(bufferedSequenceReader, sequencesToRead, numAlleles, alleles, commentCharacters, missingCharacters, false)), hapFactory);
	}
	
	@SuppressWarnings("unchecked")
	public static <H extends FSAHaplotype> List<H> makeSequences (List<int[]> rawSequences, FSAHapFactory<H> hapFactory) throws IOException {
		List<H> sequenceList = new ArrayList<H>();
		for (int[] rawSequence : rawSequences) {
			
			// account for empty ones
			if (rawSequence == null) {
				sequenceList.add (null);
				continue;
			}
			
			// check here
			if ((hapFactory instanceof FixedReferenceFSAXFactory) && (((FixedReferenceFSAXFactory<HFSAXFullHaplotype>)hapFactory).getNumSegSites() != rawSequence.length)) {
				throw new IOException("List of alleles for haplotype doesn't agree with number of segregating sites in reference.");
			}
			// make haplotype
			sequenceList.add(hapFactory.getHaplotype(rawSequence));
		}
		return sequenceList;
	}
	
	public static List<byte[]> readRawSequences (BufferedReader bufferedSequenceReader, List<Boolean> sequencesToRead, int numAlleles, char[] alleles, char[] commentCharacters, char[] missingCharacters, boolean readRef) throws IOException {
		
		assert (!readRef || (sequencesToRead == null));

		Set<Character> missingCharSet = new TreeSet<Character>();
		for (char c : missingCharacters) missingCharSet.add(c);
		
		// make comment characters
		Set<Character> commentCharacterSet = new TreeSet<Character>();
		for (char c : commentCharacters) commentCharacterSet.add(c);
		
		// what we will return
		List<byte[]> returnAlleleConfigs = new ArrayList<byte[]>();		
		
		// read
		String line = null;
		int hapIdx = -1;
		int firstNonEmptyIdx = -1;
		while ((line = bufferedSequenceReader.readLine()) != null) {
			// always ignore comment and empty lines
			if (!line.trim().equals("") && !commentCharacterSet.contains(line.charAt(0))){

				hapIdx++;
				
				// we don't want some of them
				// if sequencesToRead == null, read all of the sequences
				if (sequencesToRead != null && hapIdx >= sequencesToRead.size()) {
					assert (!readRef);
					throw new IOException("More haplotypes in inputfile than in config file."); 
				}
				if ((sequencesToRead != null) && !sequencesToRead.get(hapIdx)) {
					// not this one
					returnAlleleConfigs.add (null);
					continue;
				}

				// remember the first non-empty
				if (firstNonEmptyIdx < 0) {
					firstNonEmptyIdx = hapIdx;
				}

				// get rid of nasty stuff
				String haplotypeString = line.trim();



				int numLoci = haplotypeString.length();

				// and copy into int array
				byte[] alleleConfig = new byte[numLoci];
				returnAlleleConfigs.add (alleleConfig);

				// assert all haplotypes have same number of loci
				if (!readRef) {
					if (alleleConfig.length != returnAlleleConfigs.get(firstNonEmptyIdx).length) {
						throw new IOException("two sequences in input file have unequal length");
					}
				}

				Arrays.fill(alleleConfig, (byte)ReadSequences.SEG_SITE_IN_REFERENCE);
				for (int l = 0; l < numLoci; l++) {

					alleleConfig[l] = (byte)getAlleleWithoutMap (Character.toUpperCase(haplotypeString.charAt(l)), alleles, missingCharSet);

					assert (alleleConfig[l] >= ReadSequences.MISSING_ALLELE);
					assert (alleleConfig[l] < numAlleles);				
				}
			}
		}

		// return what you got
		return returnAlleleConfigs;
	}

	private static int getAlleleWithoutMap (char c, char[] alleles, Set<Character> missingCharSet) {
		if (missingCharSet.contains (c)) {
			return ReadSequences.MISSING_ALLELE;
		} else {
			if (alleles != null) {
				for (int a = 0; a < alleles.length; a++) {
					if (c == alleles[a]) {
						return a;
					}
				}
				return ReadSequences.SEG_SITE_IN_REFERENCE;
			} else {
				return c - '0';
			}
		}
	}

	public static List<HFSAXFullHaplotype> readVcf (Reader reader, List<Boolean> haplotypesToRead, int numInternalAlleles, char[] commentCharacters, boolean acceptUnphasedAsMissing, String pathToVcf, int offset, Reader bedReader, boolean verbose, String filterPassString, String commandLineReferenceFileName, boolean vcfIgnoreDoubleEntries) throws IOException {
		
		// hard code this here for now since it is the vcf standard (or my version of it anyways)
		char[] alleles = new char[] {'A', 'C', 'G', 'T'};
		int numRealAlleles = alleles.length;
		char[] refAlleles = alleles;
		char[] realMissingChars = new char[] {'n', 'N', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V'};
		char[] indMissingChars = new char[] {'.'};
		
		int numHapsToRead = 0;
		for (Boolean thisOne : haplotypesToRead) {
			if (thisOne) numHapsToRead++;
		}
		System.out.println ("# total haps in config-file: " + haplotypesToRead.size());
		System.out.println ("# num haps we want: " + numHapsToRead);
		System.out.println ("# number of internal alleles we want: " + numInternalAlleles);
		System.out.println ("# vcf offset: " + offset);
		
		
		List<Interval> bedRegions = null;
		// read bed stuff, if we have it
		if (bedReader != null) {
			bedRegions = new ArrayList<Interval>();
			BufferedReader buffReader = new BufferedReader (bedReader);
			String line = buffReader.readLine();
			while(line != null){
				String[] lineArray = line.split("\\s+");
				if (lineArray.length != 3) {
					throw new IOException("A row in the bed-file doesn't have 3 columns.");
				}
				Interval toAdd = new Interval(Integer.parseInt(lineArray[1]), Integer.parseInt(lineArray[2]));
				bedRegions.add(toAdd);
				// make sure that the order works
				if (bedRegions.size() >= 2) {
					if (bedRegions.get(bedRegions.size()-2).endPoint > bedRegions.get(bedRegions.size()-1).startPoint) {
						throw new IOException("Error in bed-file: region ending at " + (int)bedRegions.get(bedRegions.size()-2).endPoint + " followed by region starting at " + (int)bedRegions.get(bedRegions.size()-1).startPoint);
					}
				}
				line = buffReader.readLine();
			}
			buffReader.close();	//be nice
		}

		// missing set
		// TODO apparently the vcf standard says this is only '.'
		TreeSet<Character> indMissingCharSet = new TreeSet<Character>();
		for (char c : indMissingChars) indMissingCharSet.add(c);
		// this can be more though
		TreeSet<Character> realMissingCharSet = new TreeSet<Character>();
		for (char c : realMissingChars) realMissingCharSet.add(c);
		// and add this for good measure
		for (char c : indMissingChars) realMissingCharSet.add(c);

		
		// make comment characters
		TreeSet<Character> commentCharacterSet = new TreeSet<Character>();
		for (char c : commentCharacters) commentCharacterSet.add(c);
		
		// and a character map for the actual alleles
		// TODO the vcf standard is more specific too
		Map<Character,Integer> alleleToIdx =  getAlleleToIdxMap (alleles, numRealAlleles);
		
		// open the sequence file
		BufferedReader bSequenceReader = new BufferedReader(reader);

		// need some buffer first, cause we don't know the segregating sites right away
		ArrayList<Integer> segSites = new ArrayList<Integer>();
		// [site][individuals]
		ArrayList<TByteList> variants = new ArrayList<TByteList>();
		// also need some reference, I guess
		byte[] referenceSequence = null;
		
		int filterNotPassed = 0;
		
		// go through every line 
		// read things
		String line = null;
		// I think we don't have to parse the header
		Integer numCol = null;
		int lineNumber = 0;
		int tooManyAlleles = 0;
		int currRegion = 0;
		int numWithMissing = 0;
		int nonSegSites = 0;
		int numVariationLines = 0;
		int numDoubleVcfEntries = 0;
		int previousPosForCheck = -1;
		
		// we need them here
		boolean expectReference;
		boolean expectHeader;
		boolean expectSNPs;
		// now figure out what to do with the reference
//		boolean useExternalReference = (commandLineReferenceFileName != null) || (referenceMode != ReferenceMode.NONE);
		boolean useExternalReference = (commandLineReferenceFileName != null);
		if (useExternalReference) {
			if (commandLineReferenceFileName != null) {

				// here comes the wrapper
				referenceSequence = handleReferenceFile  (false, commandLineReferenceFileName, bedRegions, refAlleles, commentCharacters, realMissingChars, verbose);
			}
//			else {
//				// all the internal ones need this
//
//				// how long?
//				referenceSequence = new int[numLoci];
//				
//				if (referenceMode == ReferenceMode.FIXED) {
//					// fixed allele at every locus
//					// use mod to make it a bit more interesting
//					for (int i=0; i<referenceSequence.length; i++) {
//						referenceSequence[i] = alleles[i % 4]; 
//					}
//				}
//				else if (referenceMode == ReferenceMode.RANDOM) {
//					// build a reference by sampling from the stationary distribution of the mutation process
//					// which unfortunately might not have 4 alleles
//					// TODO so for now we just use uniform, maybe make command line argument for other distribution in the future
//					// does not matter too much, since it is only at non-seg sites anyways
//					
//					// and random guy at every locus
//					for (int i=0; i<referenceSequence.length; i++) {
//						referenceSequence[i] = alleles[myRandom.nextInt (alleles.length)]; 
//					}
//					// that's it?
//				}
//				else if (referenceMode == ReferenceMode.MISSING) {
//					// just put missing everywhere
//					for (int i=0; i<referenceSequence.length; i++) {
//						referenceSequence[i] = 'N'; 
//					}
//				}
//				else {
//					assert (false);
//				}
//			}
			
			assert (referenceSequence != null);
			
			// reference should is not in the file
			expectReference = false;
			expectHeader = true;
			expectSNPs = false;
		}
		else {
			// reference should be in the file
			expectReference = true;
			expectHeader = false;
			expectSNPs = false;
		}
		
		// count the unphased genotypes
		MutableInt unphasedCount = new MutableInt(0);
		
		while ((line = bSequenceReader.readLine()) != null) {
			lineNumber++;
			
			// get rid white space
			line = line.trim();

			if (expectReference) {
				if (line.startsWith("##reference")) {
					// get the name of the reference file
					String[] fields = line.split ("file://");
					if (fields.length != 2) {
						throw new IOException ("Invalid url-protocol for reference in vcf file.");
					}
					// get the name
					String referenceFileName = fields[1];

					// if the path to the reference is absolute, we can just take it
					if (!(new File(referenceFileName)).isAbsolute()) {
						// otherwise, we need to adjust it
						referenceFileName = (new File(new File(pathToVcf), referenceFileName)).getAbsolutePath();
					}
					
					// here comes the wrapper
					referenceSequence = handleReferenceFile  (true, referenceFileName, bedRegions, refAlleles, commentCharacters, realMissingChars, verbose);
										
					// set flags
					expectReference = false;
					expectHeader = true;
				}
			}
			else if (expectHeader) {
				if (line.startsWith("#CHROM")) {
					// basically we only want to know how many columns for individuals to expect
					String[] fields = line.split ("\\s+");
					
					if (!fields[ReadSequences.REF_COLUMN].equals("REF")) {
						throw new IOException ("Invalid field name in vcf file: " + fields[ReadSequences.REF_COLUMN] + ". REF field expected in column " + ReadSequences.REF_COLUMN + ".");
					}
					if (!fields[ReadSequences.ALT_COLUMN].equals("ALT")) {
						throw new IOException ("Invalid field name in vcf file: " + fields[ReadSequences.ALT_COLUMN] + ". ALT field expected in column " + ReadSequences.ALT_COLUMN + ".");
					}
					if (!fields[ReadSequences.FILTER_COLUMN].equals("FILTER")) {
						throw new IOException ("Invalid field name in vcf file: " + fields[ReadSequences.FILTER_COLUMN] + ". FILTER field expected in column " + ReadSequences.FILTER_COLUMN + ".");
					}
					if (!fields[ReadSequences.FORMAT_COLUMN].equals("FORMAT")) {
						throw new IOException ("Invalid field name in vcf file. FORMAT field expected.");
					}
					// the number of individual columns should be all columns minus 9 
					numCol = fields.length - ReadSequences.VCF_NUMBER_DEFAULT_COLUMNS;
					if (verbose) System.out.println ("# VCF has " + numCol + " columns with data for individuals.");

					// no assertion on user input
//					assert (numCol <= numHapsToRead <= numCol*2);
					// we have to check later that this works
//					if (numCol*2 != haplotypesToRead.size()) throw new IOException ("Number of haplotypes differs between vcf file and config file (" + numCol*2 + ", " + haplotypesToRead.size() + ").");
						
					// set flags
					expectHeader = false;
					expectSNPs = true;
				}
			}
			else if (expectSNPs) {
				
				if (line.equals("") || commentCharacterSet.contains(line.charAt(0))) continue;
				
				// we have a real snp-line
				numVariationLines++;
				
				// keep 'em coming
				String[] fields = line.split("\\s+");
				if (fields.length != (ReadSequences.VCF_NUMBER_DEFAULT_COLUMNS + numCol)) {
					throw new IOException ("Line " + lineNumber + " in vcf file has invalid number of columns.");
				}
				
				// everything fine, so get the entries
				// first the pos [position is in the second field]
				int snpPos = Integer.valueOf(fields[ReadSequences.SNP_COLUMN]);
				// remember, vcf starts at 1
				int pos = snpPos - 1 - offset;
				if ((pos < 0) || (pos >= referenceSequence.length)) {
					throw new IOException ("Position in vcf (corrected for offset) is outside of reference sequence.");
				}
				
				// we cannot have two vcf entries for the same position
				if (pos == previousPosForCheck) {
					if (vcfIgnoreDoubleEntries) {
						// remember it
						numDoubleVcfEntries += 1;
						continue;
					}
					else {
						throw new IOException("More than one entry provided for position " + (pos+1+offset) + " on line " + lineNumber + " in vcf-file. Use --vcfIgnoreDoubleEntries to ignore mutliple entries.");
					}
				}
				// we don't use it again
				previousPosForCheck = pos;
				
				
//				if (pos % 1000000 == 0) {
//					MetaOptimization.synchronizedPrintln(String.valueOf(pos) + '\t' + String.valueOf(variants.size()) + '\t' + String.valueOf(segSites.size()));
//				}
				
				
				// if bed file says ignore this, we ignore it
				if (bedRegions != null) {
					while ((currRegion < bedRegions.size()) && (pos >= bedRegions.get(currRegion).endPoint)) {
						currRegion++;
					}
					
					if ((currRegion < bedRegions.size()) && (pos >= bedRegions.get(currRegion).startPoint)) {
						// don't take this (reference should already be cleared)
						assert (referenceSequence[pos] == ReadSequences.MISSING_ALLELE);
						
						continue;
					}
				}
				
				// then get the reference allele
				String refAllele = fields[ReadSequences.REF_COLUMN];
				if (refAllele.length() != 1) {
					throw new IOException ("No structural variation for reference allele allowed in VCF file (line " + lineNumber + "): " + fields[ReadSequences.REF_COLUMN]);
				}

				// get char in upper case
				char refChar = Character.toUpperCase(refAllele.charAt(0));
				
				if (!(alleleToIdx.containsKey(refChar) || realMissingCharSet.contains(refChar))) {
					throw new IOException ("Unexpected reference allele in VCF file (line " + lineNumber + "): " + fields[ReadSequences.REF_COLUMN]);
				}
				
				int referenceAllele = ReadSequences.MISSING_ALLELE;
				if (alleleToIdx.containsKey(refChar)) {
					referenceAllele = alleleToIdx.get(refChar);
				}
				
				
				// this one has to change
				// if we provide a reference, it has to match,
				// otherwise we have to set it
//				boolean modifyReference = useExternalReference && (commandLineReferenceFileName == null);
//				assert (!modifyReference || (referenceMode != ReferenceMode.NONE));
//				
//				if (modifyReference) {
//					// set the reference to match vcf
//					referenceSequence[pos] = referenceAllele;
//				}
//				else {
					// otherwise, it should already match
				if (!(referenceAllele == referenceSequence[pos])) {
					throw new IOException ("Reference allele provided in VCF file at position " + (pos+1+offset) + " is different from reference sequence (line " + lineNumber + ").");
				}
//				}
				
				//  and the alternative allele(s)
				int[] vcfEntryToAlleleIdx = null;
				String[] fieldings = fields[ReadSequences.ALT_COLUMN].split(",");
				if(fields[ReadSequences.FILTER_COLUMN].equals (filterPassString)){
					if (fieldings.length > numInternalAlleles-1) {
						// too many alleles
						tooManyAlleles++;
						
						// mark reference as missing
						referenceSequence[pos] = ReadSequences.MISSING_ALLELE;
						
						// and move on
						continue;
					}
					// one more for reference
					vcfEntryToAlleleIdx = new int[1 + fieldings.length];
					for (int i=0; i<fieldings.length; i++) {
						String alt = fieldings[i];
						if (alt.length() != 1) throw new IOException ("No structural variation allowed in VCF file (line " + lineNumber + ").");
						char altChar = Character.toUpperCase(alt.charAt(0));
						if (alleleToIdx.containsKey(altChar)) {
							if (numInternalAlleles == 4) {
								// +1 here important
								vcfEntryToAlleleIdx[i+1] = alleleToIdx.get(altChar);
							}
							else {
								// +1 here important
								vcfEntryToAlleleIdx[i+1] = i+1;
							}
						}
						else if (realMissingCharSet.contains(altChar)) {
							vcfEntryToAlleleIdx[i+1] = ReadSequences.MISSING_ALLELE;
						}
						else {
							throw new IOException ("Invalid alternative allele in VCF file at pos " + (pos+1) + "(line " + lineNumber + ").");
						}
					}
					
					// put the reference in
					if (numInternalAlleles == 4) {
						vcfEntryToAlleleIdx[0] = referenceAllele;
					}
					else if (realMissingCharSet.contains(refChar)) {
						vcfEntryToAlleleIdx[0] = ReadSequences.MISSING_ALLELE;
					}
					else {
						vcfEntryToAlleleIdx[0] = 0;
					}
					// now check that they are distinct
					if (vcfEntryToAlleleIdx.length > 1) {
						TIntHashSet mySet = new TIntHashSet (vcfEntryToAlleleIdx);
						if (mySet.size() != vcfEntryToAlleleIdx.length) {
							throw new IOException ("Reference and alternative alleles in VCF file at pos " + (pos+1) + " have to be different (line " + lineNumber + ").");
						}
					}
			
				} else {	//this site didn't pass some filter, so let's go ahead and replace this one with missing
					// filter not passed
					filterNotPassed += 1;
					
					// mark as missing and go to next line
					// mark reference as missing
					referenceSequence[pos] = ReadSequences.MISSING_ALLELE;
					continue;
//					// all missing
//					vcfEntryToAlleleIdx = null;
				}
				

				// make sure there is only genotypes in the columns for the individuals
				if (!fields[ReadSequences.FORMAT_COLUMN].equals ("GT")) {
					throw new IOException("FORMAT column in VCF file at pos " + (pos+1) + " is not equal to 'GT'. Please reformat your VCF-file such that only genotypes are provided.");
				}

				// then the actual variants
				TByteList currVariants = new TByteArrayList();

				int runningFirstHapIdx = 0;
				
				// and go through the columns (each column has one or two haps)
				for (int indVcf=0; indVcf<numCol; indVcf++) {
					
					// first read the entry from file to see how many we want
					int colIdx = ReadSequences.VCF_NUMBER_DEFAULT_COLUMNS + indVcf;
					String thisSNP = fields[colIdx];

					// how many?
					boolean twoHaps = true;
					int thisLen = thisSNP.length();
					if (thisLen == 1) {
						// only one hap
						twoHaps = false;
					}
					else if (thisLen == 3) {
						// two haps
						assert (twoHaps);
					}
					else {
						// not allowed
						throw new IOException("Genotype entry in vcf-file on line " + lineNumber + " in column " + colIdx + " is invalid.");
					}
		
					int firstHapIdx = runningFirstHapIdx;
					
					// more haplotypes in VCF than in config?
					if (firstHapIdx >= haplotypesToRead.size())  {
						throw new IOException("Genotype entry in vcf-file on line " + lineNumber + " has more haplotypes than specified in config-file (" + haplotypesToRead.size() + ")." );
					}

					if (twoHaps) {
						runningFirstHapIdx += 2;
					}
					else {
						runningFirstHapIdx += 1;
					}
										
					// now we know how many we want
					// continue if we don't want (them)
					if (!haplotypesToRead.get(firstHapIdx)) {
						if (twoHaps) {
							if (!haplotypesToRead.get(firstHapIdx+1)) {
								continue;
							}
						}
						else {
							continue;
						}
					}

					
					byte[] thisAlleles = null;
					// everything missing here?
					if (vcfEntryToAlleleIdx == null) {
						if (twoHaps) {
							thisAlleles = new byte[] {ReadSequences.MISSING_ALLELE, ReadSequences.MISSING_ALLELE};
						}
						else {
							thisAlleles = new byte[] {ReadSequences.MISSING_ALLELE};
						}
					}
					else {
						// and see what you get as the actual alleles
						thisAlleles = getAllelesFromVcfEntry (thisSNP, vcfEntryToAlleleIdx, indMissingCharSet, acceptUnphasedAsMissing, unphasedCount);
					}
					
					assert (thisAlleles.length >= 1);
					assert (thisAlleles.length <= 2);
					
					// and then store them in the appropriate places
					// only the ones we want
					if (haplotypesToRead.get(firstHapIdx)) currVariants.add (thisAlleles[0]);
					if (twoHaps && haplotypesToRead.get(firstHapIdx+1)) currVariants.add (thisAlleles[1]);
				}
				
				// check that we have right number of variants
				if (runningFirstHapIdx !=  haplotypesToRead.size()) throw new IOException ("Number of haplotypes in line " + lineNumber + " differs between vcf file and config file (" + runningFirstHapIdx + ", " + haplotypesToRead.size() + ").");
				if (currVariants.size() != numHapsToRead) throw new IOException ("Number of haplotypes read in line " + lineNumber + " inconsistent with config file (" + currVariants.size() + ", " + haplotypesToRead.size() + ").");


				// only keep the site if it is really segregating
			
				// default no
				boolean segregating = false;
				// see about whether this site is actually segregating
				// therefore collect all the alleles in a set
				TByteHashSet set = new TByteHashSet();

				assert (currVariants.size() == numHapsToRead);
				for (int j=0; j<currVariants.size(); j++) {
					set.add(currVariants.get(j));
				}
				if (set.size() > 1) {
					// it is actually segregating
					segregating = true;
				}
				
				// if it is not segregating, we don't keep it
				if (!segregating) {
					// set the allele in the reference
					byte[] bytes = set.toArray();
					assert (bytes.length == 1);
					referenceSequence[pos] = bytes[0];
					
					// keep a count
					nonSegSites++;
				}
				else {
					// if it really is segregating, we have to mark this in the reference
					referenceSequence[pos] = (byte) ReadSequences.SEG_SITE_IN_REFERENCE;
					
					// add the current variants
					variants.add (currVariants);

					// also add the segregating site now
					segSites.add (pos);

					// but see if one at least one of them is missing
					// (we do have other alleles, if not, we wouldn't be here)
					if (set.contains ((byte) ReadSequences.MISSING_ALLELE)) {
						numWithMissing++;
					}
				}
			}
		}
		if (!expectSNPs) {
			if (useExternalReference) {
				throw new IOException("No header given in the vcf file.");
			}
			else {
				throw new IOException("No reference sequence or header given in the vcf file. You can provide an external reference file with --vcfReferenceFile.");
			}
		}
		
		// just make sure
		int numSites = segSites.size();
		assert (numVariationLines == numSites + tooManyAlleles + filterNotPassed + nonSegSites + numDoubleVcfEntries);
		
		System.out.println("# VCF: Read " + numVariationLines + " lines with genetic variation from the vcf-file.");
		if (vcfIgnoreDoubleEntries) {
			System.out.println ("# VCF: " + numDoubleVcfEntries + " (out of " + numVariationLines + ") lines in vcf-file were double entries (same position) and ignored.");
		}
		System.out.println ("# VCF: " + nonSegSites + " (out of " + numVariationLines + ") sites were not segregating (for the chosen haplotypes).");
		System.out.println ("# VCF: " + tooManyAlleles + " (out of " + numVariationLines + ") sites had more than two alleles, and were ignored.");
		
		System.out.println ("# VCF: " + filterNotPassed + " (out of " + numVariationLines + ") sites did not pass the filter.");
		
		// some checking that everything looks good
		// check for ordering and duplicate entries, mayhaps do this checking earlier
		Integer prev = -1;
		for (Integer curr : segSites) {
			// vcf convention
			if (curr < prev) throw new IOException ("VCF file contains wrongly ordered sites: (" + (prev+1) + ", " + (curr+1) + ").");
			if (curr == prev) throw new IOException ("VCF file contains duplicate entry for site " + (curr+1) + ".");
			// very important to add this
			prev = curr;
		}
		
		// having no segregating sites is actually totally fine
		if (segSites.size() > 0) {
			// check first pos
			if (segSites.get(0) < 0) throw new IOException ("First site in VCF is smaller than 1.");
			// check last pos
			if (segSites.get(segSites.size()-1) >= referenceSequence.length) throw new IOException ("Last site in VCF is too big for given reference.");
		}

		
		//tell everybody how many unphased sites you threw away
		System.out.println("# VCF: " + unphasedCount.getValue() + " unphased genotypes were converted into missing alleles.");

		
		// then get the actual compressed haplotypes
		// somehow we need to get the actual number of segregating sites
		numSites = 0;
		for (Integer segSite : segSites) {
			if (segSite != null) numSites++;
		}
		System.out.println("# VCF: " + numSites + " segregating site(s) left after pre-processing and filtering.");
		
		// if no segregating sites left, complain
		if (numSites < 1) {
			throw new IOException ("No segregating sites left in the input vcf after pre-processing and filtering. Perhaps you did not specify the correct --vcfFilterPassString.");
		}
		
		// what about partially missing
		System.out.println("# VCF: At " + numWithMissing + " site(s) (out of " + numSites + " segregating site(s)) at least one haplotype has a missing allele and will potentially be ignored if loci are grouped together (in mutli-locus step handler).");
		
		// and now get the actual haplotypes out into an array
		// 2*numIndividualsToRead, because now it gets haploid
		byte[][] rawHaplotypes = new byte[2*numHapsToRead][numSites];
		TIntIntMap referenceIdxMap = new TIntIntHashMap();
		int currCompressedSite = 0;
		for (int i=0; i<variants.size(); i++) {
			// do nothing if it is null
			TByteList currVariants = variants.get(i);
			if (currVariants != null) {
				// just some checking
				assert (referenceSequence[segSites.get(i)] == ReadSequences.SEG_SITE_IN_REFERENCE);
				
				// also remember it in map
				referenceIdxMap.put (segSites.get(i), currCompressedSite);
				
				// copy it
				for (int j=0; j<currVariants.size(); j++) {
					rawHaplotypes[j][currCompressedSite] = currVariants.get(j);
				}
				// and increase the current compressed site
				currCompressedSite++;
			}
		}
		// we should have kind of arrived at the same value
		assert (currCompressedSite == numSites);
		

		// now we finally can create the haplotypes
		List<HFSAXFullHaplotype> hapList = new ArrayList<HFSAXFullHaplotype>();
		// modify reference if wanted
		if (numInternalAlleles < 4) {
			for (int i=0; i<referenceSequence.length; i++) {
				if ((referenceSequence[i] != ReadSequences.MISSING_ALLELE) && (referenceSequence[i] != ReadSequences.SEG_SITE_IN_REFERENCE)) {
					referenceSequence[i] = 0;
				}
			}
		}
		// one reference to rule them all
		FSAReference ref = new SimpleFSARef (referenceSequence, referenceIdxMap); 
		// iterate through the haplotypes
//		assert (numCol <= haplotypesToRead.size());
//		assert (haplotypesToRead.size() <= numCol*2);

		// should work
		int rawIdx = 0;
		for (int i=0; i<haplotypesToRead.size(); i++) {
			HFSAXFullHaplotype hap = null;
			// did we read something here?
			if (haplotypesToRead.get(i)) {
				// yes, so make haplotype with reference
				hap = new HFSAXFullHaplotype(ref, rawHaplotypes[rawIdx]);
				rawIdx++;
			} // if not, null is ok
			// add either the things we read or null for we did not read anything
			hapList.add (hap);
		}
		// seems to be fine
		
		// and give them away
		return hapList;
	}

	private static byte[] handleReferenceFile (boolean fromVCF, String referenceFileName, List<Interval> bedRegions, char[] refAlleles, char[] commentCharacters, char[] realMissingChars, boolean verbose) throws IOException {
		// some initial reference 
		byte[] referenceSequence = null;
		
		boolean referenceExists = (new File(referenceFileName)).exists();
		if (verbose) {
			if (fromVCF) {
				System.out.println("# VCF: reference file: " + referenceFileName + ". Exists? " + (new File(referenceFileName)).exists());
			}
			else {
				System.out.println("# Reference file: " + referenceFileName + ". Exists? " + (new File(referenceFileName)).exists());
			}
		}
		if (!referenceExists) {
			if (fromVCF) {
				throw new IOException ("Error in VCF: Cannot access reference file specified in VCF-file: " + referenceFileName + ". Try providing it via --referenceFile.");
			}
			else {
				throw new IOException ("Error reading reference: specified reference file does not exist: " + referenceFileName);
			}
		}
		
		
		// read reference
		BufferedReader refReader = new BufferedReader(new FileReader(referenceFileName));
		Integer refNumAlleles = refAlleles.length;
		// reference might have special alleles
		// read it fully, even if we modify it later
		referenceSequence = readVcfReferenceSequence (refReader, bedRegions, refNumAlleles, refAlleles, commentCharacters, realMissingChars); 
		
		if (verbose) {
			System.out.println("# Length of reference sequence is " + referenceSequence.length + " bp.");
		}

		// that should be it
		return referenceSequence;
	}

	private static byte[] readVcfReferenceSequence(BufferedReader refReader, List<Interval> bedRegions, int numAlleles, char[] alleles, char[] commentCharacters, char[] missingCharacters) throws IOException {
		// read whatever we got
		assert ((numAlleles > 1) && (numAlleles <= 4));
		List<byte[]> stuff = readRawSequences(refReader, null, numAlleles, alleles, commentCharacters, missingCharacters, true);
		byte[] refSeq = null;
		if (stuff.size() != 1) {
			// first get total size
			int refLen = 0;
			for (int i=0; i<stuff.size(); i++) {
				refLen += stuff.get(i).length;
			}
			
			// then make a new long sequence and copy everything
			refSeq = new byte[refLen];
			int lastDstPos = 0;
			for (int i=0; i<stuff.size(); i++) {
				byte[] srcArray = stuff.get(i);
				System.arraycopy (srcArray, 0, refSeq, lastDstPos, srcArray.length);
				lastDstPos += srcArray.length;
			}
		}
		else {
			refSeq = stuff.get(0);
		}
		
		// apply the bed
		if(bedRegions != null){
			for(Interval currInterval : bedRegions){
				for(int i = (int)currInterval.startPoint; i < currInterval.endPoint; i++){
					refSeq[i] = -1;	//this one lives in an excluded region, so set him to missing!
				}
			}
		}
		
		return refSeq;
	}

	private static Map<Character, Integer> getAlleleToIdxMap (char[] alleles, int numAlleles) {
		Map<Character, Integer> returnMap = new TreeMap<Character, Integer>();
		// fill it
		if (alleles == null) {
			assert (numAlleles < 10);
			// just map the integers to a number
			for (int i=0; i<numAlleles; i++) {
				returnMap.put (String.valueOf(i).charAt(0), i);
			}
		}
		else {
			// we actually have some alleles
			assert (numAlleles == alleles.length);
			for (int i=0; i<numAlleles; i++) {
				returnMap.put(alleles[i], i);
			}
		}
		// should be good
		return returnMap;
	}

	private static byte[] getAllelesFromVcfEntry (String vcfEntry, int[] vcfEntryToAlleleIdx, TreeSet<Character> missingCharSet, boolean acceptUnphasedAsMissing, MutableInt unphasedCount) throws IOException {
		
		// split it
		char[] letters = vcfEntry.toCharArray();
		
		if (letters.length == 1) {
			// haploid, no question about phasing
			return new byte[] { (byte) getAllele (letters[0], vcfEntryToAlleleIdx, missingCharSet) };
		}
		else if (letters.length == 3) {
			// diploid
			byte[] theAlleles = new byte[] { (byte) getAllele (letters[0], vcfEntryToAlleleIdx, missingCharSet), (byte) getAllele (letters[2], vcfEntryToAlleleIdx, missingCharSet) };
			
			// there might be questions about phasing
			// if the two alleles are equal, there is no question about phasing
			if (theAlleles[0] != theAlleles[1]) {
				// no we look at whether it's phased
				if (letters[1] != '|') {
					if (acceptUnphasedAsMissing) {
						//for now we just mark this site as missing in both haps!
						unphasedCount.increase();
						theAlleles[0] = ReadSequences.MISSING_ALLELE;
						theAlleles[1] = ReadSequences.MISSING_ALLELE;
					}
					else {
						// raise hell
						throw new IOException ("Invalid entry in vcf file: \"" + vcfEntry + "\" . Genotype is not phased. If you want to accept unphased genotypes as missing alleles, use --acceptUnphasedAsMissing.");
					}
				}
			}
			return theAlleles;
		}
		else {
			throw new IOException ("Invalid entry in vcf file: \"" + vcfEntry + "\" . Genotype has to consist of either one allele (haploid) or two alleles and a divider.");
		}
	}

	private static int getAllele (char c, int[] vcfEntryToAlleleIdx, TreeSet<Character> missingCharSet) throws IOException {
		if (missingCharSet.contains(c)) {
			return ReadSequences.MISSING_ALLELE;
		}
		else {
			// should be a valid allele
			int index = c - '0';
			if ((index < 0) || (index >= vcfEntryToAlleleIdx.length)) throw new IOException("Invalid allele in VCF file.");
			int daAllele = vcfEntryToAlleleIdx[index];
			return daAllele;
		}
	}
	
}
