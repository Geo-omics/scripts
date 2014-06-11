## Script Descriptions
* **addFileName2header.pl**	-	Add the name of the fasta/fastq file (without the extension) to the header. Useful when merging multiple fasta/fastq files and wish to easily keep track of sequences
* **batchBlast.pl**	-	Run multiple Blasts at once in an "embarassingly" parallel manner
* **calcN50.pl**	-	Calculate N50 and L50 values for a fasta file
* **chopper.pl**	-	Chop a file (fasta/fastq/tab-delimited/multiple-line) into multiple parts.
* **createFastq.pl**	-	Use the fasta and quality files to produce fastq files
* **CRISPR\_spacer\_extractor.pl**	-	Given a fasta file with repeat sequences and a contig fasta file, get the positions of these repeats in the contigs and find the coordinates of spacers
* **dereplicate.pl**	-	Remove sequences that are exactly the same and maintain a record of the clusters.
* **extractSeqs.pl**	-	Given a list of sequence names, extract the sequences from a fasta/fastq file
* **findStretchesOfNs.pl**	-	Go through the fasta file and find the sequences with a 100 or more N's at a stretch.
* **gcSkew.pl**	-	Calculate the GC skew for for a fasta file
* **getRandomData.pl**	-	Extract a random % of data (fasta/fastq)
* **interleave.pl**	-	Take the forward(1) and reverse(2) fasta/fastq files and arrange them such that all odd sequences are forward and evens are reverse.
* **kmerFreq.pl**		-	Calculate any kmer frequency from the given fasta
* **length+GC.pl**	-	Calculate the length and  GC content from the fasta file.
* **limit2Length.pl**	-	Remove sequences from a fasta file that do not pass the length threshold set by the user. Print length distribution to the screen
* **separateInterleaved.pl**	-       Separate interleaved files.
* **usageStats.pl**	-	Monitor a given process and email report when process finishes
