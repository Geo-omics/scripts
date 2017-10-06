#! usr/bin/perl

=head1 DESCRIPTION

	chopper.pl - chops a Fasta/Fastq file into a defined number of portions.

=head1 USAGE

	perl chopper.pl [-f fasta_Filename.fasta]
	OR
	perl chopper.pl [-fq fastq_Filename.fastq]

=head2 Optional

	-p OR -parts:	integer	Split the given fasta file into these many parts.
	-s OR seqs:	integer	Maximum number of Sequences allowed in a file.
	-t OR totalseqs:	integer	If you already know the total number of sequences in the file, mention here, this will greatly speed up the script.
	-avgsize:	integer	defines the average size of a fasta file. Only used when you don't specify '-p' OR '-s'.
	-maxbase:	integer	chop the fasta file upto a certain sequence length
	-overlap:	integer	chop a fasta sequence such that all sequences are AT MOST "maxbase" long and if a sequence has been chopped there is an overlap of 'overlap' bases (default = 0) between the two sequences.

=head1 AUTHOR

	Sunit Jain, Aug, 2011.
	sunitj [AT] umich [DOT] edu.
	 
=cut

use strict;
use Getopt::Long;
use File::Spec;

my $fasta;
my $fastq;
my $parts;
my $totalNumSeqs;
my $numSeqs; # max number of seqs allowed in a file
my $avgSize= 1000000;
my $maxBases;
my $overlap;
my $version="0.1.7";

GetOptions(
	'f|fasta:s'=>\$fasta,
	'fq|fastq:s'=>\$fastq,
	'p|parts:i'=>\$parts,
	't|totalseqs:i'=>\$totalNumSeqs,
	's|seqs:i'=>\$numSeqs,
	'avgsize:i'=>\$avgSize,
	'maxbase:i'=>\$maxBases,
	'overlap:i'=>\$overlap,
	'h|help'=>sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print $version."\n"; exit;},
);

my ($totalSeqs, $filePartSize);
&checkArgs;
#print "Seq:$numSeqs\n";
#print "Parts:$parts\n";
my $lastFileName="";
my $file = $fasta ? $fasta : $fastq;
&chopFile($file);
&checkLastFile;

sub checkArgs{
	if (! $fasta && ! $fastq){
		system('perldoc', $0); exit;
	}

	if ($overlap && ! $maxBases){
		die "ERROR: '-maxbase' parameter required with '-overlap'\n";
	}
	return if ($maxBases);

	if (! $totalNumSeqs){
		$totalSeqs= $fasta ? `grep -c \'>\' $fasta` : `wc -l $fastq`;
		$totalSeqs=$totalSeqs/4 if ($fastq);
	}
	else{
		$totalSeqs=$totalNumSeqs;
	}
	print "Total Sequences:\t$totalSeqs\n";

	$filePartSize=$totalSeqs/$avgSize;
	if ((! $parts || $parts == 0 ) && (! $numSeqs || $numSeqs == 0)){
		($filePartSize >= 1) ? &bigFile : &smallFile;
	}
	elsif ($parts > 0 && $numSeqs > 0){
		warn "[WARNING] You may only use one of the two options (Number of Parts OR Max Seq) at a time.\n";
		warn "[WARNING] Only using the 'Number of Parts' cut-off\n";
		&calcNumSeqs;
	}
	elsif ($parts > 0 && (! $numSeqs || $numSeqs == 0)){
		&calcNumSeqs;		
	}
	elsif ((! $parts || $parts == 0) && $numSeqs > 0){
		return;
	}
	return;
}

# Get number of seqs to chop the file into; if the user doesn't mention it.
sub bigFile{
	exit if (! $totalSeqs);
	$parts= int($filePartSize + 0.5);
	warn "[Warning]This file seems like it's larger than my definition of an average file ($avgSize sequences).\n";
	print "Since you didn't mention a preference I'll split the Fasta file into \~ $parts parts\n";
	print "You may change the 'average file size' by using the '-avgsize' flag.\n";
	$numSeqs= int(($totalSeqs/$parts) + 0.5);
}

sub smallFile{
	exit if (! $totalSeqs);
	my $p= int($filePartSize + 0.5);
	$parts= ($p == 0) ? 1 : $p;
	warn "[Warning]This file seems like it's smaller than my definition of an 'average file($avgSize sequences).\n";
	print "Since you didn't mention a preference I'll split the Fasta file into \~ $parts parts\n";
	print "You may change the 'average file size' by using the '-avgsize' flag.\n";
	$numSeqs= int(($totalSeqs/$parts) + 0.5);
}

sub calcNumSeqs{
	exit if (! $totalSeqs);
	$numSeqs= int(($totalSeqs/$parts) + 0.5);
	print "Each file will have about $numSeqs sequences\n";
}

sub chopFile{
	my $file=shift;

	open (FILE, $file) || die "[error] $! : $file\n";
	my $numOfSeqs=0;
	my $partCount=1;
	my $totalSeqs=0;
	my ($volume,$directories,$fileName)=File::Spec->splitpath( $file );
	my ($fName, @ext)= split(/\.([^\.]+)$/, $fileName);
	my $restOfName=join(".", @ext);
	my $out=$fName.".0".$partCount.".".$restOfName;
	my $fh;
	open ($fh, ">".$out) || die "WTF!! : $!";
	$/= $fasta ? "\>" : "\n";
	my $delim= $fasta ? "\>" : "\@";
	my %seen=();

	while (my $b = <FILE>) {
		$b= &trim($b);
		next unless $b;
		my ($name, @sequence, $seq, $qName, $qual);
		if ($fasta){
			($name, @sequence)= split (/\n/, $b);
			$seq=join("", @sequence);
			if (length($seq)==0){
				my $line=<FILE>;
				(my $restName, @sequence)=split(/\n/, $line);
				$seq=join("", @sequence);
				$name=$name."^".$restName;
			}
			$b=$name."\n".$seq."\n";
		}
		elsif($fastq){
			if ($b=~/^@(\S+)/){
				$b=~ s/^@//;
				$name=$b;
				$qName="+".$b;

				$b=<FILE>;
				$b= &trim($b);
				$seq=$b;

				$b=<FILE>;

				$b=<FILE>;
				$b= &trim($b);
				$qual=$b;

				$b=$name."\n".$seq."\n".$qName."\n".$qual."\n";
			}
			else{ 
				die "ERROR: Script Borked! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n";
			}
		}

		if($overlap && $fasta){
			&overlapChop($name, $seq, $fh);
			next;
		}		
		exit if($overlap && $fasta);

		if ($maxBases && $fasta){
			&chopSeq($name, $seq, $fh);
			next;
		}
		exit if ($maxBases && $fasta);

		$numOfSeqs++;
		$totalSeqs++;
		if ($numOfSeqs < $numSeqs){
			print $fh $delim.$b;
		}
		elsif ($numOfSeqs == $numSeqs){
			print $fh $delim.$b;
			close $fh;

			$numOfSeqs=0;
			$partCount++;
			my $out=$fName.".0".$partCount.".".$restOfName;
			$lastFileName= $out;
			open ($fh, ">".$out);
		}
	}
}

sub trim{
	my $item=shift;
	chomp($item);
	$item=~ s/\r//g;
	$item=~ s/^\s+//;
	$item=~ s/\s+$//;
	return $item;
}

sub chopSeq{
	my ($name, $seq, $fh)=@_;
	my $seqLen=length($seq);
	if ($seqLen > $maxBases){
		my $p= $seqLen/$maxBases;
		my $diff= $p - int($p);
		my $seqParts;
		$diff==0 ? $seqParts=$p : $seqParts= $p + 1;
		my $pos=0;
		for (my $i=1; $i<$seqParts; $i++){
			print $fh ">".$name."_Part_".$i."/".int($seqParts)."\n";

			my $partSeq=substr($seq, $pos, $maxBases);
			$pos=$pos+$maxBases;
			print $fh $partSeq."\n";
		}
	}
	else{
		print $fh ">".$name."_Part_1/1\n";
		print $fh $seq."\n";
	}
	return;
}

sub overlapChop{
	my ($header, $seq, $fh)=@_;
	my $seq_len= length $seq;
	my $counter=0;
	if ($seq_len <= $maxBases) {
		print $fh ">".$header."_Part_".$counter."\n".$seq."\n";
		return $counter;
	}

	for (my $i = 0; $i < $seq_len; $i += $maxBases - $overlap) {
		my $contig = substr $seq, $i, $maxBases;
		if ((length $contig) == $overlap) { # for this case, don't output tc on second line acctcacctcc 5
			return $counter;
		}
		#print "$contig $seq $i $maxBases $seq_len \n";
		$counter++;
		print $fh ">".$header."_Part_".$counter."\n".$contig."\n";
		if ((length $contig) < $maxBases ) { # handle boundary conditions, two case acctcc 5
			# reached the end
			return $counter;
		}
	}
}

sub checkLastFile{
	unlink $lastFileName if (-z $lastFileName);
}
