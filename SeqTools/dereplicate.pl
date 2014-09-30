#!/user/bin/perl -w
=head1 DESCRIPTION

	Dereplicates a Fasta file at 100% identity over 100% coverage. Picks a first sequence as the representative for the cluster.
	Dereplicates a Fastq file at 100% identity over 100% coverage. Picks the sequence with the best 'avg quality score' as the representative for the cluster.
	

=head2 USAGE
	
	perl dereplicate.pl -f fasta_File -out output.fasta
	OR
	perl dereplicate.pl -fq fastq_File -out output.fastq 
	OR
	perl dereplicate.pl -fq fastq_File -out output.fasta -outfmt fasta
	
=head2 Options

	-f	or	-fasta	[characters]	fasta file input
	-fq	or	-fastq	[characters]	fastq file input
	-o	or	-out	[characters]	dereplicated fasta/fastq file. Note that an output file will only be created if this option is provided.
	-outfmt			["fasta"]	Use this if you wish that the output for your fastq input be fasta. See example #3 above.
	-n	or	-top	[integer]	Print top N sequence clusters. Default 10.

=head3 NOTE

	## NOTE1: The Script DOES NOT look for sub-string matches. 
	## NOTE2: Default Phred offset is 33
	## Note3: An output file will only be created if the '-o' option is provided.

=head2 OUTPUT

	The script will create the following as its output, the ".clust", the ".clust.list" and optionaly the fasta/fastq file. These files are explained below.
	The ".clust" file:
	This file contains the cluster number, cluster size and the names of all sequences in that particular cluster. The third column in the file is the representative for each cluster. The file is tab-delimited.

	The ".clust.list" file:
	This file contains a list of all representative sequences in the fastq file. The file is tab-delimited.
	
=head2 AUTHOR

	Sunit Jain, Oct, 2011
	sunitj [at] umich [dot] edu

=head2 Contributors

	Chris Taylor, Nov, 2011

=cut


#######################
## MODULES
#######################
use strict;
use Getopt::Long;
use List::Util 'sum';
use File::Basename;

#######################
## PARAMETERS
#######################
my $fasta;
my $fastq;
my $usearch;
my $setScore=0;
my $phredOffset=33;
my $out;
my $outfmt;
my $n=10;
my $version="dereplicate.pl v0.6.2";
GetOptions(
	'f:s'=>\$fasta,
	'fq:s'=>\$fastq,
	's:i'=>\$setScore,
	'o|out:s'=>\$out,
	'outfmt:s'=>\$outfmt,
	'usearch'=>\$usearch,
	'n|top:i'=>\$n,
	'p|phred_offset:i' => \$phredOffset,
	'h|help'=>sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print $version."\n"; exit;},
);

#######################
## CHECKS
#######################
my $file;
if (!$fasta && !$fastq){system('perldoc', $0); exit;}
elsif($fasta && !$fastq){ $file=$fasta; }
elsif(!$fasta && $fastq){ $file=$fastq; }
elsif($fasta && $fastq){die "Choose either fasta or fastq file at a time\n";}

my $seqType= $fastq ? "fastq" : "fasta";
if (! $outfmt){$outfmt = $seqType;}

#######################
## GLOBAL
#######################
my ($f,$d)=fileparse($file);
my (@fileName)=split(/\./, $f);
my $ext=pop @fileName;
my $name=join(".", @fileName);
my $clust = $d.$name.".clust";
my $list = $d.$name.".clust.list";

$setScore=$setScore+$phredOffset;
$/=$fasta ? "\>" : "\n";

my %seen;
my %numbers;
my %qual;
my $totalSequences;
my $derepSequences;
#######################
## MAIN
#######################
open(FILE, $file)|| die $!;
open(OUT, ">".$out)|| die $! if ($out);
while(my $line=<FILE>){
	$line=&trim($line);
	next unless $line;

	$fasta ? &parseFasta($line) : &parseFastq($line);
}
close FILE;

$fasta ? &fastaClustering : &fastqClustering;
close OUT if ($out);
# print top N duplicated sequences
print "# Original # of Sequences in $file:\t".$totalSequences."\n";
print "# Sequences after duplicate removal:\t".keys(%numbers)."\n";
if($n){
	my @sorted=sort{$numbers{$b} <=> $numbers{$a}} keys %numbers;
	my $counter=0;
	foreach my $s(@sorted){
		$counter++;
		print $s."\t".$numbers{$s}."\n";
		last if $counter==$n;
	}
}
exit 0;

#######################
## SUB-ROUTINES
#######################

sub trim{
	my $line=shift;
	chomp($line);
	$line=~ s/\r//;
	$line=~ s/^\s+//;
	$line=~ s/\s+$//;
	return $line;
}

## For Fastq ##

sub parseFastq{
	my 	$line=shift;

	if ($line=~ /^@\w+/){
		my $seqDesc=$line; # Sequence Header
		$seqDesc=~ s/^@//;
		
		$line=<FILE>; # Sequence
		$line=&trim($line);
		my $seq=$line;
		push(@{$seen{$seq}}, $seqDesc);

		$line=<FILE>; # Quality Header

		$line=<FILE>; # Quality
		$line=&trim($line);
		die "[ERROR: line $.] Script Borked! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n" if (length($seq)!= length($line));

		&seqScore($seqDesc, $line);
		$totalSequences++;
		return $seqDesc;
	}
	else{ die "[ERROR: line $.] Script Borked! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n"; }
}

sub seqScore{
	my ($seqDesc, $qual)=@_;

	my @scores=split(//,$qual);

	### SCORING #####################################
	@scores=map{ord($_)} @scores;
	#################################################

	my $totalScore=sum(@scores);
	# If score is lower than the threshold; get rid of the seq.
	return if ($totalScore < $setScore);
	$qual{$seqDesc}{"Total"}=$totalScore;
	$qual{$seqDesc}{"Qual"}=$qual;
	return;
}

sub fastqClustering{
	open(CLUST, ">".$clust);
	open(LIST, ">".$list);
	my $clustNum=0;
	print CLUST "#ClusterNumber\tSize\tRepresentative\tSeqHeaders\n";
	print LIST "#Representatives\n";
	foreach my $seq(keys %seen){
		$clustNum++;
		my $size=@{$seen{$seq}};
		$numbers{$seq}=$size;
		$derepSequences+=$size;
		
		print CLUST "c".$clustNum."\t".$size."\t";
		my $bestSeq="";
		my $bestSeqQual=-100;
		# Pick the Seq with the Best Score
		foreach my $seqDesc(@{$seen{$seq}}){
			my $seqQual=$qual{$seqDesc}{"Total"};		
			if ($seqQual > $bestSeqQual){
				$bestSeq=$seqDesc;
				$bestSeqQual=$seqQual;
			}
		}
		print CLUST $bestSeq."\t";
		print LIST $bestSeq."\n";
		
		# print dereplicated file.
		if($out){
			$outfmt eq "fasta" ? &printFasta($bestSeq, $seq, $size) : &printFastq($bestSeq, $seq);
		}
		
		foreach my $h(@{$seen{$seq}}){
			print CLUST $h."\t" unless ($h eq $bestSeq);
		}
		print CLUST "\n";

		delete $seen{$seq};
	}
	close CLUST;
	return;
}

## For Fasta ##

sub parseFasta{
	my 	$line=shift;
	my($seqDesc,@sequence)=split(/\n/, $line);
	my $seq = join ("", @sequence);
	$seq=~ s/\r//g;
	$seqDesc=~ s/^>//;

	push(@{$seen{$seq}}, $seqDesc);
	$totalSequences++;

	return;
}

sub fastaClustering{
	open(LIST, ">".$list);
	open(CLUST, ">".$clust);
	my $clustNum=0;
	print CLUST "#ClusterNumber\tSize\tRepresentative\tSeqHeaders\n";
	print LIST "#Representatives\n";
	foreach my $seq(keys %seen){
		$clustNum++;
		my $size=@{$seen{$seq}};
		$numbers{$seq}=$size;
		$derepSequences+=$size;

		# print dereplicated file.
		if ($out){
			&printFasta($seen{$seq}[0], $seq, $size);
		}
		
		print CLUST "c".$clustNum."\t".$size."\t";
		foreach my $h(@{$seen{$seq}}){
			print CLUST $h."\t";
		}
		print LIST $seen{$seq}[0]."\n";
		print CLUST "\n";

		delete $seen{$seq};
	}
	close LIST;
	close CLUST;
	return;
}

## Print to files ##
sub printFasta{
	my ($header, $seq, $size)=@_;
	if($usearch){
		$header=~ s/ /\_/g;
		print OUT ">".$header.";size=".$size.";\n";
	}
	else{
		print OUT ">".$header."\n";
	}
	print OUT $seq."\n";
}

sub printFastq{
	my($header, $seq)=@_;
	print OUT "@".$header."\n";
	print OUT $seq."\n";
	print OUT "\+\n";
	print OUT $qual{$header}{"Qual"}."\n";
}
