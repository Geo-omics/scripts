#!/user/bin/perl -w
=head1 DESCRIPTION

	Dereplicates a Fasta file at 100% identity over 100% coverage. Picks a first sequence as the representative for the cluster.
	Dereplicates a Fastq file at 100% identity over 100% coverage. Picks the sequence with the best 'avg quality score' as the representative for the cluster.
	

=head2 USAGE
	
	perl dereplicate.pl -f fasta_File
	OR
	perl dereplicate.pl -fq fastq_File

=head3 NOTE

	## NOTE1: The Script DOES NOT look for sub-string matches. 
	## NOTE2: Default Phred offset is 33

=head2 OUTPUT

	The script will create the following as its output, the ".clust" and the ".clust.list". These files are explained below.
	The ".clust" file:
	This file contains the cluster number, cluster size and the names of all sequences in that particular cluster. The third column in the file is the representative for each cluster. The file is tab-delimited.

	The ".clust.list" file:
	This file contains a list of all representative sequences in the fastq file. The file is tab-delimited.
	
=head2 AUTHOR

	Sunit Jain, Oct, 2011

=head2 Contributors

	Chris Taylor, Nov, 2011

=cut


#######################
## MODULES
#######################
use strict;
use Getopt::Long;
use List::Util 'sum';
#use Digest::MD5 'md5';
use File::Basename;

#######################
## PARAMETERS
#######################
my $fasta;
my $fastq;
my $setScore=0;
my $phredOffset=33;
my $version="0.4.5";
GetOptions(
	'f:s'=>\$fasta,
	'fq:s'=>\$fastq,
	's:i'=>\$setScore,
	'p|phred_offset:i' => \$phredOffset,
	'h'=>sub{system('perldoc', $0); exit;},
);

#######################
## CHECKS
#######################
my $file;
if (!$fasta && !$fastq){system('perldoc', $0); exit;}
elsif($fasta && !$fastq){ $file=$fasta; }
elsif(!$fasta && $fastq){ $file=$fastq; }
elsif($fasta && $fastq){print "Choose either fasta or fastq file at a time\n";}

#######################
## GLOBAL
#######################
my ($f,$d)=fileparse($file);
my (@fileName)=split(/\./, $f);
my $ext=pop @fileName;
my $name=join(".", @fileName);
my $clust = $d.$name.".clust";
my $list = $d.$name.".clust.list";
my $fastaOut= $d.$name.".fasta";

$setScore=$setScore+$phredOffset;
$/=$fasta ? "\>" : "\n";

my %seen;
my %numbers;
my %qual;

#######################
## MAIN
#######################
open(FILE, $file);
while(my $line=<FILE>){
	$line=&trim($line);
	next unless $line;

	$fasta ? &parseFasta($line) : &parseFastq($line);
}
close FILE;

$fasta ? &fastaClustering : &fastqClustering;

# print top 10 duplicated sequences
my @sorted=sort{$numbers{$b} <=> $numbers{$a}} keys %numbers;
my $counter=0;
foreach my $s(@sorted){
	$counter++;
	print $s."\t".$numbers{$s}."\n";
	last if $counter==10;
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
		
		$line=<FILE>; # Sequence
		$line=&trim($line);
#		my $seqMD5=md5($line);
#		push(@{$seen{$seqMD5}}, $seqDesc);
		my $seq=$line;
		push(@{$seen{$seq}}, $seqDesc);

		$line=<FILE>; # Quality Header

		$line=<FILE>; # Quality
		$line=&trim($line);
		&seqScore($seqDesc, $line);
		return $seqDesc;
	}
	else{ die "ERROR: Script Borked at line $.! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n"; }
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
	$qual{$seqDesc}=$totalScore;
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
		
		print CLUST "c".$clustNum."\t".$size."\t";
		my $bestSeq="";
		my $bestSeqQual=-100;
		# Pick the Seq with the Best Score
		foreach my $seqDesc(@{$seen{$seq}}){
			my $seqQual=$qual{$seqDesc};		
			if ($seqQual > $bestSeqQual){
				$bestSeq=$seqDesc;
				$bestSeqQual=$seqQual;
			}
		}
		print CLUST $bestSeq."\t";
		print LIST $bestSeq."\n";
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
#	my $seqMD5=md5($seq);
#	push(@{$seen{$seqMD5}}, $seqDesc);
	push(@{$seen{$seq}}, $seqDesc);

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


