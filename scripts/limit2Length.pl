#!/user/bin/perl -w
=head1 DESCRIPTION

	Check length of each sequence and only print those that clear the minimum threshold of x (default=60) bases;

=head1 USAGE
	
	perl limit2Length.pl -f fasta_File -len # [-o output file]
	OR
	perl limit2Length.pl -fq fastq_File -len # [-o output file]
	OR
	perl limit2Length.pl -below -len 2001 -f fasta_file -len # [-o output file]

=head2 OPTIONS

    -len    [INTEGER] : Get sequences of (length+1).
            Default =  Sequences greater than this parameter will be included. To change this behavior see '-below'.
	-b  [INTEGER]:	bin size for length distribution; default=2000
	-below [BOOLEAN]:	get sequences equal to or less than the 'len' threshold. 
	        Default = Only greater than -len.
	-scripts    [CHARACTER]: change the default directory for script dependencies: default: /geomicro/data1/COMMON/scripts/
	
=head3 Boolean Flags

	-quick:	don't get the distribution.
	-n50: also calculate the N50 and L50 stats for a *fasta* file

=head1 AUTHOR

	Sunit Jain, Mar, 2012
	sunitj-at-umich-dot-edu

	last updated: June 2013

=cut


#######################
## MODULES
#######################
use strict;
use Getopt::Long;
use File::Spec;
#######################
## PARAMETERS
#######################
my $version="0.1.8";
my $fasta;
my $fastq;
my $setLen=2000;
my $out=$$.".fa";
my $quick;
my $binSize=2000;
my $lower;
my $getN50;
my $scripts="/geomicro/data1/COMMON/scripts/SeqTools";
GetOptions(
	'f:s'=>\$fasta,
	'fq:s'=>\$fastq,
	'len:i'=>\$setLen,
	'o|out:s'=>\$out,
	'b|bin:i'=>\$binSize,
	'quick'=>\$quick,
	'below'=>\$lower,
	'n50'=>\$getN50,
	'scripts:s'=>\$scripts,
	'h|help'=>sub{system("perldoc", $0); exit;},
	'v|version'=>sub{print STDERR "$0 version: ".$version."\n"; exit;}
);

#######################
## CHECKS
#######################
my $file;
if (!$fasta && !$fastq){system('perldoc', $0); exit;}
elsif($fasta && !$fastq){ $file=$fasta; }
elsif(!$fasta && $fastq){ $file=$fastq; }
elsif($fasta && $fastq){print "Choose either fasta or fastq file at a time\n";}

die "[ERROR: $0] Minimum length threshold MUST be a positive integer greater than 1.\n" if ($setLen <= 1);

#######################
## GLOBAL
#######################

if(! $lower){	$setLen--; }

open (OUT, ">".$out)|| die "[ERROR: $0]: $!";

$/=$fasta ? "\>" : "\n";
my $count=0;
my $total=0;
my $sumLen=0;
my $revisedSumLen=0;
my (@sequence, %bins, $shortest, $longest);
my ($A, $T, $G, $C, $N)=(0,0,0,0,0);
#######################
## MAIN
#######################
unless($quick){
	@sequence= map { $binSize * $_ } 1 .. 50; # number of intervals
	@sequence=sort{$a <=> $b} @sequence;
	push(@sequence, (($sequence[-1])*10)); # capture REALLY long contigs.

	$shortest=100000000000000000000000; # ridiculously high number to capture the smaller ones
	$longest=-1;
}

open(FILE, $file);
while(my $line=<FILE>){
	$line=&trim($line);
	next unless $line;

	$count+= $fasta ? &parseFasta($line) : &parseFastq($line);
	$total++;
}
close FILE;
close OUT;

my $meanLenBefore= int(($sumLen/$total) + 0.5);
my $meanLenAfter= int(($revisedSumLen/$count) + 0.5);

print "# $0 version: $version\n";
print "# Total Sequences:\t$total\n";
print "# Total Bases:\t$sumLen\n";
$lower ? print "# Total Bases in sequences LESS THAN to $setLen:\t$revisedSumLen\n" : print "# Total Bases in sequences GREATER THAN $setLen:\t$revisedSumLen\n";

print "# Longest Sequence Length:\t$longest\n" unless($quick);
print "# Shortest Sequence Length:\t$shortest\n" unless($quick);

print "# Mean Length before setting the $setLen base limit:\t$meanLenBefore\n";
print "# Mean Length after setting the $setLen base limit:\t$meanLenAfter\n";

$lower ? print  "# Number of Sequences with Length LESS THAN $setLen bases:\t$count\n" : print "# Number of Sequences with Length GREATER THAN $setLen bases:\t$count\n";


my $perc_A = sprintf( "%.4f",(($A/$revisedSumLen)* 100));
my $perc_T = sprintf( "%.4f",(($T/$revisedSumLen)* 100));
my $perc_G = sprintf( "%.4f",(($G/$revisedSumLen)* 100));
my $perc_C = sprintf( "%.4f",(($C/$revisedSumLen)* 100));
my $perc_N = sprintf( "%.4f",(($N/$revisedSumLen)* 100));

print "# Nucleotide Distribution after setting the $setLen base limit:\n";
print "A\t:\t$perc_A %\n";
print "T\t:\t$perc_T %\n";
print "G\t:\t$perc_G %\n";
print "C\t:\t$perc_C %\n";
print "N\t:\t$perc_N %\n";

unless($quick){
	my $i=0;
	print "\nSequence Length Distribution:\n";
	foreach my $b(@sequence){
		my $from= $i > 0 ? ($sequence[$i-1])+1 : 1;
		print $from."-".$b."\t".$bins{$b}."\n" if $bins{$b};
		$i++;
	}
}

if ($getN50 && $fasta){
	my $calcN50=File::Spec->catfile( $scripts, "calcN50.pl");
	if (-e $calcN50){
		print "# N50 stats before setting the $setLen limit:\n";
		system("perl $calcN50 $fasta");
		print "\n";
	
		print "# N50 stats after setting the $setLen limit:\n";
		system("perl $calcN50 $out");
		print "\n";
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
		
		$line=<FILE>; # Sequence
		my $seq=&trim($line);
		my $len=length($seq);
		&binIt($len) unless ($quick);
		
		$line=<FILE>; # Quality Header

		$line=<FILE>; # Quality
		my $qual=&trim($line);

		$sumLen+=$len;

		if(($lower) && ($len < $setLen)){
			print OUT "$seqDesc\n$seq\n\+\n$qual\n";
			$revisedSumLen+=$len;
			&nucTally($seq);
			return 1;
		}
		elsif((! $lower) &&  ($len > $setLen)){
			print OUT "$seqDesc\n$seq\n\+\n$qual\n";
			&nucTally($seq);
			$revisedSumLen+=$len;
			return 1;
		}
		else{
			return 0;
		}
	}
	else{ die "ERROR: Script Borked! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n"; }
}

## For Fasta ##

sub parseFasta{
	my 	$line=shift;
	my($seqDesc,@sequence)=split(/\n/, $line);
	my $seq = join ("", @sequence);
	$seq=~ s/\r//g;

	my $len=length($seq);
	&binIt($len) unless ($quick);
	$sumLen+=$len;

	if(($lower) && ($len <= $setLen)){
		print OUT ">".$line."\n";
		$revisedSumLen+=$len;
		&nucTally($seq);	
		return 1;
	}
	elsif((! $lower) &&  ($len > $setLen)){
		print OUT ">".$line."\n";
		&nucTally($seq);
		$revisedSumLen+=$len;
		return 1;
	}
	else{
		return 0;
	}
}

sub nucTally {
	my $seq=shift;

	while ( $seq =~ /A/ig ) { $A++ }
    while ( $seq =~ /T/ig ) { $T++ }
    while ( $seq =~ /G/ig ) { $G++ }
    while ( $seq =~ /C/ig ) { $C++ }
    while ( $seq =~ /N/ig ) { $N++ }

	return;
}

sub binIt{
	my $len=shift;
	$shortest = $len > $shortest ? $shortest : $len;
	$longest = $len > $longest ? $len : $longest;
	foreach my $val (@sequence){
		if ($len <= $val){
			$bins{$val}++;
			last;
		}
	}
	return;
}
