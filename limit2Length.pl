#!/user/bin/perl -w
=head1 DESCRIPTION

	Check length of each sequence and only print those that clear the minimum threshold of x (default=60) bases;

=head2 USAGE
	
	perl limit2Length.pl -f fasta_File -len # [-o output file]
	OR
	perl limit2Length.pl -fq fastq_File -len # [-o output file]

=head2 AUTHOR

	Sunit Jain, Mar, 2012

=cut


#######################
## MODULES
#######################
use strict;
use Getopt::Long;

#######################
## PARAMETERS
#######################
my $fasta;
my $fastq;
my $setLen=2000;
my $out=$$.".fa";
my $quick;
my $binSize=2000;
my $lower;

GetOptions(
	'f:s'=>\$fasta,
	'fq:s'=>\$fastq,
	'len:i'=>\$setLen,
	'o|out:s'=>\$out,
	'b|bin:i'=>\$binSize,
	'quick'=>\$quick,
	'below'=>\$lower,
	'h|help'=>sub{system("perldoc", $0); exit;}
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

open (OUT, ">".$out)|| die "[ERROR: $0]: $!";

$/=$fasta ? "\>" : "\n";
my $count=0;
my $total=0;
my $sumLen=0;
my $revisedSumLen=0;
my (@sequence, %bins, $shortest, $longest);
#######################
## MAIN
#######################
unless($quick){
	@sequence= map { $binSize * $_ } 1 .. 50;
	@sequence=sort{$a <=> $b} @sequence;
	push(@sequence, (($sequence[-1])*10)); # capture REALLY long contigs.

	$shortest=1000000000000000000000; # ridiculously high number to capture the smaller ones
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


print "# Total Sequences:\t$total\n";
print "# Total Bases assembled:\t$sumLen\n";
$lower ? print "# Total Bases assembled in contigs <= $setLen:\t$revisedSumLen\n" : print "# Total Bases assembled in contigs > $setLen:\t$revisedSumLen\n";
print "# Longest Sequence Length:\t$longest\n" unless($quick);
print "# Shortest Sequence Length:\t$shortest\n" unless($quick);
print "# Mean Sequence Length before setting the $setLen base limit:\t$meanLenBefore\n";
$lower ? print  "# Number of Sequences with Length equal to or less than $setLen bases:\t$count\n" : print "# Number of Sequences with Length greater than $setLen bases:\t$count\n";
print "# Mean Length after setting the $setLen base limit:\t$meanLenAfter\n";


unless($quick){
	my $i=0;
	foreach my $b(@sequence){
		my $from= $i > 0 ? ($sequence[$i-1])+1 : 1;
		print $from."-".$b."\t".$bins{$b}."\n" if $bins{$b};
		$i++;
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

		if(($lower) && ($len <= $setLen)){
			print OUT "$seqDesc\n$seq\n\+\n$qual\n";
			$revisedSumLen+=$len;
			return 1;
		}
		elsif((! $lower) &&  ($len > $setLen)){
			print OUT "$seqDesc\n$seq\n\+\n$qual\n";
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
		return 1;
	}
	elsif((! $lower) &&  ($len > $setLen)){
		print OUT ">".$line."\n";
		$revisedSumLen+=$len;
		return 1;
	}
	else{
		return 0;
	}
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
