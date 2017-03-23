#! /usr/bin/perl

=head1 Description

	This program takes a class file, a names file, a fasta file and extracts the seqs for the contigs presents in the desired class.

=head1 Usage

	perl getClassFasta.pl -cls <CLASS File> -names <NAMES File> -fasta <esom FASTA File> -num <CLASS NUMBER>

=head2 Options

	-cls	<STRING>	the .cls file
	-names	<STRING>	the .names file
	-fasta	<STRING>	The Concatenated fasta file that was used with the tetramer script.
	-num	<INTEGER>	the class number you're interested in.
	-id	<STRING>	An identifier that is unique to the set of contigs you're interested in.

=head2 Experimental Feature

	-loyal	<REAL: 1-100>	Bin loyalty value; Only get a contig if loyal% bins in this class; Default = 0

=head1 Questions/Comments/Suggestions/Accolades/Beer

	Sunit Jain, sunitj [AT] umich [DOT] edu

=cut


use strict;
use Getopt::Long;
use File::Basename;

my ($class, $names, $fasta);
my $classNum=0;
my $confidence=0;
my $version="getClassFasta.pl\tv0.1.2";
my ($no_conf, $unique_id);

GetOptions(
	"cls=s"=>\$class,
	"names=s"=>\$names,
	"fasta=s"=>\$fasta,
	"num=i"=>\$classNum,
	"loyal=f"=>\$confidence,
	"no_conf"=>\$no_conf,
	"id:s"=>\$unique_id,
	'v|version'=>sub{print $version."\n"; exit;},
	"help"=>sub{system("perldoc", $0); exit;},
);

print "# $version\n";

print "CLASS: $class\nNAMES: $names\nFASTA: $fasta\nCLASS_NUM: $classNum\n";
die if (! $class or ! $names or ! $fasta);

$confidence=($confidence/100);

# Parse the *.names file to id the seq names in the fasta file using the classHash from above.
my %seqNames;
open (NAMES, $names) || die "ERROR: $names.\n".$!;
	while (my $line=<NAMES>){
		chomp($line);
		unless ($line =~ /^%/){
			my ($seqNum, $seqSplitName, $seqContigName)=split(/\t/, $line);
			push(@{$seqNames{$seqContigName}}, $seqNum); # %seqNames {Name of the whole contig} => @(SeqNums for each window)
		}
	}
close NAMES;

$classNum++;
# Parse *.cls file to get SeqID for all Seqs in the desired class
my %clsHash;
open ( CLS, $class) || die "ERROR: $class.\n".$!;
	while (my $line=<CLS>){
		chomp($line);
		unless ($line=~ /^%/){
			my ($seqNum, $cls)=split(/\t/,$line);
			$cls++;
			$clsHash{$seqNum}=$cls;	# %clsHash {Sequence Number  => Class Number}
		}
	}
close CLS;

my %currentCls;
foreach my $contig(keys %seqNames){
	my %cls;
	my $windows=0;
	
#	next if ($contig=~ /^gi/);
	if ($unique_id){
		next unless ($contig=~ /$unique_id/)
	}
	
	foreach my $seqNum(@{$seqNames{$contig}}){
		$cls{$clsHash{$seqNum}}++;
		$windows++;
	}
	
	next if (! $cls{$classNum});
	my $cls_conf=$cls{$classNum}/$windows;
	
	next if ($cls_conf < $confidence);
	$currentCls{$contig}=sprintf( "%.4f", $cls_conf);
}
undef %seqNames;
undef %clsHash;

# Parse the fasta file to get the desired sequences, using the seqNames hash from above.
$classNum--;
my $outFile=$classNum.".fasta";
my $conf_out=$classNum.".conf";
open (OUT,">".$outFile )|| die "ERROR: $outFile\n".$!;
open (FASTA, $fasta) || die "ERROR: $fasta\n".$!;
unless ($no_conf){
	warn "[WARNING] The loyalty value was set as '0', this might lead to duplication of scaffolds in your bins.
	Please refer to the $conf_out file for details about your recruited contigs\n" if ($confidence==0);
	open (CONF,">".$conf_out) || die "ERROR: $conf_out\n".$!;
	print CONF "# Min Conf: ".$confidence."\n";
	print CONF "# Bin\tContig_Name\tConfidence\n";
}

$/= ">";
while (my $line = <FASTA>) {
	chomp $line;
	next unless $line;
	my ($name, @sequence) = split (/\n/, $line);
	my $seq = uc(join ("", @sequence));
	if ($currentCls{$name}){
		print OUT ">".$name."\n";
		print OUT $seq."\n";
		print CONF $classNum."\t".$name."\t".$currentCls{$name}."\n";
	}
}
$/= "\n";
close FASTA;
close OUT;
close CONF;

