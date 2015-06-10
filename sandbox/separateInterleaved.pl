#!/usr/bin/perl

=head1 DESCRIPTION

	separateInterleaved.pl -- Separate interleaved files. Tested for illumina datasets only.

=head1 USAGE

	perl separateInterleaved.pl -fq interleaved.fastq -prefix split

=head2 Options

	-fq	OR -fastq <CHAR>	interleaved Fastq file. [REQUIRED]
	-p OR -prefix	<CHAR>	prefix for the forward and reverse files. [Default=prefix of the input fastq]

	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Wed Jan 22 16:16:19 EST 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;

my $help;
my $version="separateInterleaved.pl\tv0.0.1b";
my ($fastq, $prefix);
GetOptions(
	'fq|fastq:s'=>\$fastq,
	'o|p|prefix:s'=>\$prefix,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my @ext=qw(fastq fq);
if(! $prefix){
	$prefix=fileparse($fastq, @ext);
}

my $fwd=$prefix.".fwd.fastq";
my $rev=$prefix.".rev.fastq";

open(FASTQ, "<".$fastq) || die $!;
my $currentFile=$fwd;
while(my $line=<FASTQ>){
	my @seqData=();
	push (@seqData,$line);
	for (my $i=0; $i<3; $i++){
		$line=<FASTQ>;
		push(@seqData,$line);
	}
	
	print_to_correct_file($currentFile, \@seqData);
	
	if($currentFile eq $fwd){
		$currentFile=$rev;
	}
	elsif($currentFile eq $rev){
		$currentFile=$fwd;
	}
}
close FASTQ;

sub print_to_correct_file{
	my $currentFile=shift;
	my $seqData=shift;
	my $FH;
	open($FH, ">>".$currentFile) || die $!;
	foreach my $line(@{$seqData}){
		print $FH $line;
	}
	close $FH;
}
