#! /usr/bin/perl/

=head1 USAGE

	perl findStretchesOfNs.pl -f infile.fasta -out outfile.fasta [OPTIONS]

=head2 OPTIONS

	min	[integer]	minimum stretch of N's	(default: 10)
	max	[integer]	maximum number of N's	(default: 500)

=head1 Author

	Sunit Jain, January 2013
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my $fasta;
my $out=$$.".out";
my $min=10;
my $max=500;

GetOptions(
	"f|fasta:s"=>\$fasta,
	"o|out:s"=>\$out,
	"min:i"=>\$min,
	"max:i"=>\$max,
	"h|help"=>sub{system('perldoc', $0); exit;},
);

open(FASTA, $fasta)|| die $!;
open(OUT, ">".$out)|| die $!;
my %uniqueSeqNames;
$/=">";
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($header, @sequence)=split("\n", $line);
	my $seq=uc(join("",@sequence));
	my $seqLen=length($seq);

	my @matches;

	while($seq=~ /N{$min,$max}/ig){
#		push (@matches,[ $-[0], $+[0] ]);
		print OUT $header."\t".$-[0]."\t".$+[0]."\n";
		$uniqueSeqNames{$header}++;
	}
}
$/="\n";
close FASTA;
close OUT;
print "# Number of Sequence with at least 1 stretch of Ns >= $min:\t".keys(%uniqueSeqNames)."\n";
