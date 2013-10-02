#! /usr/bin/perl/

=head1 USAGE

	perl findStretchesOfNs.pl -f infile.fasta -coord outfile.fasta [OPTIONS]

=head2 OPTIONS

	-nuc	[character]	any character other than 'N'	(default: N)
	-coord	[character]	coordinates file, start and stop position of the stretch of N's
	-min	[integer]	minimum stretch of N's	(default: 10)
	-max	[integer]	maximum number of N's	(default: 500)
	-split	[character]	output file; split the fasta sequence at N's and write to this file.

=head1 Author

	Sunit Jain, January 2013
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my $fasta;
my $coord=$$.".out";
my $nuc="N";
my $min=10;
my $max=500;
my $split;
GetOptions(
	"f|fasta:s"=>\$fasta,
	"c|o|coord:s"=>\$coord,
	"n|nuc:s"=>\$nuc,
	"min:i"=>\$min,
	"max:i"=>\$max,
	"h|help"=>sub{system('perldoc', $0); exit;},
);

my $find=quotemeta "$nuc";
open(FASTA, $fasta)|| die $!;
open(OUT, ">".$coord)|| die $!;
my %uniqueSeqNames;
$/=">";
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($header, @sequence)=split("\n", $line);
	my $seq=uc(join("",@sequence));
	my $parts=0;
	while($seq=~ /($find){$min,$max}/ig){
		print OUT $header."\t".$-[0]."\t".$+[0]."\n";
		if($split){
			print SPLIT ">".$header."_".$parts."\n";
			print SPLIT $&."\n";
			$parts++;
		}
		$uniqueSeqNames{$header}++;
	}
}
$/="\n";
close FASTA;
close OUT;
print "# Number of Sequence with at least 1 stretch of Ns >= $min:\t".keys(%uniqueSeqNames)."\n";
