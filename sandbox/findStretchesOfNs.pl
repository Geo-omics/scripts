#! /usr/bin/perl/

=head1 USAGE

	perl findStretchesOfNs.pl -f infile.fasta -coord outfile.tsv [OPTIONS]

=head2 OPTIONS

	-nuc	[character]	any character other than 'N'	(default: N)
	-coord	[character]	outputs coordinates file, start and stop position of the stretch of N's
	-min	[integer]	minimum stretch of N's	(default: 10)
	-max	[integer]	maximum number of N's	(default: 500)
	-split	[character]	output file; split the fasta sequence at N's and write to this file.
	-minLen	[integer]	When splitting on Ns, what's the smallest sequence you wish to obtain. (Default: 100 (bases))
	
=head3 Example

	perl findStretchesOfNs.pl -f test.fasta -o test.coord -min 1 -split test_split_l20.fasta -len 20

=head1 Author

	Sunit Jain, January 2013
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my $version="v0.3.7";
my $fasta;
my $coord;
my $nuc="N";
my $min=10;
my $max=500;
my $split;
my $minLen=100;
GetOptions(
	"f|fasta:s"=>\$fasta,
	"c|o|coord:s"=>\$coord,
	"n|nuc:s"=>\$nuc,
	"split:s"=>\$split,
	"min:i"=>\$min,
	"max:i"=>\$max,
	"len:i"=>\$minLen,
	"h|help"=>sub{system('perldoc', $0); exit;},
);
print "# findStretchesOfNs.pl $version\n";
die "[ERROR] Choose an output type, either '-coord' or '-split'. See '-h' for more options.\n" if(! $coord) && (! $split);

my $find=quotemeta "$nuc";
open(FASTA, $fasta)|| die $!;
if ($coord){
	open(OUT, ">".$coord)|| die $!;
	print OUT "# Sequence_Name\tN_Start\tN_Stop\tN_Length\n";
}

if ($split){
	open(SPLIT, ">".$split)|| die $!;
}
my %uniqueSeqNames;
$/=">";
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($header, @sequence)=split("\n", $line);
	my $seq=uc(join("",@sequence));
	my ($realHeader, @description)=split(/\s+/, $header);
	my $desc=join(" ", @description);
	my $parts=0;
	my $lastOffset=0;
	while($seq=~ /($find){$min,$max}/ig){
		print OUT $header."\t".$-[0]."\t".$+[0]."\t".($+[0]-$-[0])."\n";
		if($split){
			my $currentOffset=$-[0];
			print "Last Offset:\t$lastOffset\tCurrent Offset:\t$currentOffset\n";
			my $subSeq = substr($seq, $lastOffset, ($currentOffset - $lastOffset));
			if(length($subSeq) >= $minLen){
				print SPLIT ">".$realHeader."_".$parts." ".$desc."\n";
				print SPLIT $subSeq."\n";
				$parts++;
			}
			$lastOffset=$+[0];
		}
		$uniqueSeqNames{$header}++;
	}
	print SPLIT ">".$realHeader."_".$parts." ".$desc."\n";
	print SPLIT substr($seq, $lastOffset)."\n";
}
$/="\n";
close FASTA;
close OUT if($coord);
close SPLIT if($split);
print "# Number of Sequence with at least 1 stretch of Ns >= $min:\t".keys(%uniqueSeqNames)."\n";
