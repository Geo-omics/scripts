#!/usr/bin/perl

=head1 DESCRIPTION

reverse_complement.pl -- Do this.

=head1 USAGE

perl reverse_complement.pl -fasta input.fasta -out output.fasta

=head2 Options

	-fasta -f	<CHAR>	Fasta file
	-out  -o	<CHAR>	Output file name
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Tue Jan 20 11:00:25 EST 2015)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;

my $help;
my $version="reverse_complement.pl\tv0.0.1b";
my ($fasta, $out);
GetOptions(
	'f|fasta:s'=>\$fasta,
	'o|out:s'=>\$out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

open(FASTA, "<".$fasta) || die $!;
open(OUT, ">".$out)|| die $!;
$/=">";
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;

	my($header, @seq)=split(/\n/, $line);
	my $s=join("",@seq);
	my $rcSeq=reverseComplement($s);

	print OUT ">".$header."\n".$rcSeq."\n";
}
$/="\n";
close FASTA;
close OUT;

sub reverseComplement{
	my $seq=shift;
	chomp $seq;
	my $rSeq=uc(reverse($seq));
	$rSeq=~ tr/GTCA/CAGT/;
	return $rSeq;
}
