#!/usr/local/bin/perl

=head1 DESCRIPTION

gbk2fna.pl -- Read Genbank to Nucleotide Fasta file.

=head1 USAGE

perl gbk2fna.pl

=head2 Options

	-in	<CHAR>	FASTA File
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Thu Oct  1 11:10:47 EDT 2015)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use FileHandle;
use File::Basename;

my $help;
my $version=fileparse($0)."\tv0.0.1b";
my $infile;
my $outfile;
GetOptions(
	'in:s'=>\$infile,
	'out:s'=>\$outfile,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

use Bio::SeqIO;
my $seq_in = Bio::SeqIO->new(
	 -file   => "<$infile",
	-format => "genbank",
);
my $seq_out = Bio::SeqIO->new(
	-file   => ">$outfile",                                                                                                                                              -format => "fasta",                                                                                                );
while (my $inseq = $seq_in->next_seq) {
                                                                                                                                                                                        $seq_out->write_seq($inseq);
                                                                                                                                                                                        }
