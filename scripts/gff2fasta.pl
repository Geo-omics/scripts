#!/usr/bin/perl

=head1 DESCRIPTION

gff2fasta.pl -- Do this.

=head1 USAGE

	perl gff2fasta.pl -gff file.gff -contig file.fasta -genes output_file.fasta

=head2 Options
	-prefix	-p	<CHAR>	Uniq project identifier added before sequence header.
	-genes	-g	<CHAR>	output fasta file
	-gff		<CHAR>	input GFF file
	-contig	-c	<CHAR>	Contigs/Scaffolds Input
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Fri Mar  6 12:18:22 EST 2015)
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

use Bio::SeqIO;

my $help;
my $version=fileparse($0)."\tv0.1.1";
my $verbose = 0;
my ($quiet, $contigFile, $prefix, $gffFile, $outFile);
GetOptions(
	'prefix:s'=>\$prefix,
	'g|genes:s'=>\$outFile,	
	'gff:s'=>\$gffFile,
	'c|contig:s'=>\$contigFile,
	'q|quiet'=>\$quiet,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

## read in gff 
warn "reading GFF\n";
 
my %gff;
my $GFF=FileHandle->new();
open ($GFF, '<', $gffFile) || die "fail\n";
 
while(<$GFF>){
  my ($seqid, undef, undef, $start, $end,
      undef, undef, undef, $attrs) = split(/\t/);
 
  push @{$gff{$seqid}}, [$start, $end, $attrs];
}
 
warn "OK\n";
 
 
 
## Do the fasta
my $GENES=FileHandle->new();
open( $GENES, ">", $outFile) || die $!;
$verbose=-1 if ($quiet);
my $seqio = Bio::SeqIO->new(	-file => $contigFile,
								-format => 'fasta',
								-verbose=> $verbose )
							|| die "double fail\n";

while(my $sobj = $seqio->next_seq){
  my $seqid = $sobj->id;
 
  unless(defined($gff{$seqid})){
    warn "no features for $seqid\n";
    next;
  }
 
  my $seq = $sobj->seq;
 
  for(@{$gff{$seqid}}){
    my ($start, $end, $attrs) = @$_;
 
    warn join("\t", $start, $end, $attrs), "\n"
      if $verbose > 0;
 
    my %attrs = split(/=|;/, $attrs);
 
    print $GENES ">".($prefix ? $prefix."_" : "").$seqid."-". $attrs{"ID"}.
      "/$start-$end (". ($end-$start+1). ")\n";
 
    print $GENES substr($seq, $start, $end-$start+1), "\n";
  }
 
  #exit;
}
close $GENES; 
warn "OK\n";
