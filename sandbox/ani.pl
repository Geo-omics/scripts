#!/usr/bin/perl

=head1 DESCRIPTION

ani.pl -- This program calculates Average Nucleotide Identity (ANI) for a pair of bacterial genomes.
NOTE: Requires an older version of blast. Tested on version 2.2.20
=head1 USAGE

perl ani.pl -blast path to blast binaries -query genome -subj the other strain genome -outdir output directory
	
=head2 Example

	perl ANI.pl -blast ./blast-2.2.20/bin/ -query strain1.fa -subj strain2.fa -outdir result

=head2 Options

    query	<CHAR>	Query strain genome sequence in FASTA format. [REQUIRED]
    subj	<CHAR>	Subject strain genome sequence in FASTA format. [REQUIRED]
    
    outdir	<CHAR>	output directory
    threads	<INT>	number of processors to use for blast

    -version -v	<BOOLEAN>	version of the current script
    -help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Tue May 13 16:07:06 EDT 2014)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use strict;
use List::Util qw(min);
use Cwd qw(abs_path);

my $help;
my $version="ani.pl\tv0.0.3b";
my ($q, $s, $outDir, $path_to_blast_bin, $num_threads);
GetOptions(
    "query:s" => \$q,
    "subj:s" => \$s,
    "o|outdir:s" => \$outDir,
    "t|threads:s"=>\$num_threads,
    'v|version'=>sub{print $version."\n"; exit;},
    'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

if(! $q || ! $s){ system("perldoc $0 \| cat"); exit; }

my $query=abs_path($q);
my $subj=abs_path($s);
$outDir=$q."_vs_".$s if (! $outDir);
unless(-d $outDir){`mkdir $outDir`;}
$num_threads=4 if ( ! $num_threads);

#Split query genome and write segments in $outDir/Query.split
my $chop_len = 1020; 
$/ = "\n>";
open (QUERY, "<".$query) || die "$query $!\n";
open (CHOP,">".$outDir."/Query.split") || die "$outDir/Query.Split";
while(<QUERY>){
    chomp;
    s/>//g;
    my ($scaf,$seq) = (split /\n/,$_,2);
    my $scaf_name = (split /\s+/,$scaf)[0];
    $seq =~ s/\s+//g;
    my @cut = ($seq =~ /(.{1,$chop_len})/g);
    for my $cut_num (0..$#cut){
        next if length($cut[$cut_num]) < 100; 
        my $sgmID = "$scaf_name\_$cut_num";
        print CHOP ">$sgmID\n$cut[$cut_num]\n";
    }
}
close QUERY;
close CHOP;
$/ = "\n";

#BLAST alignment
`ln -sf $subj $outDir/Subject.fa`;
`formatdb -i $outDir/Subject.fa -p F`;
#`makeblastdb -dbtype nucl -in $outDir/Subject.fa`;
#`blastn -task blastn -db $outDir/Subject.fa -query $outDir/Query.split -evalue 1e-15 -xdrop_gap 150 -num_threads $num_threads -penalty -1 -out $outDir/raw.blast -outfmt 6 -dust no`;
`blastall -i $outDir/Query.split -d $outDir/Subject.fa -X 150 -q -1 -F F -e 1e-15 -m 8 -a $num_threads -o $outDir/raw.blast -p blastn`;

#Set Identity and Alignment Percentage cut off following paper of JSpecies
my $pid = 30;
my $cvg_cut = 70;

my ($ANI,%query_best,$sumID,$count);
open (BLAST,"<".$outDir."/raw.blast") || die "raw.blast $!\n";
open (ANI, ">".$outDir."/ani_raw.blast") || die "ani_raw.blast $!\n";
while(<BLAST>){
    chomp;
    my @t = split;
    next if $t[3] < 100;
    next if exists $query_best{$t[0]}; 
    $query_best{$t[0]} = 1; #only use best hit for every query segments
    next if $t[2]<=$pid;
    next if ($t[3]*100/1020) < $cvg_cut;
    print ANI $t[0]."_".$q."\t".$t[2]."\t".$q."_vs_".$s."\n";
    $sumID += $t[2];
    $count++;
}
close BLAST;
close ANI;
$ANI = $sumID/$count;
print "$q\t$s\t$ANI\n";
