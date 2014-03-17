#!/usr/bin/perl

=head1 DESCRIPTION

	blastDensityPlot.pl -- Given multiple blasts to the same reference database, make density plot based on tabular blast outputs.
	
=head2 NOTE

	* Requires R with library ggplot2
	* Recommended that the blast outputs only be top hits. (Use "top5.pl")

=head1 USAGE

	perl blastDensityPlot.pl -b path_to_all_blast_files -e extension_of_all_blast_files -prefix output_prefix

=head2 Options

	-blast	-b	<CHAR>	path to blast folder
	-ext	-e	<CHAR>	extension used
	-prefix	-p	<CHAR>	prefix for the output files
	
	-perc		<BOOLEAN>	Use Percent Identity [DEFAULT]
	-evalue		<BOOLEAN>	Use E-Values
	-aln		<BOOLEAN>	Use Alignment Length
	-bs			<BOOLEAN>	Use Bit Scores
	-scaled		<BOOLEAN>	Scale the graphs [0,1], such that all the peaks have the same height.
	-tmp		<BOOLEAN>	Keep intermediate files. Deletes the r-script and temporary output unless this flag is provided.
	
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Fri Mar 14 11:49:03 EDT 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;

my($blast, $ext, $prefix);
my($perc,$bs,$aln,$eval, $scaled, $keep_tmp);
my $help;
my $version="blastDensityPlot.pl\tv0.1.0";
GetOptions(
	'b|blast:s'=>\$blast,
	'e|ext:s'=>\$ext,
	'p|prefix:s'=>\$prefix,
	
	'perc'=>\$perc,
	'bs'=>\$bs,
	'evalue'=>\$eval,
	'aln'=>\$aln,
	'scaled'=>\$scaled,
	'tmp'=>\$keep_tmp,
	
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

die "# [FATAL] See '-h' for help on how to use the script\n" if (! $blast);
my ($col, $col_name);
if($perc){ $col = 0; $col_name="Percent_ID"}
elsif($aln){ $col = 1; $col_name="Alignment_Length"}
elsif($bs){$col = -1; $col_name="Bit_Score"}
elsif($eval){$col = -2; $col_name="E_Value"}
else{$col = 0; $col_name="Percent_ID"}

my @files=<$blast/*.$ext>;
die "# [FATAL] Can't find \"$blast\"\n# Please check that the path exist or that you have sufficient privilages.\n" if (scalar(@files)==0);

my %REF;
my $table=$prefix.".out";
open(OUT, ">".$table) || die $!;
print OUT "NAME\t$col_name\tSamples\n";
foreach my $file(@files){
	my $FH;
	print "# Reading $file...\t";
	my $name2print=fileparse($file, ".".$ext);
	open($FH, "<".$file) || die $!;
	while(my $line=<$FH>){
		next if ($line=~ /^#/);
		chomp $line;
		next unless $line;
		
		my($query, $subj, @data)=split(/\t/, $line);
		
		my $column=$data[$col];
		print OUT $subj."__".$name2print."\t".$column."\t".$name2print."\n";
	}
	close $FH;
	print "DONE!\n";
}
close OUT;

&plot;

sub plot{
	my $rScript=$prefix.".r";
	my $pdf=$prefix.".pdf";
	my ($script, $ggplot_cmd);	
	if($scaled){
		$ggplot_cmd="ggplot(myTable, aes(x=".$col_name.", color=Samples, y=..scaled..), size=1) + stat_density(position=\"identity\", fill=NA) + theme_bw()";
	}
	else{
		$ggplot_cmd="ggplot(myTable, aes(x=".$col_name.", color=Samples)) + stat_density(position=\"identity\", fill=NA) + theme_bw()";
	}
	$script.=<<EOF;
library(ggplot2)	
myTable=read.table(\"$table\", header = TRUE, row.names=NULL)
pdf(\"$pdf\", paper="a4r")
$ggplot_cmd
dev.off()
EOF

	#write and execute R script
	open(RFILE, ">".$rScript)|| die "[ERROR] $!\n";
	print RFILE $script;
	close RFILE;
	print "# Starting R...\t";
	system("R --no-save < $rScript > /dev/null");
	print "DONE!\n";
	
	# Remove intermediate files.
	unless($keep_tmp){
		unlink $rScript;
		unlink $table;
	}
}
