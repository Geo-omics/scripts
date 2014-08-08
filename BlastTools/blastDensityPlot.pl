#!/usr/bin/perl

=head1 DESCRIPTION

	blastDensityPlot.pl -- Given multiple blasts to the same reference database, make density plot based on tabular blast outputs.
	
=head2 NOTE

	* Requires Rscript
	* Requires R package 'ggplot2'
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
	-verbose	<BOOLEAN>	Get R log file
	
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

my($blast, $perc,$bs,$aln,$eval, $scaled, $keep_tmp, $verbose, $table, $help, $ani);
my $ext="blastn";
my $prefix=$$;
my $version="blastDensityPlot.pl\tv0.1.5";
GetOptions(
	'b|blast:s'=>\$blast,
	't|table:s'=>\$table,
	'e|ext:s'=>\$ext,
	'p|prefix:s'=>\$prefix,
	'ani:f'=>\$ani,
	
	'perc'=>\$perc,
	'bs'=>\$bs,
	'evalue'=>\$eval,
	'aln'=>\$aln,
	'scaled'=>\$scaled,
	'tmp'=>\$keep_tmp,
	'verbose'=>\$verbose,
	
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

die "# [FATAL] See '-h' for help on how to use the script\n" if (! $blast && ! $table);
my ($col, $col_name);
if($perc){ $col = 0; $col_name="Percent_ID"}
elsif($aln){ $col = 1; $col_name="Alignment_Length"}
elsif($bs){$col = -1; $col_name="Bit_Score"}
elsif($eval){$col = -2; $col_name="E_Value"}
else{$col = 0; $col_name="Percent_ID"}

if($table){
	&plot;
	exit 0;
}

my @files=<$blast/*.$ext>;
die "# [FATAL] Can't find \"$blast\"\n# Please check that the path exist or that you have sufficient privilages.\n" if (scalar(@files)==0);

my %REF;
$table=$prefix.".out";
open(OUT, ">".$table) || die $!;
print OUT "NAME\t$col_name\tSamples\n";
foreach my $file(@files){
	my ($FH, %seen);
	print "# Reading $file...\t";
	my $name2print=fileparse($file, ".".$ext);
	open($FH, "<".$file) || die $!;
	while(my $line=<$FH>){
		next if ($line=~ /^#/);
		chomp $line;
		next unless $line;
		
		my($query, $subj, @data)=split(/\t/, $line);
		
		next if $seen{$query};
		print OUT $subj."__".$name2print."\t".$data[$col]."\t".$name2print."\n";
		$seen{$query}++;
	}
	close $FH;
	print "DONE!\n";
}
close OUT;

&plot;
exit 0;

sub plot{
	my $pdf=$prefix."_".lc($col_name).($scaled ? "_scaled" : "").($ani ? "_withANI" : "").".pdf";
	if (-e $pdf){ $pdf=$prefix."_".lc($col_name).($scaled ? "_scaled" : "")."_".$$.".pdf"}
	my $ggplot_cmd;
	if($scaled){
		$ggplot_cmd="ggplot(myTable, aes(x=".$col_name.", color=as.factor(Samples), y=..scaled..), size=1) + stat_density(position=\"identity\", fill=NA) + theme_bw() + scale_colour_discrete(name  =\"Samples\")";
	}
	else{
		$ggplot_cmd="ggplot(myTable, aes(x=".$col_name.", color=as.factor(Samples))) + stat_density(position=\"identity\", fill=NA) + theme_bw() + scale_colour_discrete(name  =\"Samples\")";
	}

	if ($ani){
		$ggplot_cmd.="+ geom_vline(aes(xintercept=$ani), show_guide=F) +geom_text(aes(x=$ani, y=0, label=\"ANI\"), color=\"black\", angle=90, vjust=0, hjust=0)";
	}
	my $save_pdf="ggsave(file=\"$pdf\", width = 210, height = 297, units = \"mm\")";

	#Create command line Rscript
	print "# Creating the density plot...\t";
	my $script="-e 'library(ggplot2)' -e 'myTable=read.table(\"$table\", header = TRUE, row.names=NULL)' -e 'myPlot=$ggplot_cmd' -e '$save_pdf'";
	if ($verbose){
		$script.=" &> ".$prefix.".log";
	}
	else{
		$script.=" &> /dev/null";
	}
	
	# Execute Rscript
	system("Rscript ".$script);
	print "DONE!\n";
	
	# Remove intermediate files.
	unless($keep_tmp){
		unlink $table;
	}
}
