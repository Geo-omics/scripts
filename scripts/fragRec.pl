#!/usr/bin/perl

=head1 Description

	fragRec.pl : use multiple blast outputs in tabular format to generate a fragment recruitment plot.
	-- x-axis = ref genome(2nd column of blast tabular output);
	-- y-axis = percent-id. This script parses the blast output;
	-- Generates an R script and executes it. 
	-- All the blasts must have the same reference database

=head1 Usage

	fragRec.pl [-f reference genome fasta] [-b path to folder containing the blast outputs] [-e extension for the blast outputs] [OPTIONS]

=head2 Options

	-prefix	:	output R-script and image file name; default= Process_id.r and iProcess_id.png
	-start	:	define a specific region in the reference genome
	-stop	:	define a specific region in the reference genome; not required if stop is the end of the genome.
	-p	:	min percent id; default= 0%
	-l	:	min alignment length; default= 100
	-s	:	min bit score; default= 40
	-details	:	To get a list of query names.
	-h	:	help; this screen

=head1 Author

	Sunit Jain, (Fri July 25 11:49:03 EDT 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;

my ($wgs, $ext,$blast, $details,$start, $stop, $by_aln);
my $prefix=$$;
my $minPer=0;
my $minLen=100;
my $minBS=40;
my $version="fragRec.pl version 0.6.5";
GetOptions(
	'f|fasta:s'=>\$wgs,
	'b|blast:s'=>\$blast,
	'e|ext:s'=>\$ext,
	'prefix:s'=>\$prefix,
	'start:f'=>\$start,
	'stop:f'=>\$stop,
	'p|percent:f'=>\$minPer,
	'l|length:i'=>\$minLen,
	's|bitscore:i'=>\$minBS,
	'details'=>\$details,
	'by_aln'=>\$by_aln,
	'h|help'=>sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print "# $version\n"; exit;},
);

## Sanity Checks ##

my ($min, $max,%seqCoords, %fragments, @limitValues);

if ($start){
	$min=$start;
}
else{
	$min=0;
}

if($stop){
	$max=$stop;
	&num_bases($wgs);
}
else{
	$max=&num_bases($wgs);
}

my $rFile=$prefix.".r";
my $pdf=$prefix.".pdf";

my @files=<$blast/*.$ext>;
die "# [FATAL] Can't find \"$blast\"\n# Please check that the path exist or that you have sufficient privilages.\n" if (scalar(@files)==0);
foreach my $file(@files){
	&parseBlastOut($file);
}

&plot;

## Sub-Routines ##

sub num_bases{
	my $fasta=shift;
	my $pos=1;
	open(FASTA, $fasta) || die "[ERROR]: $!\n";
	$/=">";
	while(my $line=<FASTA>){
		next if ($line=~ m/^\#/);
		chomp($line);
		$line=~ s/\r//g;

		my($header, @sequence)=split(/\n/, $line);
		my $seq= join("", @sequence);
		my $seqLen=length($seq);
		my ($head, @etc)=split(/ /, $header);
		$seqCoords{$head}=$pos;
		$pos+=$seqLen;
	}
	$/="\n";
	return($pos);
	close FASTA;
}

sub parseBlastOut{
	my $blastOut=shift;
	my (%seen, $FH);
	my $name2print=fileparse($blastOut, ".".$ext);
	print "# Reading\t".$name2print."\t...\n";
#	my $fileLimitValue= ($by_aln) ? 0 : 101;
	my $color=getRandomColor();
	open($FH, "<".$blastOut) || die "[ERROR]: $!\n";
	while(my $line=<$FH>){

		next if ($line=~ m/^\#/);
		chomp($line);
		$line=~ s/\r//;

		my($query, $ref, $pid, $alnLen, $mm, $gaps, $qStart, $qEnd, $sStart, $sEnd, $eval, $bs)=split(/\t/, $line);
		next if ($ref eq $query);
		next if ($pid < $minPer);
		next if ($bs < $minBS);
		next if ($alnLen < $minLen);
		next if $seen{$query};
		$seen{$query}++;

		my $value;
		if($by_aln){
			$value=$alnLen;	
		}
		else{
			$value=$pid;
		}
		push(@limitValues, $value);

		my($START, $END)=sort{ $a <=> $b }($sStart, $sEnd);
		$START= $START + $seqCoords{$ref};
		$END= $END + $seqCoords{$ref};
		my $name=$query."__".$name2print; # so the query name is unique for each file.

		if(($END <= $min) || ($max <= $START)){
#		Case4: When alignment length is out of scope; ignore
			next;
		}
		elsif(($START <= $min) && ($max <= $END)){
#		Case1: When alignment is longer than the selected region; span the entire region
			$START=$min;
			$END=$max;
		}
		elsif(($START <= $min) && ($END <= $max)){
#		Case2: When alignment starts before the selected region but ends within the 
#			selected region; change the start to match the start of selected region
			$START=$min;
		}
		elsif(($min <= $START) && ($max <= $END)){
#		Case3: When alignment ends after the selected region but starts within the 
#			selected region; change the end to match the end of selected region
			$END=$max;
		}	

		print $query."\n" if ($details);
		$fragments{$name}{'start'}=$START;
		$fragments{$name}{'end'}=$END;
		$fragments{$name}{'y_value'}=$value;
		$fragments{$name}{'color'}=$color;
	}
	close $FH;
}

sub getRandomColor{
    my ($r, $g, $b) = map { int rand 256 } 1 .. 3;
	my $color=sprintf ("#%2.2X%2.2X%2.2X",$r, $g, $b);
    return ($color);
}

sub plot{
	#generate R script
	print STDERR "\nCreating R script...";
	my $script;
	$script .= <<EOF;
library("maptools")
library("sp")
gpclibPermit()
EOF

	#R: initialise the png output
	$script .= <<EOF;
pdf(\"$pdf\", paper="a4r")
EOF

	#R: set plot title and initialise plot area
	my($plot_title, $ylim, $minY);
	my @sorted = sort { $a <=> $b } @limitValues;
	if ($by_aln){
		$plot_title = "Alignment Length";
		$ylim = $sorted[-1];
		$minY = $sorted[0];
	}
	else{
		$plot_title="% ID";
		$ylim=100;
		$minY=$sorted[0];
	}
	$script .= <<EOF;
par(xpd=TRUE)
newX=c(1,2)
newY=c(3,4)
plot(newX, newY, type="n", ylab = c("$plot_title"), xlab = c("Reference Genome"), xlim=c($min, $max), ylim=c($minY, $ylim))
EOF

	#R: draw a line for each read
	my %colors;
  foreach my $frags (keys (%fragments)) {
    $script .= <<EOF;
lines(c($fragments{$frags}{'start'}, $fragments{$frags}{'end'}), c($fragments{$frags}{'y_value'}, $fragments{$frags}{'y_value'}), col=c("$fragments{$frags}{'color'}"))
EOF
	my($fragment, $file)=split(/__/, $frags);
	$colors{$file}=$fragments{$frags}{'color'};
  }

	my ($fileNames, $fileColors);
	foreach my $f(keys %colors){
		$fileNames.="\'".$f."\',";
		$fileColors.="\'".$colors{$f}."\',";
	}
	$fileNames=~ s/\,$//;
	$fileColors=~ s/\,$//;

#	$script.=<<EOF;
#legend("topright", legend=c($fileNames), col=c($fileColors), bty="n", inset=c(($max+10), $ylim), cex=2, pch=15, xpd=TRUE)
#EOF

  #R: close device
  $script .= <<EOF;
dev.off()
EOF

	#write and execute R script
	open(RFILE, ">".$rFile)|| die "[ERROR] $!\n";
	print RFILE $script;
	close RFILE;
	print STDERR "\nStarting R...";
	system("R --no-save < $rFile > /dev/null");

}
