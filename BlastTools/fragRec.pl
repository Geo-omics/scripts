#!/user/bin/perl

=head2 Description

	fragRec.pl : use blast tabular output to generate a fragment recruitment plot x-axis = ref genome(2nd column of blast tabular output); y-axis = percent-id. This script parses the blast output,  generates an R script and executes it.

=head2 Usage

	fragRec.pl [-f reference genome fasta] [-b tabular blast output] [OPTIONS]

=head3 Options

	-prefix	:	output R-script and image file name; default= Process_id.r and iProcess_id.png
	-start	:	define a specific region in the reference genome
	-stop	:	define a specific region in the reference genome; not required if stop is the end of the genome.
	-p	:	min percent id; default= 0%
	-l	:	min alignment length; default= 100
	-s	:	min bit score; default= 40
	-details	:	To get a list of query names.
	-h	:	help; this screen

=head2 Suggestions/Feedback/Beer

	Sunit Jain, 2016 CCL (sunitj@umich.edu)

=cut

use strict;
use Getopt::Long;

my $wgs;
my $blastOut;
my $prefix=$$;
my $start;
my $stop;
my $minPer=0;
my $minLen=100;
my $minBS=40;
my $details;
GetOptions(
	'f|fasta:s'=>\$wgs,
	'b|blast:s'=>\$blastOut,
	'prefix:s'=>\$prefix,
	'start:i'=>\$start,
	'stop:i'=>\$stop,
	'p|percent:f'=>\$minPer,
	'l|length:i'=>\$minLen,
	's|bitscore:i'=>\$minBS,
	'details'=>\$details,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

my ($min, $max,%seqCoords, %fragments);

if ($start){
	$min=$start;
	$max= ($stop) ? $stop : &num_bases($wgs);
}
else{
	$min=1;
	$max=&num_bases($wgs);
}
my $rFile=$prefix.".r";
my $png=$prefix.".png";


&parseBlastOut;

&plot;

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
}

sub parseBlastOut{
	open(BLAST, $blastOut) || die "[ERROR]: $!\n";
	my $prevQuery="";
	while(my $line=<BLAST>){

		next if ($line=~ m/^\#/);
		chomp($line);
		$line=~ s/\r//;

		my($thisQuery, $ref, $pid, $alnLen, $mm, $gaps, $qStart, $qEnd, $sStart, $sEnd, $eval, $bs)=split(/\t/, $line);
		next if ($ref eq $thisQuery);
		next if ($pid < $minPer);
		next if ($bs < $minBS);
		next if ($alnLen < $minLen);
		next if ($prevQuery eq $thisQuery);
		$prevQuery=$thisQuery;

		my($START, $END)=sort{ $a <=> $b }($sStart, $sEnd);
		$START= $START + $seqCoords{$ref};
		$END= $END + $seqCoords{$ref};

		if ($start){
			if ($START > $min && $END < $max){
				print $thisQuery."\n" if ($details);
				$fragments{$thisQuery}{'start'}=$START;
				$fragments{$thisQuery}{'end'}=$END;
				$fragments{$thisQuery}{'perc'}=$pid;
			}
			else{
				next;
			}
		}
		else{
			print $thisQuery."\n" if ($details);
			$fragments{$thisQuery}{'start'}=$START;
			$fragments{$thisQuery}{'end'}=$END;
			$fragments{$thisQuery}{'perc'}=$pid;
		}
	}
}

sub plot{
	#generate R script
	print STDERR "\nCreating R script...";
	my $script;
	$script .= <<EOF;
library("maptools")
EOF

	#R: initialise the png output
	$script .= <<EOF;
png(\"$png\", width = 1920, height = 1200, res=200)
EOF

	#R: set plot title and initialise plot area
	my $plot_title = "% ID";
	my $ylim = 100;
	$script .= <<EOF;
par(xpd=TRUE)
newX=c(1,2)
newY=c(3,4)
plot(newX, newY, type="n", ylab = c("$plot_title"), xlab = c("Reference Genome"), xlim=c($min, $max), ylim=c($minPer, $ylim))
EOF

	#R: draw a line for each read
  foreach my $frags (keys (%fragments)) {
    $script .= <<EOF;
lines(c($fragments{$frags}{'start'}, $fragments{$frags}{'end'}), c($fragments{$frags}{'perc'}, $fragments{$frags}{'perc'}), col=c("firebrick1"))
EOF
  }

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
