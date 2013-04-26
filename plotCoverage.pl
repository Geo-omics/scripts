#!/usr/bin/perl

#generate a binned coverage plot for blast hits on a full genome
#test that the user has R installed and available

use strict;
use Getopt::Long;

my $Rversion = `R --version`;
unless ($Rversion =~ /R\sversion/) {
	print "\nYou must have R installed and available to use this program.";
	exit;
}


my $USAGE=q/ USAGE:
plot_coverage.pl	-r <filename> Reference genome in fasta format.
			-b <filename> Blast output, in -m8 hit table format. You may plot multiple samples by using the -b flag multiple times, e.g. "-b firstsample.blast_output -b secondsample.blast_output". The coverage maps for the different samples will be displayed on the sample plot, overlayed in different colours. A maximum of 5 different samples can be displayed on the same plot.
			-p <string> Prefix for output files.
			-w <integer> Window size for binning. A window size of 0 will simply plot base-by-base coverage.
			-i <number> Plot width, in inches
			-h <number> Plot height, in inches

OPTIONAL:		-f <filename> Features file: a comma-separated file, one feature per line, fields [start pos] [end pos] [colour] [feature name]. If no colour is specified, default is firebrick. If no feature name is specified, default is blank.
			-y Plot the y axis on a log scale.
/;

my ($reference_genome, @blast_output, $window_size, $minPer, $minLen, $minBS, $features_file, $logY);
my $output_prefix=$$;
#get and check options

GetOptions (
'r=s' => \$reference_genome,
'b=s' => \@blast_output,
'prefix=s' => \$output_prefix,
'w=i' => \$window_size,
'p|percent:f'=>\$minPer,
'l|length:i'=>\$minLen,
's|bitscore:i'=>\$minBS,
'f=s' => \$features_file,
'y!' => \$logY,
) or die $USAGE;
die $USAGE if !$reference_genome or scalar(@blast_output) == 0;
die ("ERROR - a maximum of 5 different samples can be displayed on the same plot") if @blast_output > 5;

##BODY
my $plot_width= 11; #in Inches
my $plot_height= 8.5; #in Inches
my ($reference_genome_length, $script, %seqCoords, %window_average, $maxCoverage, %coverage );

&num_bases;
$window_size=100 if (! $window_size);
&get_coverage;
&do_moving_window_average;
&plot;
exit;
##END BODY

##SUBS

#get the length of the reference genome
sub num_bases{
	my $fasta=$reference_genome;
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
	$reference_genome_length=$pos;
	$/="\n";
}

sub get_coverage {

	#loop over multiple blast outputs
	my $coverageTotal;
	foreach my $blastOut (@blast_output) {
		my $totalCoverageForThisBlast;

		die ("ERROR - could not open blast output $blastOut\n") unless open(BLAST, "<$blastOut");
		print STDERR "\n";

		while (my $line = <BLAST>) {

			print STDERR "\nProcessed $. lines of $blastOut..." if $. % 10000 == 0;

			#get the start and end positions of the hit
			chomp $line;
			my($thisQuery, $ref, $pid, $alnLen, $mm, $gaps, $qStart, $qEnd, $sStart, $sEnd, $eval, $bs)=split(/\t/, $line);
			next if ($pid < $minPer);
			next if ($bs < $minBS);
			next if ($alnLen < $minLen);

			my ($START,$END) = sort {$a <=> $b} ($sStart, $sEnd);
			$START= $START + $seqCoords{$ref};
			$END= $END + $seqCoords{$ref};		

			#make sure the start and end positions exist and are numbers
			die ("ERROR - malformed line in blast output $blastOut at line $.\n") unless $START =~ /^\d+$/ && $END =~ /^\d+$/;

			#foreach base i in the reference genome, increment its coverage count and %id count if needed
			for (my $i = $START; $i <= $END; ++$i) {
				++$coverage{$blastOut}{$i};
				++$coverageTotal;
				++$totalCoverageForThisBlast;
			}
		}

		close BLAST;

		my $coverageMeanForThisBlast = $totalCoverageForThisBlast / $reference_genome_length;
		print STDERR "\nMEAN COVERAGE FOR $blastOut: $coverageMeanForThisBlast";
	}

	my $coverageMean = $coverageTotal / $reference_genome_length;
	print STDERR "\nMEAN COVERAGE: $coverageMean";
}

sub do_moving_window_average {

	#loop over multiple blast outputs
	foreach my $blastOut (@blast_output) {

		print STDERR "\nComputing moving window average for $blastOut...";
		for (my $i = 1; $i <= $reference_genome_length; $i += $window_size) {
			my ($sum, $n);
			for (my $j = 1; $j <= $window_size; ++$j) {
				my $pos_to_count = $i - $j;
        		next unless exists $coverage{$blastOut}{$pos_to_count};
		        $sum += $coverage{$blastOut}{$pos_to_count};
        		++$n;
			}
			if ($n == 0) {
				$window_average{$blastOut}{$i} = 0;
				next;
			}
			$window_average{$blastOut}{$i} = $sum / $n;
			$maxCoverage = $sum / $n unless $maxCoverage >= $sum / $n;
		}

	}
}

sub plot {
	print STDERR "\nWriting values to temporary file for plotting...";

	#each of multiple blast outputs gets a seperate plotfile
	my $j = 0;
	foreach my $blastOut (@blast_output) {
		die ("ERROR - could not create temporary plotfile $output_prefix-$j-plotfile.tmp\n") unless open(PLOT, ">$output_prefix-$j-plotfile.tmp");
		print PLOT "\"pos\",\"value\"";

		for (my $i = 1; $i <= $reference_genome_length; $i += $window_size) { 
			if ($logY && ! $window_average{$blastOut}{$i} == 0) {
				my $log = log($window_average{$blastOut}{$i});
				print PLOT "\n$i,$log";
			} else {
				print PLOT "\n$i,$window_average{$blastOut}{$i}";
			}
		}

		close PLOT;
		++$j;
	}

	#generate R script
	print STDERR "\nWriting R script to $output_prefix.tmp...";
	die ("ERROR - cannot create R script at $output_prefix.tmp\n") unless open(R, ">$output_prefix.tmp");
	#our $script;
	$script .= <<EOF;
library("maptools")
EOF

	#R: read in plotfiles for each blast output
	$j = 0;
	foreach my $blastOut (@blast_output) {
		$script .= <<EOF;
coverage$j = read.csv("$output_prefix-$j-plotfile.tmp", head=TRUE)
EOF
		++$j;
	}

	#R: initialise the pdf output
  $script .= <<EOF;
pdf("$output_prefix-w$window_size.pdf", width = $plot_width, height = $plot_height)
EOF

	#R: set up seperate par row for features if needed
	#unless (!$features_file) { #unless no features file has been specified
		#$script .= <<EOF;
#par(mfrow = c(2,1))
#EOF
	#}
	
	#R: set plot title and initialise plot area
	my $plot_title;
	if ($logY) {
		$plot_title = "ln(Coverage)";
	}
	else {
		$plot_title = "Coverage";
	}
	my $ylim;
	$ylim = $maxCoverage * 1.4;
	if ($logY) {
		$ylim = log($ylim);
	}
	my $yAxis;
	$script .= <<EOF;
par(xpd=TRUE)
plot(coverage0\$pos, coverage0\$value, type="n", ylab = c("$plot_title"), xlab = c(""), xlim=c(0, $reference_genome_length), ylim=c(0, $ylim)$yAxis)
EOF

	#R: draw a polygon representing coverage for each blast output
	my @rcolours = qw(#104E8B70 #B2222270 #228B2270 #8B0A5070 #CDAD0070);
	my $j = 0;
	foreach my $blastOut (@blast_output) {
	    $script .= <<EOF;
polygon_coords = rbind(coverage$j, c($reference_genome_length, 0), c(0,0))
polygon(polygon_coords, col=c("@rcolours[$j]"), lty=0)
EOF
		++$j;
	}

  #R: add a legend to the coverage plot my $legendText;
	
	my ($legendText, $legendCols);
	my $j = 0;
	foreach my $blastOut (@blast_output) {
		$legendText .= "\"$blastOut\",";
		$legendCols .= "\"@rcolours[$j]\",";
		++$j;
	}
	chop $legendText;
	chop $legendCols;
	$script .= <<EOF;
legendtext = c($legendText)
legendcols = c($legendCols)
legend(c("topright"), legend=legendtext, fill=legendcols)
EOF

	#R: draw on the features
	if (!$features_file) {
		$script .= <<EOF;
dev.off()
EOF
	} else {

		#$script .= <<EOF;
#plot(1, type="n", axes=F, xlim = c(0, $reference_genome_length), ylim = c(0, 10), ylab = "", xlab = "")
#EOF

		print STDERR "\nAdding features...";
		die ("ERROR - could not open features in $features_file\n") unless open(FEATURES, "<$features_file");
			our @xs;
			our @ys;
			our @labels;
			while (my $line = <FEATURES>) {
				chomp $line;
				my @fields = split(/,/, $line);
				die ("ERROR - malformed features line on line $. of $features_file - all lines should have four fields - if you want to leave a field empty, you still have to put in commas to show it's there\n") unless @fields == 4;
				print STDERR "\nAdding feature @fields[3]...";
				my ($x, $y, $l)= &shape(@fields);
				@xs= @$x;
				@ys=@$y;
				@labels=@$l;
			}
		close FEATURES;

		my $xs = join(",", @xs);
		my $ys = join(",", @ys);
		my $labels = join(",", @labels);
		$script .= <<EOF;
pointLabel(c($xs), c($ys), pos=1, labels=c($labels), xpd=NA, allowSmallOverlap = FALSE)
dev.off()
EOF

	}
	print R $script;
	close R;

	#execute R script
	print STDERR "\nExecuting R script...";
	system("R --no-save < $output_prefix.tmp");

}

sub shape {

	my @fields = @_;
	my (@xs, @ys, @labels);
	my $start_pos = @fields[0];
	my $end_pos = @fields[1];
	my $gene_label = @fields[3];
	
	my $text_pos = ($start_pos - $end_pos) * 0.5 + $end_pos;
	my $arrow_colour;
	if (@fields[2] eq "") {
		$arrow_colour = "firebrick";
	} else {
		$arrow_colour = @fields[2];
	}

	push(@xs, $text_pos);
	push(@ys, 0);
	push(@labels, "\"$gene_label\"");

	$script .= <<EOF;
draw_me = rbind(
c($end_pos,0.7),
c($start_pos, 0.7),
c($start_pos, 0.3),
c($end_pos, 0.3),
c($end_pos,0.7)
)

polygon(draw_me, col=c("$arrow_colour"), lty=1, xpd=NA)
EOF

	return (\@xs, \@ys, \@labels);
}
