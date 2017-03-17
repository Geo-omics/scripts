#!/usr/bin/perl

=head1 Description

	esomTrain.pl - Given tetramer frequency script output files, the script will normalize the files and train the ESOM.

=head2 Dependencies

	ESOM 1.1 or above

=head1 Usage

	perl esomTrain.pl -lrn esom.lrn -cls esom.cls -rows #1 -cols #2
	OR
	perl esomTrain.pl -lrn esom.lrn -cls esom.cls -rows #1 -cols #2 -info file.info -names esom.names -scripts /path/to/other/scripts

=head2 Required options

	-lrn	[characters]	*.lrn file produced by the tetramer frequency script
	-cls	[characters]	*.cls file produced by the tetramer frequency script
	-rows	[integer]	number of rows in the ESOM grid; vertical size of the map
	-cols	[integer]	number of columns in the ESOM grid; horizontal size of the map

	NOTE:	The product of rows and columns, i.e. the number of neurons should be at least 1K neurons.
	The ratio of rows and columns should be significantly different from unity.

=head2 More Options

	-info	[characters]	file with information associated with the sequences; -names flag also required
				Column1=Column1 name of scaffold (as it appears in col3 of names file); Column2=%GC(example)...ColumnN=any_numerical_metric;
	-names	[characters]	*.names file produced by the tetramer frequency script; required with the -info flag
	-norm	[characters]	Normalization strategy; possible choices:
				'RZT'	- to robust-zt the dataset (replaces the mean and std dev with median and median absolute deviation); [default]
				'Scale' - to scale the dataset from 0-1;
				'ZT'	- to Z-transform the dataset;

=head3 Advance Options

	-scripts	[characters]	path to the "addInfo2lrn.pl" script. [default=current working directory]
	-k	[float]		k for k-Batch as percentage of data set.[default is 15% written as = 0.15]
				Percentage of all data vectors that is processed between map updates.
	-algorithm	[characters]	training algorithm; possible choices:
					'kbatch' - k-batch training [default]
					'online' - standard online training
	-startRadius	[integer]	start value for radius, default: half the smaller length of the grid.
	-epochs		[integer]	The number of iterations over the training data; default=20
	-cmd	[boolean]	tells the script to print the ESOM training command to the screen but not to execute it.
	-debug	[boolean]	same as '-cmd' but also print some basic stats per k-mer.
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message. press "q" to exit this screen.

=head1 Author

	Sunit Jain, (Thu Aug 15 10:56:11 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;
use File::Spec::Functions;
use File::Basename;

# Required
my ($lrn, $COLS, $ROWS, $CLS);

# Optional
my $bmSearch="standard";
my $ALGO="kbatch";
my $K=0.15;
my $epochs=20;
my $dist="euc";
my ($startRadius,$endRadius,$norm,$info,$names,$DEBUG, $CMD,$scripts);
my $version="esomTrain.pl\tv0.0.17b";
GetOptions(
	'lrn=s'=>\$lrn,
	'cls=s'=>\$CLS,
	'r|rows=i'=>\$ROWS,
	'c|cols|columns=i'=>\$COLS,
	'a|algorithm:s'=>\$ALGO,
	'k:s'=>\$K,
	'bms|bmsearch:s'=>\$bmSearch,
	'e|epochs:s'=>\$epochs,
	'rs|start-radius:i'=>\$startRadius,
	're|end-radius:i'=>\$endRadius,
	'd|dist:s'=>\$dist,
	'norm:s'=>\$norm,
	'i|info:s'=>\$info,
	'n|names:s'=>\$names,
	'scripts:s'=>\$scripts,
	'cmd'=>\$CMD,
	'debug'=>\$DEBUG,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

### Sanity Checks ###
$norm=lc($norm);
$ALGO=lc($ALGO);

if (! $lrn ||  ! $COLS || ! $ROWS || ! $CLS){
	die "[FATAL] Missing required files. See `-h' for help on how to use the script\n"
}

if ($info && ! $names){
	die "[FATAL] Using the '-info' flag requires that you provide the provide the *.names file as well, using '-names'.\n Info File='$info'\nNames File='$names' not found!\n";
}


### Global Variables ###
if (! $startRadius){ 
	my @compare=sort{$a <=> $b} ($ROWS, $COLS);
	$startRadius= int(($compare[0] / 2) + 0.5);
}

my $local_bm_radius=int(($startRadius/4) + 0.5);

my($SCALE, $ZT, $RZT, $LRN, $OUT, $bmOut);
my ($prefix,$path,$suffix) = fileparse($lrn,".lrn");
if ($info){	$prefix.="_with_Info"}

if ($norm eq "zt"){
	$ZT++;
	$LRN=$prefix.".zt_norm.lrn";
	$OUT=$prefix.".zt_norm_".$ROWS."x".$COLS."e".$epochs.".wts";
	$bmOut=$prefix.".zt_norm_".$ROWS."x".$COLS."e".$epochs.".bm";
}
elsif($norm eq "scale"){
	$SCALE++;
	$LRN=$prefix.".01_norm.lrn";
	$OUT=$prefix.".01_norm_".$ROWS."x".$COLS."e".$epochs.".wts";	
	$bmOut=$prefix.".01_norm_".$ROWS."x".$COLS."e".$epochs.".bm";
}
else{
#	$RZT++;
	$LRN=$prefix.".no_norm.lrn";
	$OUT=$prefix.".".$ROWS."x".$COLS."e".$epochs.".wts";	
	$bmOut=$prefix.".".$ROWS."x".$COLS."e".$epochs.".bm";
	warn "[WARNING]\tA normalization scheme was either not mentioned or not valid, no normalization has been applied.\n" if ((! $norm) || ($norm != "rzt"));
}


### MAIN ###
if (-e $info){
	if (! $scripts){
		$scripts="./";
	}
	my $info2lrn=catfile($scripts, "addInfo2lrn.pl");
	
	die "[FATAL] Script 'addInfo2lrn.pl' not found. Please use the '-scripts' flag to point to the location.\n" unless (-e $info2lrn);
	my $tmpLrn=$prefix.".lrn";
	system("perl $info2lrn -names $names -lrn $lrn -info $info -o $tmpLrn");
	$lrn=$tmpLrn;
}

my (%matrix, %stats, %colStats) ;
my @key;
if ( ! $CMD){
open(LRN, "<".$lrn)|| die "ERROR with reading the lrn file: $lrn\t$!\n";
open(NORM, ">".$LRN)|| die "ERROR with writing to the lrn file: $LRN\t$!\n";
while(my $line=<LRN>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;

	if($.<=3){
		print NORM $line."\n";
	}
	elsif($.==4){
		(my $header, @key)=split(/\t/, $line);
		print NORM $line."\n";
	}
	else{
		&getStats($line, \@key);
	}
}
close LRN;
}
#push(@key,@infoKey);
## Average and Std Dev ##
if($ZT){
	print "#KMER\tAverage\tMedian\tStdDev\n" if ($DEBUG);
	foreach my $col(keys %colStats){
		$stats{"Avg"}{$col}= &average(\@{$colStats{$col}});
		$stats{"Stdev"}{$col}= &stdev(\@{$colStats{$col}}, $stats{"Avg"}{$col});
		my $median=&median(\@{$colStats{$col}}) if ($DEBUG);
		print $col."\t".$stats{"Avg"}{$col}."\t".$median."\t".$stats{"Stdev"}{$col}."\n" if ($DEBUG);
	}
}
elsif($RZT){
	print "#KMER\tAverage\tMedian\tMedian_Absolute_Dev\n" if ($DEBUG);
	foreach my $col(keys %colStats){
		my $average= &average(\@{$colStats{$col}}) if ($DEBUG);
		$stats{"Avg"}{$col}=&median(\@{$colStats{$col}});
		$stats{"Stdev"}{$col}= &median_absolute_deviation(\@{$colStats{$col}},$stats{"Avg"}{$col});
		print $col."\t".$average."\t".$stats{"Avg"}{$col}."\t".$stats{"Stdev"}{$col}."\n" if ($DEBUG);
	}
}
undef %colStats;

my @sortedContigList= sort {$a <=> $b} keys %matrix;
foreach my $contig(@sortedContigList){
	last if($CMD);
	my $line.=$contig."\t";
	foreach my $kmer(@key){
		my $normalized;
		if ($SCALE){
			$normalized= &scale01($contig, $kmer);
		}
		elsif(($ZT) || ($RZT)){
			$normalized= &zt($contig, $kmer);
		}
		$line.=$normalized."\t";
	}
	$line=~ s/\t$/\n/;
	print NORM $line;
}
close NORM;

### Let the training begin... ###
my $train="esomtrn --permute --out $OUT -b $bmOut --cls $CLS --lrn $LRN --algorithm $ALGO --rows $ROWS --columns $COLS -bmc $local_bm_radius --start-radius $startRadius"
	.($epochs ? " --epochs $epochs" : "")
	.($ALGO eq "kbatch" ? " -k $K" : "")
	.($bmSearch ? " --bmsearch $bmSearch" : "")
	.($endRadius ? " --end-radius $endRadius" : "")
	.($dist ? " --dist $dist" : "");

print "#Executing...\n".$train."\n";

system($train) unless($DEBUG || $CMD);

### Sub Routines ##
sub getStats{
	my $line=shift;
	my $keyArr=shift;
	
	my @featureKey=@$keyArr;
	my($name,@data)=split(/\t/, $line);
	for(my $i=0; $i<= $#featureKey; $i++){
		$matrix{$name}{$featureKey[$i]}=$data[$i];
## Maximum ##
		if($stats{"Max"}{$featureKey[$i]}){
			if($data[$i] > $stats{"Max"}{$featureKey[$i]}){
				$stats{"Max"}{$featureKey[$i]}=$data[$i];
			}
		}
		else{
			$stats{"Max"}{$featureKey[$i]}=$data[$i];
		}
## Minimum ##
		if($stats{"Min"}{$featureKey[$i]}){
			if($data[$i] < $stats{"Min"}{$featureKey[$i]}){
				$stats{"Min"}{$featureKey[$i]}=$data[$i];
			}
		}
		else{
			$stats{"Min"}{$featureKey[$i]}=$data[$i];
		}
## For Average and Std Dev ##
		push(@{$colStats{$featureKey[$i]}}, $data[$i]);
	}
	return;
}

sub scale01{
	my $contig=shift;
	my $kmer=shift;
	
	my $value=$matrix{$contig}{$kmer};
	my $min=$stats{"Min"}{$kmer};
	my $max=$stats{"Max"}{$kmer};
	
	my $s;
	if ($max==$min){
		$s=0.5;
	}
	else{
		$s=(($value-$min)/($max-$min));
	}
	return $s;
}

sub zt{
	my $contig=shift;
	my $kmer=shift;
	
	my $value=$matrix{$contig}{$kmer};
	my $avg=$stats{"Avg"}{$kmer};
	my $stdev=$stats{"Stdev"}{$kmer};

	my $zscore=0;
	unless($stdev == 0){
		$zscore=(($value-$avg)/$stdev);
	}
	return $zscore;
}

sub average{
	my $data  = shift;
	die("Empty array in sub-routine 'average'\n") if (! @$data);
	my $total = 0;
	foreach (@$data) {
		$total += $_;
	}
	my $average = $total / scalar(@$data);
	return $average;
}

sub stdev{
	my $data = shift;
	my $average = shift; # &average($data);
	die("Empty array in sub-routine 'std'\n") if (! @$data);
	
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = sqrt($sqtotal / (scalar(@$data)));
	return $std;
}

sub median { 
	my $data = shift;
	die("Empty array in sub-routine 'median'\n") if (! $data);
	my $count = scalar(@$data);

	# Sort a COPY of the array, leaving the original untouched 
	my @array = sort { $a <=> $b } @$data;

	if (($count % 2) == 0) {
		my $mid=$count/2;
		return ($array[$mid] + $array[$mid + 1])/2;
	}
	else{
		return ($array[($count+1)/2]);
	} 
}

sub median_absolute_deviation{
	my $data = shift;
	my $MEDIAN = shift;
	die("Empty array in sub-routine 'remove_outliers'\n") if (! @$data);
	
	my($Q1, $Q3)=&get_quartiles($data, $MEDIAN);
	my $constant=1/$Q3;
#	my $constant=1.4826; # assuming the the data is normally distributed.
	my @abs_data;
	push (@abs_data, abs($_ - $MEDIAN)) foreach (@$data);
	my $abs_median=&median(\@abs_data);
	
	my $mad= $constant * $abs_median;
	
	return $mad;
}

sub get_quartiles{
	my $data = shift;
	my $median = shift;
	die("Empty array in sub-routine 'significant_data_range'\n") if (! @$data);
	
	my @array = sort { $a <=> $b } @$data;

	my (@Q1, @Q3);
	foreach my $arr(@array){
		if($arr > $median){
			push(@Q3, $arr);
		}
		elsif($arr < $median){
			push(@Q1, $arr);		
		}
	}
	my $q1 = &median(\@Q1);
	my $q3 = &median(\@Q3);

	return ($q1, $q3);
}
