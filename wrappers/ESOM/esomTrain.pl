#!/usr/bin/perl

=head1 Description

	esomTrain.pl - Given tetramer frequency script output files, the script will normalize the files and train the ESOM.

=head2 Dependencies

	ESOM 1.1 or above

=head1 Usage

	perl esomTrain.pl -lrn esom.lrn -cls esom.cls -rows #1 -cols #2
	OR
	perl esomTrain.pl -lrn esom.lrn -cls esom.cls -rows #1 -cols #2 -info file.info

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
				'Scale' - to scale the dataset from 0-1; [default]
				'ZT'	- to Z-transform the dataset;

=head3 Advance Options

	-k	[float]		k for k-Batch as percentage of data set.[default is 15% written as = 0.15]
				Percentage of all data vectors that is processed between map updates.
	-algorithm	[characters]	training algorithm; possible choices:
					'kbatch' - k-batch training [default]
					'online' - standard online training
	-startRadius	[integer]	start value for radius, default: half the smaller length of the grid.
	-epochs		[integer]	The number of iterations over the training data; default=20
	
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message. press "q" to exit this screen.

=head1 Author

	Sunit Jain, (Thu Aug 15 10:56:11 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;
use File::Spec;
use File::Basename;

# Required
my ($lrn, $COLS, $ROWS, $CLS);

# Optional
my $info;
my $names;
my $bmSearch="standard";
my $ALGO="kbatch";
my $K=0.15;
my $epochs=20;
my $dist="euc";
my ($startRadius,$endRadius,$norm,$info);
my $version="esomTrain.pl\tv0.0.9";
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
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);
print "\# $version\n";

### Sanity Checks ###
$norm=lc($norm);
$ALGO=lc($ALGO);

if (! $lrn ||  ! $COLS || ! $ROWS || ! $CLS){
	die "Missing required files. See `-h' for help on how to use the script\n"
}

if (! $startRadius){ 
	my @compare=sort{$a <=> $b} ($ROWS, $COLS);
	$startRadius= int(($compare[0] / 2) + 0.5);
}

my($SCALE, $ZT, $LRN, $OUT, $bmOut);
my ($prefix,$path,$suffix) = fileparse($lrn,"lrn");
if ($norm eq "zt"){
	$ZT++;
	$LRN=$prefix."zt_norm.lrn";
	$OUT=$prefix."zt_norm_".$ROWS."x".$COLS."e".$epochs.".wts";
	$bmOut=$prefix."zt_norm_".$ROWS."x".$COLS."e".$epochs.".bm";
}
elsif($norm eq "scale"){
	$SCALE++;
	$LRN=$prefix."01_norm.lrn";
	$OUT=$prefix."01_norm_".$ROWS."x".$COLS."e".$epochs.".wts";	
	$bmOut=$prefix."01_norm_".$ROWS."x".$COLS."e".$epochs.".bm";
}
else{
	$SCALE++;
	$LRN=$prefix."01_norm.lrn";
	$OUT=$prefix."01_norm_".$ROWS."x".$COLS."e".$epochs.".wts";	
	$bmOut=$prefix."01_norm_".$ROWS."x".$COLS."e".$epochs.".bm";
	warn "[WARNING:]\tA normalization scheme was not mentioned, applying default normalization (scaling column values to [0,1]).\n";
}


### MAIN ###

my %lrnHeaders;
my %matrix;
my (@key,@infoKey);
my %stats;
my %colStats;
my %NAMES;
if(-e $info){
	## Add ability to concatenate more columns before normalizing ##
	die "names file required with the info flag, use -names flag\n" if (! $names);
	open (NAMES, "<".$names)|| die "ERROR with names file: $names\t$!\n";
	while(my $line=<NAMES>){
		next if $line=~ m/^\%/;
		chomp $line;
		next unless $line;
		
		my($num, $annotation, $contig)=split(/\t/, $line);
		push(@{$NAMES{$contig}}, $num);
	}
	open(INFO, "<".$info)|| die "ERROR with info file: $info\t$!\n";
	while(my $line=<INFO>){
		chomp $line;
		if ($line=~ /^\#/){
			my($name, @moreHeaders)=split(/\t/, $line);
			push(@key,@moreHeaders); #Ordered list of additional features
		}
		elsif($.==1){
			my (@cols)=split(/\t/, $line);
			my $i=0;
			my $header;
			while ($i < $#cols){
				$i++;
				push(@key,"Feature".$i);
			}
			
			&for_all_parts_of_contig($line, \@key);
		}
		else{
			&for_all_parts_of_contig($line, \@key);
		}
	}
	close INFO;
	@infoKey=@key;
}

sub for_all_parts_of_contig{
	my $line=shift;
	my $keyArr=shift;
	my($name, @content)	=split(/\t/, $line);
	foreach my $part(@{$NAMES{$name}}){
		$line=~ s/^($name)/$part/;
		&getStats($line); # the value of the array @key changes when this function is called, happens only here.
	}
	return;
}

my @key=();
open(LRN, "<".$lrn)|| die "ERROR with lrn file: $LRN\t$!\n";
open(NORM, ">".$LRN)|| die "ERROR with lrn file: $lrn\t$!\n";
while(my $line=<LRN>){
	next if $line=~ /^#/;
	chomp $line;
	next unless $line;
	if ($line=~ /^\%/){
		if($.==1){
			print NORM $line."\n";
		}
		if($.==2){
			$line=~ s/\%//;
			my $numHeaders=$line+scalar(@infoKey);
			print NORM "\% $numHeaders\n";
		}
		if($.==3){
			for(my $i=0; $i<= $#infoKey; $i++){
				$line.="\t1";
			}
			print NORM $line."\n";
		}
		if($.==4){
			my($name, @kmers)=split(/\t/, $line);
			push(@key,@kmers); # Ordered list of KMer Headers.
			for(my $i=0; $i<= $#infoKey; $i++){
				$line.="\t$infoKey[$i]";
				push(@key,$infoKey[$i]);
			}
			print NORM $line."\n";
		}
	}
	else{
		&getStats($line);
	}
}
close LRN;

#push(@key,@infoKey);
## Average and Std Dev ##
if($ZT){
	foreach my $key(keys %colStats){
		$stats{"Avg"}{$key}= &average(\@{$colStats{$key}});
		$stats{"Stdev"}{$key}= &stdev(\@{$colStats{$key}});
	}
}
undef %colStats;

my @sortedContigList= sort {$a <=> $b} keys %matrix;
foreach my $contig(@sortedContigList){
	my $line.=$contig."\t";
	foreach my $kmer(@key){
		my $normalized;
		if ($SCALE){
			$normalized= &scale01($contig, $kmer);
		}
		elsif($ZT){
			$normalized= &zt($contig, $kmer);
		}
		$line.=$normalized."\t";
	}
	$line=~ s/\t$/\n/;
	print NORM $line;
}
close NORM;

exit;
### Let the training begin... ###
my $train="esomtrn --permute --out $OUT -b $bmOut --cls $CLS --lrn $LRN --algorithm $ALGO --rows $ROWS --columns $COLS"
	.($epochs ? " --epochs $epochs " : "")
	.($ALGO eq "kbatch" ? " -k $K " : "")
	.($bmSearch ? "--bmsearch $bmSearch " : "")
	.($startRadius ? "--start-radius $startRadius " : "")
	.($endRadius ? "--end-radius $endRadius " : "")
	.($dist ? "--dist $dist " : "");

print "Executing...\n".$train."\n";
system($train);

### Sub Routines ##

sub getStats{
	my $line=shift;
	my($name,@data)=split(/\t/, $line);
	for(my $i=0; $i<= $#key; $i++){
		$matrix{$name}{$key[$i]}=$data[$i];
## Maximum ##
		if($stats{"Max"}{$key[$i]}){
			if($data[$i] > $stats{"Max"}{$key[$i]}){
				$stats{"Max"}{$key[$i]}=$data[$i];
			}
		}
		else{
			$stats{"Max"}{$key[$i]}=$data[$i];
		}
## Minimum ##
		if($stats{"Min"}{$key[$i]}){
			if($data[$i] < $stats{"Min"}{$key[$i]}){
				$stats{"Min"}{$key[$i]}=$data[$i];
			}
		}
		else{
			$stats{"Min"}{$key[$i]}=$data[$i];
		}
## For Average and Std Dev ##
		push(@{$colStats{$key[$i]}}, $data[$i]);
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

	my $zscore;
	unless($stdev == 0){
		$zscore=(($value-$avg)/$stdev);
	}
	else{
		$zscore=0;
	}
	return $zscore;
}

sub average{
	my($data) = @_;
	if (not @$data) {
		die("Empty array\n");
	}
	my $total = 0;
	foreach (@$data) {
		$total += $_;
	}
	my $average = $total / @$data;
	return $average;
}

sub stdev{
	my($data) = @_;
	if(@$data == 1){
		return 0;
	}
	my $average = &average($data);
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
}

