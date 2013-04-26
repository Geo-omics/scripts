#!/usr/bin/perll

=head1 Motivation

Perform vlookups on large tab-delimited files

=head1 Usage


=head2 This may get a little confusing so read carefully!

	perl vlookup.pl [-base: File you wish to add stuff 'TO'] [-from: File you wish to add stuff 'FROM']
	[-baseCol: which column in the 'BASE' file do you want matched] [-c: which column for the 'FROM' file do you want matched]
	[-cols: which column for the 'FROM' file do you want matched; adds more than one columns to the 'BASE' file]

=head1 Author

Sunit Jain, Sept 2011.
sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;

my $getDataFromThisFile;
my $matchToThisFile;
my $baseCol=1;
my $matchCol=1;
my $out=$$.".out";
my $cols;
my $col;

GetOptions(
	'base:s'=>\$matchToThisFile,
	'baseCol:i'=>\$baseCol,
	'from:s'=>\$getDataFromThisFile,
	'cols:s'=>\$cols,
	'c|col:i'=>\$col,
	'o:s'=>\$out,
	'h|help'=> sub{system('perldoc', $0); exit;}
);

open (APP, $getDataFromThisFile) || die "[ERROR] $getDataFromThisFile: $!\n";
my %index;
my $colm=&parseColInput;
$baseCol--;
$matchCol--;
print "printing column $col\n";

open (META, $matchToThisFile) || die "[ERROR] $matchToThisFile: $!\n";
while (my $line=<META>){
	next if ($line=~ m/^#/);
	chomp $line;
	$line=~ s/\r//;
	next unless $line;

	my (@stuff)=split(/\t/, $line);
	$index{$stuff[$baseCol]}=$line;
}
close META;

open (OUT, ">".$out);
while (my $line=<APP>){
	next if ($line=~ m/^#/);
	chomp $line;
	$line=~ s/\r//;
	next unless $line;
	
	my (@content)=split(/\t/, $line);
	
	if ($index{$content[$matchCol]}){
		my $matched;
		if ($cols){
			my @tmpArray=split(/\,/,$colm);
			$matched=join("\t", @content[@tmpArray]);
		}
		else{
			$matched=join("\t", @content[$colm]);
		}
		print OUT $index{$content[$matchCol]}."\t".$matched."\n";
	}
	else{
		next;
	}
}
close APP;
close OUT;


sub parseColInput{
	my $column;
	int($cols);
	if ($cols){
		chomp $cols;
		my (@tempArray1,@tempArray);
		if($cols=~ /\-/i){
			my($a1, $b1)=split(/\-/, $cols);
			@tempArray1=($a1..$b1);
		}
		foreach my $c(@tempArray1){
			$c--;
			push( @tempArray, $c);
		}
		$column=join(",", @tempArray);
	}
	elsif($col){
		$col--;
		$column=$col;
	}
	elsif(! $cols && ! $col){
		print STDERR "[Error] Missing Coloumn Parameter\n";
		exit;
	}
	return $column;
}
