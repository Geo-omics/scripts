#!/usr/bin/perl

=head1 DESCRIPTION

	changeClasses.pl -- Edit the *.cls file to change the class numbers for a given list of contig names.

=head1 USAGE

	perl changeClasses.pl -cls esom.cls -names esom.names -list contigs.list -o output.cls

=head2 Options

	-cls	<CHARACTERS>	*.cls file produced by the esomWrapper.pl script
	-names	<CHARACTERS>	*.names file produced by the esomWrapper.pl script
	-list	<CHARACTERS>	list of contig names(replace all special chracters with '_')
	-tag	<CHARACTERS>	new class name/number to be assigned.[default=next available class number]
	-out	<CHARACTERS>	new class file. default: esom_edited.cls
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Fri Dec  6 08:59:27 EST 2013)
	sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;

my($cls, $names, $list, $tag, $out);
my $help;
my $version="changeClasses.pl\tv0.1.1";
GetOptions(
	'c|cls:s'=>\$cls,
	'n|names:s'=>\$names,
	'l|list:s'=>\$list,
	'tag:s'=>\$tag,
	'o|out:s'=>\$out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

die "Missing required input files. See '$0 -h' for help on how to use the script\n" if((! $cls)||(! $names)||(! $list));
$out=fileparse($cls,".cls")."_edited.cls" if (! $out);

open(LIST, "<".$list)||die $!;
my %index;
while(my $line=<LIST>){
	$line=strip($line);
	$index{lc($line)}++;
}
close LIST;

open(NAMES, "<".$names)||die $!;
my %num_name;
while(my $line=<NAMES>){
	$line=strip($line);
	next if ($line=~ /^%/);
	
	my($number, $window, $name)=split(/\t/, $line);
	next unless ($index{lc($name)});
	$num_name{$number}++;
}
close NAMES;

open(CLS, "<".$cls)|| die $!;
open(NCLS, ">".$out)|| die $!;
my $num_classes=-1;
while(my $line=<CLS>){
	$line=strip($line);
	if ($line=~ /^%/){
		print NCLS $line."\n";
		$num_classes++;
	}
	else{
		if($.==($num_classes+2)){
			print NCLS "%".($num_classes)."\t".
			($tag ? $tag : $num_classes)."\t".
			getRandomColor()."\n";
		}
		my($number, $class)=split(/\t/, $line);
		if($num_name{$number}){
			print NCLS $number."\t".$num_classes."\n";
		}
		else{
			print NCLS $line."\n";
		}
	}
}
close CLS;
close NCLS;

sub strip{
	my $data=shift;
	chomp $data;
	$data=~ m/^\s+/;
	$data=~ m/\s+$/;
	return $data;
}

sub getRandomColor {
    my ($r, $g, $b) = map { int rand 256 } 1 .. 3;
    my $color= join("\t", $r, $g, $b);
    return ($color);
}

