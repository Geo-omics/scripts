#!/usr/bin/perl

=head1 DESCRIPTION

	itemize.pl -- You you used mapper.pl for multiple datasets and you wish to get a comparison for each dataset. Use this script.

=head1 USAGE

	perl itemize.pl -info contigs.info -mlog mapper_output.log -out output.tsv

=head2 Options

	-info	-i	<CHAR>	A file containing the dataset info for each contig. One entry per line; FORMAT: Contig_Name <TAB> Dataset_Name.
	-log	-l	<CHAR>	".log" output produced by the mapper.pl script.
	-out	-o	<CHAR>	output produced by this script

	-map	-m	<CHAR>	mapped output produced by the mapper.pl script. If you want the output from this script added to the mapper output.

	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Thu Mar 13 11:51:59 EDT 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;

my($info, $log, $map, $out);

my $help;
my $version="itemize.pl\tv0.1.0";
GetOptions(
	'i|info:s'=>\$info,
	'l|log:s'=>\$log,
	'm|map:s'=>\$map,
	'o|out:s'=>\$out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my %index;
open(INFO, "<".$info)|| die $!;
while(my $line=<INFO>){
	next if ($line=~ /^#/);
	chomp $line;
	next unless $line;

	my ($name, $set)=split(/\t/, $line);
	$index{lc($name)}=$set;
}
close INFO;

my (%distribution, %sets);
open(LOG, "<".$log)|| die $!;
while(my $line=<LOG>){
	next if ($line=~ /^#/);
	chomp $line;
	next unless $line;

	my($ref, @qList)=split(/\t/, $line);
	foreach my $query(@qList){
		my $set=$index{lc($query)} ? $index{lc($query)} : "Unassigned";
		$distribution{$ref}{$set}++;
		$sets{$set}++;
	}
}
close LOG;
my @setList=keys(%sets);

open(OUT, ">".$out) || die $!;
if(-e $map){
	open(MAP, "<".$map) || die $!;
	while(my $line=<MAP>){
		if(($. <= 5) && ($line=~ /^#/)){
			print OUT $line;
		}
		elsif(($. == 6) && ($line=~ /^#/)){
			chomp $line;
			$line.= "\t".$_ foreach @setList;
			print OUT $line."\n";
		}
		else{
			chomp $line;
			my($ref, @data)=split(/\t/, $line);
			foreach my $set(@setList){
				my $num=$distribution{$ref}{$set} ?  $distribution{$ref}{$set} : "0";
				$line.= "\t". $num;
			}
			print OUT $line."\n";
		}
	}
	close MAP;
}
else{
	print OUT "# REFERENCE_REGION";
	print OUT "\t".$_ foreach (@setList);
	print OUT "\n";
	foreach my $ref(keys %distribution){
		my $line=$ref;
		foreach my $set(@setList){
			my $num=$distribution{$ref}{$set} ?  $distribution{$ref}{$set} : "0";
			$line.="\t".$num;
		}
		print OUT $line."\n";
	}
}
close OUT;

print "# Total Regions identified in ".$_."\t".$sets{$_}."\n" foreach (keys %sets);
exit 0;
