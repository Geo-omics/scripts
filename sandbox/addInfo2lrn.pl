#!/usr/bin/perl

=head1 DESCRIPTION

	addInfo2lrn.pl -- Adds additional columns of data to the LRN files. You may want to add data like bin coverage, GC etc.

=head1 USAGE

	perl addInfo2lrn.pl -info scaffolds.info -lrn esom.lrn -names esom.names -out file_with_additional_data.lrn

=head2 Options

	-info	<CHAR>	file with additional data; FORMAT: Scaffold_Name<TAB>Feature1<TAB>Feature2<TAB>...<TAB>FeatureN
	-lrn	<CHAR>	LRN file created from the ESOM wrapper script.
	-names	<CHAR>	NAMES file created from the ESOM wrapper script.
	-out	<CHAR>	output file in the *.lrn format with the additional columns.

	-map	<CHAR>	If your contig/scaffold names were changed at some point such that the coverage names don;t match your scaffold names anymore, create a file in the format:
			scaffold_name_asIn_Coverage <TAB> scaffold_name_asIn_ESOM_input_fasta

	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Mon Nov 18 11:34:44 EST 2013)
	sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;

my $help;
my $version="addInfo2lrn.pl\tv0.0.4";
my ($info, $lrn, $names, $out, $map);
GetOptions(
	'map:s'=>\$map,
	'info:s'=>\$info,
	'lrn:s'=>\$lrn,
	'names:s'=>\$names,
	'out|o:s'=>\$out,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my %alias;
if($map){
	open(MAP, "<".$map)|| die $!;
	while(my $line=<MAP>){
		next if ($line=~ /#/);
		chomp $line;
		next unless $line;

		my($original, $JGI)=split(/\t/, $line);

		$alias{uc($original)}=$JGI;
	}
	close MAP;
}

open(INFO, "<".$info) || die $!;
my (%INFO,@newHeaders);
while(my $line=<INFO>){
	chomp $line;
	next unless $line;
	my ($name, @data)=split(/\t/, $line);

	if (($.==1) && ($line=~ /^#/)){
		push(@newHeaders, @data);
		next;
	}
	elsif($.==1){
		for(my $i=1; $i <= scalar(@data); $i++){
			push(@newHeaders, "Feature".$i);
		}
	}
	
	$name=$alias{uc($name)} if ($map);

	$INFO{uc($name)}=join("\t", @data);
}
close INFO;

open(NAMES, "<".$names) || die $!;
my %NAMES;
while(my $line=<NAMES>){
	chomp $line;
	next unless $line;
	next if ($line=~ /^#/);
	next if ($line=~ /^%/);
	
	my ($num, $ann, $contig)=split(/\t/, $line);
	
	if ($INFO{uc($contig)}){
		$NAMES{uc($num)}=$INFO{uc($contig)};
	}
}
close NAMES;
undef %INFO;

open(LRN, "<".$lrn)|| die $!;
open(OUT, ">".$out)|| die $!;
while(my $line=<LRN>){
	chomp $line;
	next unless $line;
	next if ($line=~ /^#/);
	
	my $addCols=scalar(@newHeaders);
	if ($line=~ /^\%/){
		if($.==1){
			print OUT $line."\n";
		}
		if($.==2){
			$line=~ s/\%//;
			print OUT "%".($line+$addCols)."\n";
		}
		if($.==3){
			for(my $i=0; $i < $addCols; $i++){
				$line.="\t1";
			}
			print OUT $line."\n";
		}
		if($.==4){
			for(my $i=0; $i < $addCols; $i++){
				$line.="\t".$newHeaders[$i];
			}
			print OUT $line."\n";
		}
	}
	else{
		my($num, @data)=split(/\t/, $line);
		if($NAMES{$num}){
			$line.="\t".$NAMES{$num};
		}
		else{
			$line.="\t".(join("\t",((0) x $addCols)));
		}
		print OUT $line."\n";
	}
}
close LRN;
close OUT;

exit 0;
