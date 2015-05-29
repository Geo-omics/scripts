#!/usr/bin/perl

=head1 DESCRIPTION

top5.pl -- This program gets the top 5 hits for each query in a blast tabular output file.

=head1 USAGE

	perl top5.pl -b <blast output> -t <top X hits; default: 5> -o <output_File_Name>
	
=head2 Example

perl top5.pl -b sample.blastn -o sample_topHits.blastn -t 1

=head2 Options

-b	<CHAR>	Blast output in tabular (-outfmt 6/7 OR -m 8/9). [REQUIRED]
-t	<INT>	# of top hits. [default = 5 ]
-o	<CHAR>	output file; same format as blast file
-no_self	<BOOLEAN>	remove self hits. helpful in a self blast.
	
-version -v	<BOOLEAN>	version of the current script
-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Tue May 13 16:07:06 EDT 2011)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;

my $in;
my $x=5;
my $out=$$.".top5";
my $no_self;
my $version="top5.pl\tv0.1.1";
GetOptions(
	'b:s'=>\$in,
	't:i'=>\$x,
	'o|out:s'=>\$out,
	'no_self'=>\$no_self,
	'h'=> sub{system("perldoc $0 \| cat"); exit;},
);

open (TB, $in) || die "[err] $in: $! \n";
my %file;
while (my $line=<TB>){
	next if ($line=~ m/^#/);
	chomp $line;
	$line=~ s/\r//;
	next unless $line;
	
	my($query, @etc)=split(/\t/, $line);
	if($no_self){
		my ($subj, @desc)=split(/ /, $etc[0]);
		next if $query eq $subj;
	}
	
	push(@{$file{$query}}, $line);
}
close TB;

my($nothing, $fileName)=split(/\_/, $in);
open (OUT, ">".$out)|| die "[err] $out: $! \n";;
while(my($k,$val)= each(%file)){
	my $count=0;
	foreach my $v(@{$val}){
		if ($count < $x){
			print OUT $v."\n";
		}
		$count++;
	}
}
close OUT;
