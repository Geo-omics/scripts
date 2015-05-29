#!/usr/bin/perl

=head1 USAGE

	Get me the .contig file from minimus and I'll tell you which contigs combined to form the minimus output.

	perl minimusMetadata.pl -c < .contig file from minimus> -o <output_File_Name>

=head2 OPTIONAL

	-len	minimum length
	-list	List of minimus scaffolds you wish to extract the contig names for
	-summary Don't get the read names, just the numbers. Also try the following command for the same effect
	#	system("grep "^##" filename.contig | cut -f 1,2 -d " " | sort -g -r -k 2 > summary_afg.map");

=head1 AUTHOR
	
	Sunit Jain, December 2012
	
=cut

use strict;
use Getopt::Long;

my $contigFile;
my $out=$$.".map";
my $list;
my $summary;
my $minLength=0;
GetOptions(
	'c:s'=>\$contigFile,
	'o:s'=>\$out,
	'l|list:s'=>\$list,
	'len|length:i'=>\$minLength,
	'summary'=>\$summary,
	'h'=>sub{system('perldoc', $0); exit;},
);

my @listOfHeaders=`grep '^#' $contigFile`;



##1 4 87965 bases, 00000000 checksum.
#Sequence0000000001(0) [] 87965 bases, 00000000 checksum. {1 87965} <1 87965>
#Sequence0000027928(12770) [] 403 bases, 00000000 checksum. {1 403} <12771 13173>
#Sequence0000029868(42929) [RC] 368 bases, 00000000 checksum. {368 1} <42930 43297>
#Sequence0000029936(82824) [RC] 366 bases, 00000000 checksum. {366 1} <82825 83190>

open(OUT, ">".$out) || die $!;
print OUT "#Scaffold\t# of Reads\tContigs(Position on Scaffold); * = Reverse complement of this sequence was used\n";
for(my $i=0; $i <= $#listOfHeaders; $i++){
	my $header=$listOfHeaders[$i];
	chomp $header;

	my $line;
	if($header=~ /^##/){
		$header=~ s/\#//g;
		my ($scaffold, $members, $length, @etc)=split(/\s+/, $header);
		next unless ($length >= $minLength);
		$line=$scaffold."\t".$members."\t";
		unless($summary){
		for(my $j=0; $j < $members; $j++){
			$i++;
			my ($contig, @etc)=split(/\s+/, $listOfHeaders[$i]);
			$contig=~ s/\#//g;
			$contig.="*" if ($etc[0]=~ /\[RC\]/);
			$line.=$contig."\t";
		}
		}
		$line=~ s/\t$/\n/;
		print OUT $line;
	}
}
close OUT;
