#! /usr/bin/perl

=head1 Description

	getMappedReadList: give me a sorted bam(-b) file and a list of contigs(-l) and I'll tell you which reads mapped(-o) onto those list of contigs.

=head2 NOTE: For Single Contig

	If you just have one contig, use the following command instead:
		samtools view -F0x4 sortedBamFile contigName | cut -f 1 >> outputFile.list
	Some other commands to try:
		# of unmapped reads: 
		samtools view -f0x4 -c sortedBamFile contigName 
		# of uniquely mapped reads:
		samtools view -F0x4 -c -q 1 sortedBamFile contigName
	where, -c = Instead of printing the alignments, only count them and print the total number

=head2 Preparation

	[1] Make sure your sequence headers don't have any funny characters like: (,);:
	[2] Make sure you have an indexed sorted bam file in the same folder; and
	[3] The samtools module is loaded.

=head2 Options

	-b or bam:	sorted bam file; 
	-l or list:	list of contig used for the mapping.
	-o or out:	output list

=head2 Questions/Comments/Feedback/Beer

	April 2012, Sunit Jain (sunitj [AT] umich [DOT] edu)

=cut

use strict;
use Getopt::Long;

my $bam;
my $list;
my $out=$$.".list";
my $verbose;

GetOptions(
	"b|bam:s"=>\$bam,
	"l|list:s"=>\$list,
	"o|out:s"=>\$out,
	"verbose"=>\$verbose,
	"h|help"=>sub{system('perldoc', $0); exit;},
);

if (! $bam || ! $list){ print $_." Not Found!!";}

my @contigsList;
open (LIST, $list)|| die "[ERROR: $0]: $!\n";
while(my $line=<LIST>){
	next if $line=~ /^\W/;
	chomp ($line);
	$line=~ s/\r//;

	next unless $line;
	my($contig, @etc)=split(/\s+/, $line);

	$contig=~ s/(\W)/\\$1/g;
	push(@contigsList, $contig);
}
close LIST;


my $i=0;
my $tmp=$$.".tmp";
print "Aggregating a list of mapped reads..\n";
foreach my $c(@contigsList){
	$i++;
	my $cmd="samtools view -F0x4 $bam $c | cut -f 1 >> $tmp";
	print "[$i]: ".$cmd."\n" unless ($verbose);
	system($cmd);
}

print "Removing duplicate entries..\n";
system("sort -u $tmp > $out");

print "Removing temporary files..\n";
unlink $tmp;

print "Done!\n";
exit;
