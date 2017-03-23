#! /usr/bin/perl

=head1 DESCRIPTION

	interleave.pl -- Take the fwd and rev sequence files and arrange the sequence pairs in a single file such that forward sequence is followed by it's reverse pair.

=head2 USAGE

	interleave.pl -fwd forward_file.fasta -rev reverse_file.fasta [OPTIONS]
	
=head3 Examples

	# Have 2 fasta files; interleave and print to blah_int.fasta, blah_sfwd.fasta, blah_srev.fasta
	interleave.pl -fwd forward_file.fasta -rev reverse_file.fasta -o blah
	OR
	# Have 2 fastq files; interleave and print to blah_int.fastq, blah_sfwd.fastq, blah_srev.fastq	
	interleave.pl -fastq -fwd forward_file.fastq -rev reverse_file.fastq -o blah
	OR
	# Have 2 fastq files; interleave and print as fasta to blah_int.fasta, blah_sfwd.fasta, blah_srev.fasta	
	interleave.pl -fastq -fwd forward_file.fastq -rev reverse_file.fastq -o blah -outfmt fasta

=head3 Options

	-fastq	:	Use this flag if your forward and reverse sequence files are fastq
	-outfmt	:	Use this flag if you want the output of your fastq to be in fasta
	-old	:	Use this flag if your sequence descriptors follow the old illumina naming convention; i.e. if the descriptor ends with a "/1" or "/2".
	-p or -o:	Pefix for your output files; default: process_id
	-h	:	This documentation.
	-v	:	Script version.

=head1 Suggestions/Feedback/Beer

	Sunit Jain, May 2012
	sunitj [AT] umich [DOT] edu	

=cut

#######################
## MODULES
#######################
use strict;
use Getopt::Long;
use IO::Handle;

#######################
## PARAMETERS
#######################
my ($fwd, $rev, $prefix, $fastq, $old, $outfmt);
my $version="interleave.pl v0.3.7";

GetOptions(
	"o|out|p|prefix=s"	=>	\$prefix,
	"fastq"=>	\$fastq,
	"fwd|forward=s"	=>	\$fwd,
	"rev|reverse=s"	=>	\$rev,
	"outfmt:s"	=>	\$outfmt,
	"old!" => \$old,
	"h|help"	=>	sub {system('perldoc', $0); exit;},
	"v|version"=>sub{print $version."\n"; exit;},
);

#######################
## CHECKS
#######################
die "FastQ/Fasta forward file required" if (! $fwd);
die "FastQ/Fasta reverse file required" if (! $rev);
my ($out, $fwdS, $revS);
my $seqType= $fastq ? "fastq" : "fasta";
if (! $outfmt){$outfmt = $seqType};
if (! $prefix){
	$out=$$."int.".$outfmt;
	$fwdS=$$.".sfwd.".$outfmt;
	$revS=$$.".srev.".$outfmt;
}
else{
	$out=$prefix."_int.".$outfmt;
	$fwdS=$prefix."_sfwd.".$outfmt;
	$revS=$prefix."_srev.".$outfmt;
}

$/ = $fastq ? "\n" : ">";
my $beginsWith= $fastq ? "@" : ">";
my $delim= $old ? "\/" : " ";

my (%forward);
my ($prefixFlag);
my $FH=IO::Handle->new();
#######################
## MAIN
#######################
open($FH, $fwd)|| die "[ERROR1: $0] $!\n";
while(my $line=<$FH>){
	chomp $line;
	next unless $line;
	$prefixFlag=0;
	
	my ($name, $seq, $qName, $qual)= $fastq ? &parseFastq($line, $FH) : &parseFasta($line);
	my (@headerParts)= &splitHeader($name);

	my $common=shift @headerParts;
	my $strand=shift @headerParts;
#	my $prefix=join("_", @headerParts) if ($prefixFlag==1);


#	$forward{$common}{"PREFIX"}=$prefix  if ($prefixFlag==1);
	$forward{$common}{"NAME"}=$name;
	$forward{$common}{"STRAND"}=$strand;
	$forward{$common}{"SEQ"}=$seq;
	$forward{$common}{"QUAL"}=$qual if ($fastq);
}
close $FH;

my $FH=IO::Handle->new();
open($FH, $rev)|| die "[ERROR2: $0] $!\n";
open(INT, ">".$out)|| die "[ERROR3: $0] $!\n";
open(RS, ">".$revS)|| die "[ERROR4: $0] $!\n";

while(my $line=<$FH>){
	chomp $line;
	next unless $line;

	my ($name, $seq, $qName, $qual)= $fastq ? &parseFastq($line, $FH) : &parseFasta($line);
	my (@headerParts)= &splitHeader($name);
	
	my $common=shift @headerParts;
	my $strand=shift @headerParts;

	if ($forward{$common}){ # pair found!
		print INT ($outfmt eq "fastq" ? "@" : ">").$forward{$common}{"NAME"}."\n".$forward{$common}{"SEQ"}."\n"; # forward
		print INT "\+\n".$forward{$common}{"QUAL"}."\n" if ($outfmt eq "fastq"); 
		print INT ($outfmt eq "fastq" ? "@" : ">").$name."\n".$seq."\n";  #reverse
		print INT "\+\n".$qual."\n"  if ($outfmt eq "fastq"); 
		delete $forward{$common};
	}
	else{ # singletons reverse
		print RS ($outfmt eq "fastq" ? "@" : ">").$name."\n".$seq."\n";
		print RS "\+\n".$qual."\n"  if ($outfmt eq "fastq"); 
	}
}
close $FH;
close INT;
close RS;
$/="\n";


open(FS, ">".$fwdS)|| die "[ERROR5: $0] $!\n";
foreach my $n(keys %forward){ # singletons forward
#	if ($prefixFlag==1){
		print FS ($outfmt eq "fastq" ? "@" : ">").$forward{$n}{"NAME"}."\n".$forward{$n}{"SEQ"}."\n";
		print FS "\+\n".$forward{$n}{"QUAL"}."\n" if ($outfmt eq "fastq");
#	}
#	else{
#		print FS $beginsWith.$n.$delim.$forward{$n}{"STRAND"}."\n".$forward{$n}{"SEQ"}."\n";
#	}
#	print FS "\+\n".$forward{$n}{"QUAL"}."\n" if $fastq;
}
close FS;

#######################
## SUB-ROUTINES
#######################
sub splitHeader{
	my $header=shift;
	$header=~ s/^[\@\>]//;
	my (@headerParts)=split(/ /, $header);
	if($header=~ /\//){
		@headerParts=split(/\//, $header);
		return(@headerParts);
	}
#	elsif($headerParts[0]=~ /\_/){
#		my @firstPart=split(/\_/, $headerParts[0]);
#		my @secondPart=split(/\_/, $headerParts[1]);
#		$prefixFlag=1;
#		return($firstPart[-1], $secondPart[0]);
#	}
	else{
		return(@headerParts);
	}
}


sub parseFastq{
	my $line= shift;
	my $fh=shift;
	if ($line=~ /^@\w+/){
		$line=~ s/^@//;
		$line=~ s/\r//g;
		chomp $line;
		my $seqDesc=$line;

		$line=<$fh>;
		$line=~ s/\r//g;
		chomp $line;
		my $seq=$line;

		$line=<$fh>;
		$line=~ s/\r//g;
		chomp $line;
		my $qHead=$line;

		$line=<$fh>;
		$line=~ s/\r//g;
		chomp $line;
		my $qual=$line;

		die "[ERROR: line $.] Script Borked! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n" if (length($seq)!= length($qual));
		
		return ($seqDesc, $seq, $qHead, $qual);
	}
	else{ die "[ERROR: line $.] Script Borked! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n"; }

}

sub parseFasta{
	my 	$line=shift;
	my($seqDesc,@sequence)=split(/\n/, $line);
	my $seq = join ("", @sequence);
	$seq=~ s/\r//g;
	return($seqDesc, $seq);
}
