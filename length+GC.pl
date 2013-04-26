#! /usr/bin/perl
=head1 Description

	This program takes a fasta file, extracts length and %GC information (if '-gc' is specified)

=head2 Usage

	length+GC.pl -f input.fasta

=head3 Options

	-gc	Calculate GC content.
	-len calculate for sequences abovea certain length only.	
    
=head1 Author

	Created by Gene W. Tyson
	modified by Sunit Jain

=cut


use strict;
use Getopt::Long;
my $calcGC;
my $fasta;
my $minLen=1;

GetOptions(
	"gc"=>\$calcGC,
	"f:s"=>\$fasta,
	"len:i"=>\$minLen,
	"h|help"=>sub{system('perldoc', $0); exit;},
);

&help if ! $fasta;

open (CONTIGS, $fasta) || die "Couldn't open $fasta\n";
$/= ">";
my (%sequences, @names);
while (my $b = <CONTIGS>) {
    chomp $b;
    next unless $b;
    my ($name, @sequence) = split (/\n/, $b);
    my $seq = join ("", @sequence);
    my $length = length($seq);
	next unless ($length >=$minLen);

	unless ($calcGC){
		print "$name\t$length\n" ;
	}
	else{
	    push (@names, $name);
    	$sequences{$name} = uc $seq;
	}
}
close CONTIGS;

exit unless ($calcGC);

foreach my $value (@names) {
	my ($g,$c, $size);
    my $seq = $sequences{$value};
    my @sequence = (split //, $seq);
    my $size = @sequence;

	foreach my $base (@sequence) {
		if ($base eq "G") {$g++;}
		if ($base eq "C") {$c++;}
	}
	my $GC = (($g+$c)/$size);
	my $GC_conversion = (int($GC*1000))/1000;
	print "$value\t$GC_conversion\t$size\n";
}

sub help{
	system('perldoc', $0);
	exit;
}

exit;


