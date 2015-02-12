#! /usr/bin/perl
=head1 Description

	This program takes a fasta file, extracts length and %GC information (if '-gc' is specified)

=head2 Usage

	length+GC.pl -f input.fasta

=head3 Options

	-gc	Calculate GC(%) content.
	-len calculate for sequences abovea certain length only.	
    
=head1 Author

	Sunit Jain

=cut


use strict;
use Getopt::Long;
use FileHandle;

my $calcGC;
my $fasta;
my $minLen=1;
my $version="0.1.1";

GetOptions(
	"gc"=>\$calcGC,
	"f:s"=>\$fasta,
	"len:i"=>\$minLen,
	"v|version"=>\$version,
	"h|help"=>sub{system('perldoc', $0); exit;},
);

&help if ! $fasta;

my $CONTIGS=FileHandle->new();
open ($CONTIGS, "<",$fasta) || die "Couldn't open $fasta\n";
$/= ">";
my (%sequences, @names);
while (my $b = <$CONTIGS>) {
    chomp $b;
    next unless $b;
    my ($name, @sequence) = split (/\n/, $b);
    my $seq = join ("", @sequence);
    my $length = length($seq);
	if($length < $minLen){
	    print STDERR "[WARNING: Length_less_than_minimum]\t".$name."\t".$length."\n";
	    next;
	}

	unless ($calcGC){
		print "$name\t$length\n" ;
	}
	else{
		my ($g, $c);
		$seq=uc($seq);
	    while ( $seq =~ /G/ig ) { $g++ }
	    while ( $seq =~ /C/ig ) { $c++ }

		my $GC = (($g+$c)/$length)*100;
		my $printGC = sprintf( "%.4f", $GC);
		print "$name\t$printGC\t$length\n";
	}
}
close $CONTIGS;

sub help{
	system('perldoc', $0);
	exit;
}

exit;


