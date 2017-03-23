#!/usr/bin/perl
=head1 Description

	Given a fasta file with repeat sequences and a contig fasta file, get the positions of these repeats in the contigs and find the coordinates of spacers, defined by default as sequences between two repeats which are less than 200 bases apart. This definition can be changed using a '-len' parameter.
	NOTE: This script only searches for exact matches (100% over 100% length) and not "similar" sequences (not even 99.99%)

=head1 Usage

	perl crispr_spacer_extractor.pl -c contigs.fasta -r repeats.fasta

=head2 Options

	-r or -repeats	(required)	repeats fasta file
	-c or -contigs	(required)	contigs fasta file
	-p or -prefix	output file prefix	(optional; default= process_id)
	-l or -len	maximum length of a spacer	(optional; default=200)
	-rc	or -revcomp	look for reverse complements of the repeat sequences as well
	
=head1 Author

	Sunit Jain, April 2013
	sunitj@umich.edu
	
=cut

use strict;
use Getopt::Long;

my ($contigs, $repeats, $RC);
my $prefix=$$;
my $maxSpacerLen=200;

GetOptions(
	'r|repeats=s'=>\$repeats,
	'c|contigs=s'=>\$contigs,
	's|l|len:i'=>\$maxSpacerLen,
	'p|prefix:s'=>\$prefix,
	'rc|revcomp'=>\$RC,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

my ($spacerFile,$repOut);
if ($prefix){
	$spacerFile= $prefix.".spacer.coords";
	$repOut= $prefix.".repeat.coords";
}
else{
	$spacerFile= $$.".spacer.coords";
	$repOut=$$.".repeat.coords";
}

$/=">";
my %REPEATS;
open(REP, $repeats) || die "$!";
while(my $line=<REP>){
	chomp $line;
	$line=~ s/\r//g;
	next unless $line;
	my ($repeat, $seq)=split(/\n/, $line);
	next if (! $repeat);
	$REPEATS{$seq}=$repeat;
	if ($RC){
		my $rcSeq=reverseComplement($seq);
		$REPEATS{$rcSeq}=$repeat."_(*)";
	}
}
$/="\n";
close REP;

$/=">";
my ($numRepeats, $numSpacers);
open(CONTIGS, $contigs)|| die "$!";
open(REPOUT, ">".$repOut)|| die $!;
open(SPACER, ">".$spacerFile)|| die $!;
while(my $line=<CONTIGS>){
	chomp $line;
	$line=~ s/\r//g;
	next unless $line;
	my ($header, @sequence)=split(/\n/, $line);
	my $seq=join("", @sequence);
	my %SPACERS;
	foreach my $r(keys(%REPEATS)){
		my @matches;
		push (@matches, match_all_positions($r,$seq));
		my $matchNum=scalar(@matches);
		if ($matchNum > 0){
			$numRepeats++;
			print REPOUT ">".$header."\t".$REPEATS{$r}."\t".$matchNum."\n";
			foreach my $m(@matches){
				my $pos=join("\t",@$m);
				$SPACERS{$pos}=$$m[0];
				print REPOUT $pos."\n";
			}
		}
	}
	my $prevStop=0;
	foreach my $pos(sort{$SPACERS{$a} <=> $SPACERS{$b}} keys %SPACERS ){
		my($start, $stop)=split(/\t/, $pos);
		if(($start-$prevStop) <= $maxSpacerLen){
			$numSpacers++;
			print SPACER $header."\t".$prevStop."\t".$start."\n"
		}
		$prevStop=$stop;
	}	
}
$/="\n";
close CONTIGS;
close REPOUT;
close SPACER;
print "Repeats Found:\t$numRepeats\n";
print "Spacers Found:\t$numSpacers\n";

if (-e "extractSubSeq.pl"){
	system("perl extractSubSeq.pl -f $contigs -t $spacerFile -o $prefix.spacers.fasta -start 2 -stop 3");
}
elsif(-e "/geomicro/data1/COMMON/scripts/extractSubSeq.pl"){
	system("perl /geomicro/data1/COMMON/scripts/extractSubSeq.pl -f $contigs -t $spacerFile -o $prefix.spacers.fasta -start 2 -stop 3");
}
else{
	warn "Conldn't find the extractSubSeqs.pl script.\nIf you need the fasta file for your spacer sequences, please run the script on your contigs and spacer.coords file\n";
}

exit 0;
###########################################################################################

sub match_positions {
    my ($regex, $string) = @_;
    return if not $string =~ /$regex/;
    return ($-[0], $+[0]);
}
sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, [ $-[0], $+[0] ];
    }
    return @ret
}
sub reverseComplement{
	my $seq=scalar(@_) ? shift : <IN>;
	chomp $seq;
	my $rSeq=uc(reverse($seq));
	$rSeq=~ tr/GTCA/CAGT/;
	return $rSeq;
}
