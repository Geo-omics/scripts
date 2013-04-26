#!/user/bin/perl

=head1 Usage

	perl extractSubSeq.pl -f <fasta file> -t <custom tab-delimited file>
	OR
	perl extractSubSeq.pl -f <fasta file> -blast <blast file>

=head2 Optional

	-start <start coordinate coloumn Number>; default: 7 for blast; and 2 for custom tabbed file
	-stop <stop coordinate coloumn Number>; default 8 for blast; and 3 for custom tabbed file
	-o <output fileName>; default processID.fasta
	-n <Column number>; default: 1
	-l <set minimum alignment length>; default: 0
	-p <set minimum % identity>; default: 0
	-s <set minimum bit score>; default: 0
	-top flag to only get the top hit for each query.

=cut

use strict;
use Getopt::Long;

my $fasta;
my $tabbedFile;
my $seqNameCol=1;
my $start;
my $stop;
my $setLen=0;
my $setBS=0;
my $isBlastOut;
my $setPid=0;
my $topHitOnly;
my $out = $$."_subSeqs.fasta";

GetOptions(
	'f|fasta:s'=>\$fasta,
	't|tsv:s'=>\$tabbedFile,
	'n|name:i'=>\$seqNameCol,
	'start:i'=>\$start,
	'stop:i'=>\$stop,
	'l|len:i'=>\$setLen,
	'p|pid:i'=>\$setPid,
	's|bitscore:i'=>\$setBS,
	'blast:s'=>\$isBlastOut,
	'top'=>\$topHitOnly,
	'o|out:s'=>\$out,
	'h|help'=> sub{system('perldoc', $0); exit;},
);

my $tsv;
if(! $tabbedFile && ! $isBlastOut){system('perldoc', $0); exit;}
elsif( ! $tabbedFile && $isBlastOut){ $tsv = $isBlastOut; }
elsif( $tabbedFile && ! $isBlastOut){ $tsv = $tabbedFile; }
elsif($tabbedFile && $isBlastOut){print "Decide!! Is it a blast output OR a custom tab-delimited file?!\nAssuming it's a BLAST output!\n";}

if ($tabbedFile && !$start && !$stop){print "[warn] Need start and stop columns\nAssuming Start as col 2 and Stop as col 3:\n"; $start=2; $stop=3;}
if ($isBlastOut && !$start && !$stop){print "[warn] Need start and stop columns\nAssuming Start as col 7 and Stop as col 8:\n"; $start=7; $stop=8;}

$seqNameCol--;
$start--;
$stop--;

my %coord;
open( TSV, $tsv) || die "[error] $tsv: $! \n";
while (my $line=<TSV>){
	next if $line=~ m/^#/;
	chomp $line;
	$line=~ s/\r//g;
	next unless $line;
	
	my @cols= split(/\t/, $line);
	my ($alnLen, $pid, $bs);
	if ($isBlastOut){ 
		($alnLen, $pid, $bs)=&handleBlastOut(\@cols); 
		next if ($alnLen < $setLen);
		next if ($pid < $setPid);
		next if ($bs < $setBS);
	}
	my $xy= $cols[$start]."\t".$cols[$stop];
	if($topHitOnly){
		next if ($coord{$cols[$seqNameCol]}{$cols[0]});
	}
	push(@{$coord{$cols[$seqNameCol]}{$cols[0]}},$xy);
}

open (FASTA, $fasta)|| die "[error] $fasta: $!\n";
open (OUT, ">".$out);
$/=">";
while (my $line=<FASTA>){
	next if $line=~ m/^#/;
	chomp $line;
	$line=~ s/\r//g;
	next unless $line;

	my ($name, @seq)=split(/\n/, $line);
	my $sequence=join("", @seq);
	my @wholeSeq=split(//, $sequence);
	if (($name=~ m/\s/g)&& ($isBlastOut)){ 
		my $etc;
		($name, $etc)=split(/\s/, $name);
	}
	if ($coord{$name}){
		foreach my $hh(keys %{$coord{$name}}){
			foreach my $hha(@{$coord{$name}{$hh}}){
				my ($st, $sp)=split(/\t/, $hha);
				my ($start, $stop, $strand);
				if ($st > $sp){
					$start=$sp;
					$stop=$st;
					$strand="(-)";
				}
				else{
					$start=$st;
					$stop=$sp;
					$strand="(+)";
				}
				print OUT ">".$hh."_vs_".$name."__[".$start."-".$stop."]_".$strand."\n";
				print OUT @wholeSeq[$start..$stop];
				print OUT "\n";
			}
		}
	}
}


sub handleBlastOut{
#	0		1		2	3		4		5		6		7		8		9		10		11
#	query	sbjct	%id	a_len	m_cnt	g_cnt	q_start	q_end	s_start	s_end	e_val	b_score
	my $c=shift;
	my @cols=@{$c};
	my $aLen=$cols[3];
	my $mis=$cols[4];
	my $gaps=$cols[5];
	my $alnLen= $aLen - ($mis + $gaps);
	my $pid= int($cols[2] + 0.5);
	my $bs=int($cols[-1] + 0.5);
	
	return ($alnLen, $pid, $bs);
}
