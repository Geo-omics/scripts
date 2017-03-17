#!/user/bin/perl

=head1 Usage

	perl extractSubSeq.pl -f <fasta file> -tsv <custom tab-delimited file>
	OR
	perl extractSubSeq.pl -f <fasta file> -blast <blast file> -query
	OR
	perl extractSubSeq.pl -f <fasta file> -gff <gff3 file>

=head2 Optional

=head3 For Blast Output

	-query	Do you want sub sequences of your queries?
	-subj	Do you want sub sequences of your matches/subjects?
	
=head3 For GFF3

	-map	If you wish to replace the names of the contigs, give it a "map" file; 
		column1=name you wish to change to
		column2=name of the contig in the GFF file

	-header	only use the tag from attributes column as sequence header (example: "locus_tag" or "id"); default: scaffoldName__locus_tag=####__[start-end]
	-feat or feature_type	only extract genes of a certain type (example: "cds" or "exon"); default: all.
	
=head3 General Options
	
	-start <start coordinate coloumn Number>; default: 7 for blast query; 9 for blast subject; 2 for custom tabbed file; 4 for GFF3
	-stop <stop coordinate coloumn Number>; default: 8 for blast query; 10 for blast subject; 3 for custom tabbed file; 5 for GFF3
	-o <output fileName>; default processID.fasta
	-n <Column number>; default: 1
	-l <set minimum alignment length>; default: 0
	-p <set minimum % identity>; default: 0
	-s <set minimum bit score>; default: 0
	-top flag to only get the top hit for each query.
	
=head1 Author

	Sunit Jain, sunitj-at-umich-dot-edu
	April 2010

=cut

use strict;
use Getopt::Long;

my ($fasta, $needsSubj, $needsQuery);
my ($isTabbedFile,$isBlastOut,$isGff, $map);
my ($start, $stop,$topHitOnly,$featType,$headerType);
my $setLen=0;
my $setBS=0;
my $setPid=0;
my $out = $$."_subSeqs.fasta";
my $version="extractSubSeqs.pl\tv0.3.1";

GetOptions(
	'f|fasta:s'=>\$fasta,
	't|tsv|list:s'=>\$isTabbedFile,
	'blast:s'=>\$isBlastOut,
	'gff:s'=>\$isGff,
	'map:s'=>\$map,
	'start:i'=>\$start,
	'stop:i'=>\$stop,
	'query'=>\$needsQuery,
	'subj'=>\$needsSubj,
	'l|len:i'=>\$setLen,
	'p|pid:i'=>\$setPid,
	's|bitscore:i'=>\$setBS,
	'top'=>\$topHitOnly,
	'o|out:s'=>\$out,
	'feat|feature_type:s'=>\$featType,
	'header:s'=>\$headerType,
	'h|help'=> sub{print $version."\n"; system("perldoc $0 \| cat"); exit;},
	'v|version'=> sub{print "# $version\n"; exit;},
);

print "# $version\n";

my $tsv;
if(! $isTabbedFile && ! $isBlastOut && ! $isGff){system('perldoc', $0); exit;}
if($headerType) {$headerType=lc($headerType)}
if($featType) {$featType=lc($featType)}
if ($isTabbedFile){ 
	$tsv = $isTabbedFile;
	if (!$start && !$stop){
		warn "[WARNING] Need start and stop columns\nAssuming Start as col 2 and Stop as col 3:\n";
		$start=2;
		$stop=3;
	}
}
elsif ($isBlastOut){
	$tsv = $isBlastOut;
	if ($needsQuery){
		$start=7;
		$stop=8;
	}
	elsif($needsSubj){
		$start=9;
		$stop=10;
	}
	else{
		die "[ERROR] You need to specify, which column (-query or -subj)I need to look at\nSee '-h' for help on how to use this script\n";
	}
}
elsif ($isGff){
	$tsv = $isGff;
	$start=4;
	$stop=5;
}

$start--;
$stop--;

my %coord;
open( TSV, $tsv) || die "[error] $tsv: $! \n";
while (my $line=<TSV>){
	next if $line=~ /^#/;
	chomp $line;
	$line=~ s/\r//g;
	next unless $line;
	
	my $xy;
	my @cols= split(/\t/, $line);
	if ($isBlastOut){ 
		my ($alnLen, $pid, $bs, $q, $s)=&parseBlastOut($line); 
		next if ($alnLen < $setLen);
		next if ($pid < $setPid);
		next if ($bs < $setBS);

		if ($needsQuery){
			next if ($topHitOnly && $coord{$q});
			$xy=$s."\t".$cols[$start]."\t".$cols[$stop];
			push(@{$coord{$q}},$xy);
		}
		elsif($needsSubj){
			next if ($topHitOnly && $coord{$s});
			$xy=$q."\t".$cols[$start]."\t".$cols[$stop];
			push(@{$coord{$s}},$xy);
		}
	}
	elsif($isGff){
		my $name= &parseGFF3($line);
		next unless $name;
		my $contig=$cols[0];
		$xy=$name."\t".$cols[$start]."\t".$cols[$stop];
		push(@{$coord{$contig}},$xy);
	}
	elsif($isTabbedFile){
		my $name=$cols[0];
		$xy=$cols[$start]."\t".$cols[$stop];
		push(@{$coord{$name}},$xy);
	}
}

my %mapping;
if($map){
	open(MAP, $map)|| die $!;
	while (my $line=<MAP>){
		next if $line=~ /^#/;
		chomp $line;
		$line=~ s/\r//;
		next unless $line;
	
		my($origName, $givenName)=split(/\t/, $line);
		$mapping{$givenName}=$origName
	}
	close MAP;
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
	my $seqLen=length($sequence);
	my @wholeSeq=split(//, $sequence);
	if (($name=~ /\s/g) && ($isBlastOut)){ 
		my $etc;
		($name, $etc)=split(/\s/, $name);
	}
	if ($coord{$name}){
		foreach my $pos(@{$coord{$name}}){
			my (@posStuff)=split(/\t/, $pos);
			my ($match, $posStart, $posStop, $printLine);
			if($isBlastOut){
				$match=$posStuff[0];
				$posStart=$posStuff[1];
				$posStop=$posStuff[2];
				$printLine.=">".$name."_vs_".$match;
			}
			elsif($isGff){
				if ($mapping{$name}){ $match=$mapping{$name}}else{$match=$name};
				my $attribute=$posStuff[0];
				$posStart=$posStuff[1];
				$posStop=$posStuff[2];
				if($headerType){
					$printLine.=">".$attribute
				}
				else{
					$printLine.=">".$match."__locus_tag=".$attribute;
				}
			}
			elsif($isTabbedFile){
				$posStart=$posStuff[0];
				$posStop=$posStuff[1];
				$printLine.=">".$name;
			}
			$printLine.="__[".$posStart."-".$posStop."]" unless ($headerType);
			$printLine.="\n";
			my ($begin, $end)=sort{$a <=> $b}($posStart, $posStop);
			if($end > $seqLen){
				print STDERR "The 'stop coordinate' [".$end."] is out of bounds. Length of the sequence [".$seqLen."]\n";
				die "[FATAL] $printLine\n"
			}
			$begin--;
			$end--;
			print OUT $printLine;
			print OUT @wholeSeq[$begin..$end];
			print OUT "\n";
		}
	}
}
close FASTA;
close OUT;

sub parseGFF3{
#http://gmod.org/wiki/GFF
	my $line=shift;
	my @cols=split(/\t/, $line);

	my(@attributes)=split(/\;/, $cols[-1]);
	
	if ($featType){
		return unless (lc($cols[2]) eq lc($featType));
	}

	my ($locusID,%attribs);
	foreach my $att(@attributes){
		my @a=split(/\=/, $att);
		my $value=join("=",@a[1..$#a]);
		$attribs{lc($a[0])}=$value;
	}

	if(($headerType) && ($attribs{lc($headerType)})){
		$locusID = $attribs{lc($headerType)};
	}
	elsif(! $headerType){
		if( $attribs{"locus_tag"}){
			$locusID=$attribs{"locus_tag"}
		}
		else{
			if ($attribs{"parent"}){
				$locusID=$attribs{"parent"}."_exon"
			}
			elsif($attribs{"repeat"}){
				$locusID=$attribs{"repeat"} if ($attribs{"rpt_type"}); # rpt_type=CRISPR;rpt_unit=13023..13055
				$locusID.="[ ".$attribs{"rpt_unit"}." ]" if ($attribs{"rpt_unit"});
			}
			else{
				$locusID=$attribs{"id"}."_".$cols[2];
			}
		}
	}
	
	return $locusID;
}

sub parseBlastOut{
#	0		1		2	3		4		5		6		7		8		9		10		11
#	query	sbjct	%id	a_len	m_cnt	g_cnt	q_start	q_end	s_start	s_end	e_val	b_score
	my $line=shift;
	my(@cols)=split(/\t/, $line);
	my $aLen=$cols[3];
	my $mis=$cols[4];
	my $gaps=$cols[5];
	my $alnLen= $aLen - ($mis + $gaps);
	my $pid= int($cols[2] + 0.5);
	my $bs=int($cols[11] + 0.5);
	
	return ($alnLen, $pid, $bs, $cols[0], $cols[1]);
}

