#!/user/bin/perl

=head1 Description
	
	The Script is a modification of Ryan's 'highestHitsPerQuery.pl' and 'coolScript.pl'.
	It scores subjects based on it's popularity with the queries in the blast output.

=head2 Usage

	perl mapper_v052.pl [-b: blast output (m8/m9 for blastall & outfmt 6/7 for blast 2.2.24+)] [-q: query file for blast]
	[-o: output filename] [-p: min % id] [-l: min alignment length] [-d: % deviation allowed from top bit score]
	[-s: min bit score] [-v: version]

=head2 Defaults
	
	-o	'out.tsv'
	-p	00%
	-l	00 (bases)
	-d	00%
	-s	0 (bits)

=head1 Author

	Sunit Jain (June, 2011)
	sunitj-at-umich-dot-edu

=cut


use strict;
use Getopt::Long;

## Set Defaults ##

my $bOut;
my $qFasta;
my $setPer=0;
my $setLen=0;
my $setDev=0;
my $setBscore=0;
my $out;
my $version= "0.5.2";

GetOptions(
	'b|blastout:s'=>\$bOut,
	'q|query:s'=>\$qFasta,
	'o|output:s'=>\$out,
	'p|per:f'=>\$setPer,
	'l|alen:i'=>\$setLen,
	'd|deviation:f'=>\$setDev,
	's|bit_score:f'=>\$setBscore,
	'h|help'=> sub{system('perldoc', $0); exit;},
	'v|version'=>sub{print "Version: mapper_".$version."\n"; exit;}
);

## Checks ##

if (! $bOut){ print "Blast output file not found. See help (-h) for details.\n"; exit;}
if (! $qFasta){print "Query file not found. See help (-h) for details.\n"; exit;}
if (! $out){ $out="out_pid".$$."_P".$setPer.".mapped";}

## MAIN ##

my $sd=$setDev/100;
my (%hs, %queryCount,%bestQueryScore, %tally, %vote);
my ($gd5Count, $gd6Count, $gd7Count, $gd8Count, $totalSeqs);
print "Applying Thresholds...\n";
&setThresholds;

print "Scoring...\n";
&getCounts;
&score;

## Sub-Routines ##

sub setThresholds {
## Sub-Routine Input i.e. BLASTALL Output Format: m8, & m6(blastn/p/x) ##
#	0		1		2	3		4		5		6		7		8		9		10		11
#	query	sbjct	%id	a_len	m_cnt	g_cnt	q_start	q_end	s_start	s_end	e_val	b_score
	my %scoreTrack=();
	my %queryScores=();

	open (BOUT, "$bOut") || die "ERROR: $bOut \n:".$!;

	while(my $line=<BOUT>){
		next if ($line=~ m/^\#/);
		chomp($line);
		my @blastOut=split(/\t/, $line);
		chomp(@blastOut);
		my $query=uc($blastOut[0]);
		my $subj=uc($blastOut[1]);
		my $per=$blastOut[2];
		my $aLen=$blastOut[3];
		my $score=$blastOut[-1];
		chomp($query,$subj,$per,$score, $aLen);
#		Get all hits that clear the thresholds
		my $score= int($score + 0.5);
		if(($aLen >= $setLen) && ($per >= $setPer) && ($score >= $setBscore)){
			push(@{$scoreTrack{$query}{$score}},$subj);
			push (@{$queryScores{$query}}, $score);
		}
		else{
			next;
		}
	}
print "\tBlast output parsed\n";
print "\tProcessing...\n";
# Get more useable data structures from the ones created above.
	foreach my $q(keys %queryScores){
#		get the highest score for each query; descending order
		my @sortedScores=sort{$b <=> $a } @{$queryScores{$q}};
		my $highestScore=$sortedScores[0];
		$bestQueryScore{$q}=$highestScore;
#		Count the number of times the query appears with a bit-score above the set deviation.		
		my $floor=$highestScore - ($highestScore * $sd);
		foreach my $bs(@{$queryScores{$q}}){
			if ($bs >= $floor){
				foreach my $s(@{$scoreTrack{$q}{$bs}}){
#					which queries contributed to the score of the subject; the value of the hash is just used to keep track of the highest score for the subj-query pair.
					$hs{$q}{$s}=$bs;
#					freq of subj-query pair.
					$tally{$s}++;
#					this query is present with a bit score >= the floor these many times.
					$queryCount{$q}++;
				}
			}			
			else{
				next;
			}
		}
	}
	undef %scoreTrack;
	undef %queryScores;
	close BOUT;

}

sub getCounts{
	open (QUERY, $qFasta) || die "[ERROR] $qFasta $! \n";
	$/= ">";

	while (my $b = <QUERY>) {
		chomp $b;
		next unless $b;
		my ($name, @sequence) = split (/\n/, $b);
		my ($gd, @descriptor)= split(/\_/, $name);
		$gd5Count++ if (uc($gd) eq "GD1");
		$gd6Count++ if (uc($gd) eq "GD2");
		$gd7Count++ if (uc($gd) eq "GD9");
		$gd8Count++ if (uc($gd) eq "GD10");
		$totalSeqs++;
#		$length = length($seq);
	}
	close QUERY;
	$/="\n";
}

sub score{
# coolScript.pl by Ryan

	my($qList, @etc)=split(/\./, $out);
print "\tWriting \'List\' file to:\t$qList.list\n";	

	open (LIST, ">".$qList.".list");
	my $qCount=0;
	foreach my $q(keys(%queryCount)){
		print LIST $q."\n";
		$qCount++;
	}
	close LIST;
	undef %queryCount;

	my %bestBitScores=();
	my %subjBitScores=();
	my %out=();
	my %log=();
	foreach my $query(keys %hs){
		my $num= keys %{$hs{$query}};
		my $bitScore= $bestQueryScore{$query};

		my @subj2ThisQuery=keys %{$hs{$query}};
		my %subjTally;
		foreach my $s2q(@subj2ThisQuery){
			$subjTally{$s2q}=$tally{$s2q};
		}
# sort the #tally to get the subj with the highest freq of best-hits in the dataset;
		my @sortedSubj;
		if (@subj2ThisQuery >1){
			@sortedSubj = sort{$subjTally{$b} <=> $subjTally{$a}} keys %subjTally;
		}
		else{
			$sortedSubj[0]=$subj2ThisQuery[0];
		}
		my($gd, @number);
		if ($query =~ m/\_/g){
			($gd, @number)=split(/_/, $query);
			$gd=uc($gd);
			my $highestFreq=$tally{$sortedSubj[0]};
			my $s=$sortedSubj[0];
#			foreach my $s(@sortedSubj){
#				if ($subjTally{$s} == $highestFreq){
					if ($gd=~ m/GD[1,2,9,10]/){
						if ($gd eq "GD1"){ $vote{$s}{'GD1'}++;}
						elsif ($gd eq "GD2"){ $vote{$s}{'GD2'}++;}
						elsif ($gd eq "GD9"){ $vote{$s}{'GD9'}++;}
						elsif ($gd eq "GD10"){ $vote{$s}{'GD10'}++;}
					}
					else{
						print "1a\t".$query."\n";
					}
#				}
#			}
		}
		else{ print "1b\t".$query."\n";}

		foreach my $subj(keys %{$hs{$query}}){
			my $score=1/$num;
			$out{$subj}{'Score'}+=$score;
			if ($gd eq "GD1"){ $out{$subj}{'ScoreGD1'}+=1/$num; $out{$subj}{'GD1'}++; }
			elsif ($gd eq "GD2"){ $out{$subj}{'ScoreGD2'}+=1/$num; $out{$subj}{'GD2'}++; }
			elsif ($gd eq "GD9"){ $out{$subj}{'ScoreGD9'}+=1/$num; $out{$subj}{'GD9'}++;}
			elsif ($gd eq "GD10"){ $out{$subj}{'ScoreGD10'}+=1/$num; $out{$subj}{'GD10'}++;}
			else{ print $query."\t";next;}
			push(@{$log{$subj}}, $query);
			$bestBitScores{$subj}{$bitScore}=$query;
			push (@{$subjBitScores{$subj}},$bitScore);
		}
	}
	close LOG;
	undef %hs;
	undef %bestQueryScore;


	undef %tally;

print "\tWriting \'Mapped\' file to:\t$out\n";
#print keys(%vote)."\t";
#print keys(%out)."\n";
	open (OUT, ">$out");
	print OUT "\# Blast File:\t$bOut\t\tTotalSeq:\t$totalSeqs\n";
	print OUT "\# \% ID:\t$setPer\t\tTotalGD1\t$gd5Count\n";
	print OUT "\# min_length(bases):\t$setLen\t\tTotalGD2\t$gd6Count\n";
	print OUT "\# Bit Score min:\t$setBscore\t\tTotalGD9\t$gd7Count\n";
	print OUT "\# \% BitScore deviation:\t$setDev\t\tTotalGD10\t$gd8Count\n";
	print OUT "\# Reference Region\tTally\tGD1Tally\tGD2Tally\tGD9Tally\tGD10Tally\tScoring\tTop Transcript\tTop BitScore\tTotal\tScore\t\%Total\tGD1\tScoreGD1\t\%GD1\tGD2\tScoreGD2\t\%GD2\tGD9\tScoreGD9\t\%GD9\tGD10\tScoreGD10\t\%GD10\n";
	
	foreach my $s(keys %out){
		my $bestQuery;
		my $bestBitScore;
		my @sortedBitScores= sort{$b <=> $a} @{$subjBitScores{$s}};
		$bestBitScore=$sortedBitScores[0];
		$bestQuery=$bestBitScores{$s}{$bestBitScore};
		
		my $tallyTotal=$vote{$s}{'GD1'}+$vote{$s}{'GD2'}+$vote{$s}{'GD9'}+$vote{$s}{'GD10'};
		my $gd5Tally=$vote{$s}{'GD1'};
		my $gd6Tally=$vote{$s}{'GD2'};
		my $gd7Tally=$vote{$s}{'GD9'};
		my $gd8Tally=$vote{$s}{'GD10'};

		my $total=$out{$s}{'GD1'} + $out{$s}{'GD2'} + $out{$s}{'GD9'} + $out{$s}{'GD10'};
		my $perTotal= ($total / $totalSeqs) * 100;
		my $gd5NScore=($out{$s}{'GD1'} / $gd5Count) * 100;
		my $gd6NScore=($out{$s}{'GD2'} / $gd6Count) * 100;
		my $gd7NScore=($out{$s}{'GD9'} / $gd7Count) * 100;
		my $gd8NScore=($out{$s}{'GD10'} / $gd8Count) * 100;

		print OUT $s."\t";
		print OUT $tallyTotal."\t";
		print OUT $gd5Tally."\t".$gd6Tally."\t".$gd7Tally."\t".$gd8Tally."\t\*\t";
		print OUT $bestQuery."\t".$bestBitScore."\t";
		print OUT $total."\t".$out{$s}{'Score'}."\t".$perTotal."\t";
		print OUT $out{$s}{'GD1'}."\t".$out{$s}{'ScoreGD1'}."\t".$gd5NScore."\t";
		print OUT $out{$s}{'GD2'}."\t".$out{$s}{'ScoreGD2'}."\t".$gd6NScore."\t";
		print OUT $out{$s}{'GD9'}."\t".$out{$s}{'ScoreGD9'}."\t".$gd7NScore."\t";
		print OUT $out{$s}{'GD10'}."\t".$out{$s}{'ScoreGD10'}."\t".$gd8NScore."\n";

	}
	close OUT;
	undef %out;
	undef %bestBitScores;
	undef %subjBitScores;

print "\tWriting \'Log\' file to:\t$qList.log\n";

	open (LOG, ">".$qList.".log");
	print LOG "\#Subj\tQuery1\tQuery2...etc.\n";
	foreach my $s(keys %log){
		print LOG $s."\t";
		foreach my $q(@{$log{$s}}){
			print LOG $q."\t";
		}
		print LOG "\n";
	}
	print "Total Unique Queries considered:\t".$qCount."\n";
}
