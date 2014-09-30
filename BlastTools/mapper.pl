#!/user/bin/perl

=head1 Description
	
	The Script is a modification of Ryan's 'highestHitsPerQuery.pl' and 'coolScript.pl'.
	It scores subjects based on it's popularity with the queries in the blast output.

=head2 Usage

	perl mapper.pl [-b: blast output (m8/m9 for blastall & outfmt 6/7 for blast 2.2.24+)] [-q: query file for blast]

=head2 Options
	
	-o	=	Output Prefix; default='out'
	-p	=	minimum percent ID; default=95%
	-l	=	minimum alignment length= alignment Length - (mismatches + gaps)	default=40 (bases)
	-d	=	% Deviation from the top bit score allowed, choose a value between 0-100; default=0
	-s	=	minimum bit score; default=0 (bits)

=head2 Output files

	.log file: contains the 'clusters'.
	.list file: contains a list all unique subjects that have met the criteria
	.mapped file OR output file mentioned in command line: mapped analysis of your blast output.

=head1 Author

	Sunit Jain (July, 2011)
	sunitj-at-umich-dot-edu

=cut


use strict;
use Getopt::Long;
use File::Basename;

## Set Defaults ##

my $bOut;
my $qFasta;
my $setPer=95;
my $setLen=40;
my $setDev=0;
my $setBscore=0;
my $out;
my $version= "0.4.6";

GetOptions(
	'b|blastout:s'=>\$bOut,
	'q|query:s'=>\$qFasta,
	'o|prefix:s'=>\$out,
	'p|per:f'=>\$setPer,
	'l|alen:i'=>\$setLen,
	'd|deviation:f'=>\$setDev,
	's|bit_score:f'=>\$setBscore,
	'h|help'=> sub{system("perldoc $0 \| cat"); exit;},
	'v|version'=>sub{print "Version: $0\tv$version\n"; exit;}
);

## Checks ##

if (! $bOut){ print "Blast output file not found. See help (-h) for details.\n"; exit;}
if (! $qFasta){print "Query file not found. See help (-h) for details.\n"; exit;}
if (! $out){ $out="out_pid".$$."_P".$setPer.".mapped";}

## MAIN ##

# my($qList, @etc)=fileparse($out, ".mapped");
my $qList=$out;
$out=$out.".mapped";

my $sd=$setDev/100;
my (%hs, %queryCount,%bestQueryScore);
my ($totalSeqs);
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
		my $query=$blastOut[0];
		my $subj=$blastOut[1];
		my $per=$blastOut[2];
		my $aLen=$blastOut[3];
		my $mm=$blastOut[4];
		my $gaps=$blastOut[5];
		my $score=$blastOut[11];
		chomp($query,$subj,$per,$aLen,$gaps,$mm,$score);
#		Get all hits that clear the thresholds
		my $tempLen=$aLen-$mm-$gaps;
		my $score= int($score + 0.5);
		if(($tempLen >= $setLen) && ($per >= $setPer) && ($score >= $setBscore)){
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
#					tally of which queries contributed to the score of the subject; the value of the hash is just used to keep track of the highest score for the subj-query pair.
					$hs{$q}{$s}=$bs;
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
		$totalSeqs++;
	}
	close QUERY;
	$/="\n";
}
sub score{
# coolScript.pl by Ryan
	print "\tWriting \'List\' file to:\t$qList.list\n";	# unique hits

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
		my $highestScore=-1;
		my $num= keys %{$hs{$query}};
		my $bitScore= $bestQueryScore{$query};

		foreach my $subj(keys %{$hs{$query}}){
			my $score=1/$num;
			$out{$subj}{'Score'}+=$score;
			$out{$subj}{'Total'}++;
			push(@{$log{$subj}}, $query);
			$bestBitScores{$subj}{$bitScore}=$query;
			push (@{$subjBitScores{$subj}},$bitScore);
		}
	}
	undef %hs;
	undef %bestQueryScore;

	print "\tWriting \'Mapped\' file to:\t$out\n";
	open (OUT, ">$out");
	print OUT "\# Blast File:\t$bOut\n";
	print OUT "\# \% ID:\t$setPer\n";
	print OUT "\# min_length(bases):\t$setLen\n";
	print OUT "\# Bit Score min:\t$setBscore\n";
	print OUT "\# \% BitScore deviation:\t$setDev\t\t\tTotalSeq:\t$totalSeqs\n";
	print OUT "\# Reference Region\tTop Query\tTop BitScore\tTotal\tScore\tTotal(in \%)\n";
	
	foreach my $s(keys %out){
		my @sortedBitScores= sort{$b <=> $a} @{$subjBitScores{$s}};
		my $bestBitScore=$sortedBitScores[0];
		my $bestQuery=$bestBitScores{$s}{$bestBitScore};
		
		my $total=$out{$s}{'Total'};
		my $perTotal= ($total / $totalSeqs) * 100;

		print OUT $s."\t";
		print OUT $bestQuery."\t".$bestBitScore."\t";
		print OUT $total."\t".$out{$s}{'Score'}."\t".$perTotal."\n";
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
	my $perc_queries_considered=sprintf( "%.2f",(($qCount/$totalSeqs)*100));
	print "Unique Queries considered:\t$perc_queries_considered\% [ $qCount out of $totalSeqs queries ]\n";
}
