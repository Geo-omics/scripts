 #!/user/bin/perl -w

=head1 NAME

	extract blast hits of interest

=head1 Description

	-You give list of Query/Subject names you're interested in and the blast output.
	-Script gives you a sub-set of the blast result that contain the Query or Subjects of Interest.
	-You may adjust the amount of detail you want in you output.
	-Works for tabular Blast output (format 6 and 7 for 2.2.22 and above and for 8 and 9 for the versions 2.2.20 and below)

=head1 Usage

	perl extract_Blast_Hits_Of_Interest.pl -b <blast output file> -l <list of hits of interest>

=head1 Options

	-o:		output file; 
	-c:		where to look, for the names in the list:
			'q': look in Query Column;
			's': look in Subj Column ;
			'b': look in Both Columns; 
	-brief:		level of detail you want in the output file. This is a boolean flag, its presence in the command line would mean you just want the corresponding subj/query as your output
			The default behavior (i.e -brief is not present) of the script prints the output in the blast tabular format.
	-mp:		minimum % id
	-mal:		minimum alignment length

=head3 Defaults

	-o:		<pid>_BlastHitsOfInterest.out
	-c:		q
	-mp:		0
	-mal:		0

=head1 Author

	Sunit Jain, May 2011
	sunitj-at-umich-dot-edu

=cut

use strict;
use Getopt::Long;

## Set Options ##

my($bOut,$listFile, $cNum, $noDetail);
my $out= $$."_BlastHitsOfInterest.out";
my $col= "q";
my $setPer= 0;
my $setLen= 0;
my $setBS= 0;
GetOptions(
	'b:s'	=>	\$bOut,
	'l:s'	=>	\$listFile,
	'o:s'	=> 	\$out,
	'c:s'	=>	\$col,
	'brief'		=>	\$noDetail,
	'mp:f'	=>	\$setPer,
	'mal:i'	=>	\$setLen,
	'mbs:f'	=>	\$setBS,
	'h|help'=>	sub{system('perldoc', $0); exit;},
);

## Checks ##

if (! $bOut || ! $listFile){system('perldoc', $0); exit;}

unless (lc($col) eq "q" || lc($col) eq "s" || lc($col) eq "b"){
	print "The script does not recognize \'$col\' as an input option for \'-c\' flag\n";
	exit;
}

## Main ##

open(LON, $listFile) || die "[err] $listFile: $! \n";
my %index;
while (my $l=<LON>){
	chomp $l;
    next unless $l;
	next if $l=~ m/^#/;
	$l=~ s/^>// if ($l=~ m/^>/);
	$l=uc($l);
	if ($l=~ m/GD[5,6]pt[1,2]/i){$l=~ s/pt\d//i;}
	$index{$l}++;
}
close LON;

print "#Searching for ".keys(%index)." Sequence Names.\n";

my %found;
&extract;
my $notFound=keys(%index) - keys(%found);
if ($notFound >0){
	print "# $notFound items not found.\n";
	foreach my $i(keys(%index)){
		print $i."\n" unless ($found{$i});
	}
}
print "#".keys(%found)." Sequence Names found\n";
## Sub-Routines ##

sub extract {
## Input i.e. BLASTALL Output Format: m8/m9l, & 6/7(blastn/p/x) ##
#	0		1		2	3		4		5		6		7		8		9		10		11
#	query	sbjct	%id	a_len	m_cnt	g_cnt	q_start	q_end	s_start	s_end	e_val	b_score
	
	my $count=0;
	open (OUT, ">".$out);
	print OUT "#min\%ID=$setPer\tminAL=$setLen\tminBitScore=$setBS\n";
	if ($noDetail){
		print OUT "\#Query\n" if (lc($col) eq "s");
		print OUT "\#Subj\n" if (lc($col) eq "q");
		print OUT "\#Query\tSubj\n" if (lc($col) eq "b");
	}
	open (BOUT, "$bOut") || die "[err]: $bOut: $! \n";

	while(my $line=<BOUT>){
		next if ($line=~ m/^\#/);
		chomp($line);
		next unless $line;
		my @blastOut=split(/\t/, $line);
		chomp(@blastOut);
		my $query=uc($blastOut[0]);
		my $subj=uc($blastOut[1]);
		my $per=$blastOut[2];
		my $score=$blastOut[-1];
		chomp($query,$subj,$per,$score);
# Get all hits that clear the thresholds
		my ($q,$s);
		my $flag=0;
		if (lc($col) eq "q"){
			$q= parseName($query);
			$flag = 1;
			next if (! $q);
		}
		elsif(lc($col) eq "s"){	
			$s= parseName($subj);
			$flag = 1;
			next if (! $s);
		}
		elsif(lc($col) eq "b"){
			$q= parseName($query);
			$s= parseName($subj);
			$flag = 1;
			if (! $q && ! $s){ next; }
		}
		else{
			$flag = 0;
			next;
		}

		if (($score >= $setBS) && ($flag == 1)){	
			&printStuff($line);
		}		

	}
	close BOUT;
	close OUT;
#	print "$count results found!\n";
}

sub parseName{
	my $name=shift;
	chomp $name;
	if ($index{$name}){
		$found{$name}++;	
		return $name;
	}
	$name=~ s/[^A-Z0-9_]/ /g;
	my @nameParts=split(/ /, $name);
	foreach my $n(@nameParts){
		if ($index{$n}){ $found{$n}++; return $n; }
		else{next;}
	}
}

sub printStuff{
	my $stuff=shift;
	if ($noDetail){
		my @lineParts= split(/\t/, $stuff);
		print OUT "$lineParts[0]\n" if (lc($col) eq "s");
		print OUT "$lineParts[1]\n" if (lc($col) eq "q");
		print OUT "$lineParts[0]\t$lineParts[1]\n" if (lc($col) eq "b");
	}
	else{
		print OUT $stuff."\n";
	}	
}
