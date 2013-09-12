#!/usr/bin/perl

=head2 Description

	Given the names of sequences. extract from a fasta file. Make sure each name is in a seperate line.
	By default, the script prints the sequences with the names in the list.
	This behavior can be changed by adding the '-e' flag which will ask the script to exclude the sequences present in the list

=head2 Examples
	
	extractSeqs.pl -l <list of sequences names> -f <fasta file> -o output.fasta
	OR
	extractSeqs.pl -l <list of sequences names> -fq <fastq file> -o output.fastq
	OR
	extractSeqs.pl -l <list of sequences names> -fq <fastq file> -o output.fasta -outfmt fasta
	
=head2 Optional

	-fq	fastq file
	-fuzzy If you list file contains a part of the name in the sequence file, use this.
		Default: exact match.
	-o	name of the output file; default: processID.fasta.
	-outfmt	fasta; if the output file required is a fasta file when a fastq file is the input.
	-e	exclude list.
	-h	help; this text.

=head2 Author

	Sunit Jain, March 2011
	sunitj [AT] umich [DOT] edu

=head2 Contributors

	Chris Taylor, 2011

=cut

#######################
## MODULES
#######################
use strict;
use Getopt::Long;


#######################
## PARAMETERS
#######################
my $listOfNames;
my $fasta;
my $fastq;
my $excl;
my $fuzzy;
my $out;
my $outfmt;
my $version="extractSeqs.pl v0.4.5";

GetOptions(
	"l|list:s"	=>	\$listOfNames,
	"f|fasta:s"	=>	\$fasta,
	"fq|fastq:s"=>	\$fastq,
	"o|out:s"	=>	\$out,
	"outfmt:s"	=>	\$outfmt,
	"e|exclude"	=>	\$excl,
	"fuzzy"	=> \$fuzzy,
	"h|help"	=>	sub {system('perldoc', $0); exit;},
	"v|version"	=>	sub{print $version."\n"; exit;},
);

#######################
## CHECKS
#######################
if ( ! $listOfNames){
	system('perldoc', $0);
	exit;
}

my $seqFile;
if ( $fasta && ! $fastq){ $seqFile = $fasta }
elsif ( $fastq && ! $fasta){ $seqFile = $fastq }
elsif ( $fastq && $fasta ){ print STDERR "[Warning]: I can only handle one file type at a time; ignoring fastq file..."; $seqFile = $fasta; }
elsif (! $fastq && ! $fasta){ system('perldoc', $0); exit; }

if (! $outfmt){$outfmt=$fastq ? "fastq" : "fasta";}
if (! $out){	$out=$$.".".$outfmt	};

my $printFastaOut=$outfmt eq "fasta" ? "1" : "";

#######################
## MAIN
#######################
open (LON, $listOfNames)|| die "Couldn't open $listOfNames\n";
my %index;
while (my $l=<LON>){
	chomp $l;
	$l=~ s/\r//;
    next unless $l;
	next if $l=~ m/^#/;
	my @etc;	
	($l, @etc)=split(/\s+/, $l); 
	$l=~ s/^>// if ($l=~ m/^>/);
	$l=~ s/^@// if ($l=~ m/^@/);
	$l=uc($l);
	$index{$l}++;
}
close LON;
#######################

print keys(%index)." Sequence Names found.\n";

$/=  $fastq ? "\n" : ">";
my $count=0;
my $total=0;
my $stupidNames=0;
my $prevName="";
my $prevSeq="";
my %seen=();

#######################
open (SEQ, $seqFile) || die "Couldn't open $seqFile\n";
open (OUT, ">".$out);

while (my $line = <SEQ>) {
	next if ($line=~ m/^#/);
    chomp $line;
	$line=~ s/\r//;
    next unless $line;

	my ($name, $seq, $qName, $qual);
	if ($fasta){
		($name, $seq)= &parseFasta($line);
		### To handle rogue '>' in seq headers	
		if (length($seq)==0){
			$stupidNames++;
			$prevSeq=$seq;
			$prevName=$name;
			next;
		}
		if (length($prevSeq) == 0){
			$name = $prevName.$name;
			$prevSeq=$seq;
			$prevName=$name;
		}
		###
	}
	else{
		($name, $seq, $qName, $qual)= &parseFastq($line);	
	}
	
	my $printName=$name;
	$name=uc($name);
	
	my (@nameParts);
	if ($fuzzy){
		@nameParts=&fuzzyMatch($name)
	}
	else{
		@nameParts=split(/ /, $name);
	}
	my $c=0;
	$total++;

	### This is where the magic happens
	if ($excl){
		$count += exclude($printName, \@nameParts, $seq, $qName, $qual);
	}
	else{
		$count += include($printName, \@nameParts, $seq, $qName, $qual);
	}
	###
}
close SEQ;
close OUT;
#######################

### Rant about the rogue '>' in sequence headers
print $stupidNames." Name(s) with \> sign in the header found: \n" if ($stupidNames > 0);
if (scalar(keys %index)-scalar(keys %seen) > 0){
	my $s= scalar(keys %index)-scalar(keys %seen);
	print "$s names caused an error and were not processed further.Here they are:\n";
	foreach my $key(keys %index){
		print $key."\n" unless $seen{$key};
	}
	print "Try using the '-fuzzy' flag.\nWARNING: Make sure you check your output after using this flag. This is still an experimental feature.\n"
}
else{
	print "Got 'em all!!\n";
	if ($excl){
		print $count." Sequence(s) were printed in $out\n";
		print "Exclusion was on.\n";
	}
	else{
		print $count." off ". keys(%index)." Sequence(s) were printed in $out\n";
	}
}
###

#######################
## SUB-ROUTINES
#######################
sub fuzzyMatch{
	my $name=shift;
	## If you're having any issues with the file, it's most likely due to the next line. IF YOU KNOW WHAT YOU'RE DOING try modifying the reg-ex.
	$name=~ s/[^A-Z0-9=_#\/\\\.]/ /g;
	my @nameParts=split(/ /, $name);
	return(@nameParts);
}

sub parseFastq{
	my $line= shift;
	if ($line=~ /^@\w+/){
		$line=~ s/^@//;
		$line=~ s/\r//g;
		my $seqDesc=$line;

		$line=<SEQ>;
		$line=~ s/\r//g;
		chomp $line;
		my $seq=$line;

		$line=<SEQ>;
		$line=~ s/\r//g;
		chomp $line;
		my $qHead=$line;

		$line=<SEQ>;
		$line=~ s/\r//g;
		chomp $line;
		my $qual=$line;

		return ($seqDesc, $seq, $qHead, $qual);
	}
	else{ die "ERROR: Script Borked! Get Sunit (sunitj [ AT ] umich [ DOT ] edu)\n"; }

}

sub parseFasta{
	my 	$line=shift;
	my($seqDesc,@sequence)=split(/\n/, $line);
	my $seq = join ("", @sequence);
	$seq=~ s/\r//g;
	return($seqDesc, $seq);
}

sub exclude{
	my ($name, $np, $seq, $qName, $qual)=@_;
	my @nameParts=@{$np};
	my $c=0;

	foreach my $n(@nameParts){
		$n=uc($n);
		if ($index{$n}){
			$c++;
			$seen{$n}++;
			last;
		} 
		else{
			$c=0;
		}
	}
	unless ($c>0){
		if ($printFastaOut){
			&printFasta($name, $seq);
		}
		else{
			&printFastq($name, $seq, $qName, $qual);
		}
		return 1;
	}
	else { return 0; } 
}

sub include{
	my ($name, $np, $seq, $qName, $qual)=@_;
	my @nameParts=@{$np};
	my $c=0;

	foreach my $n(@nameParts){
		$n=uc($n);
		if ($index{$n}){
			$c++;
			$seen{$n}++;
			last;
		} 
		else{
			$c=0;
		}
	}
	if ($c>0){
		if ($printFastaOut){
			&printFasta($name, $seq);
		}
		else{
			&printFastq($name, $seq, $qName, $qual);
		}
		return 1;
	}
	else { return 0; } 
}

sub printFasta{
	my ($name, $seq)=@_;
	print OUT ">".$name."\n";
	print OUT $seq."\n";
}

sub printFastq{
	my ($name, $seq, $qName, $qual)=@_;
	print OUT "@".$name."\n";
	print OUT $seq."\n";
	print OUT "$qName\n$qual\n";
}

