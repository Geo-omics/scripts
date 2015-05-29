#! /usr/bin/perl

=head1 NAME

	tetramer_freqs_esom.pl 
		calculates Tetramer frequencies for the given fasta file and produces 4 output files

=head1 USAGE
	
	perl tetramer_freqs_esom.pl -f fastaFile -a annotationFile [OPTIONS]

=head1 OPTIONS
	
	-f		Required	Fasta File, may include X:s and N:s
	-a		Required	Annotation File (3 columns; 1. full contig name, 2.annotation, 3.Class Number (Your metagenome has class#0, everything else 1+))
	-k		Optional	default=4;	any reasonable k-mer size; ideally don't change this but if you do keep it between 3-6.
	-min	Optional	default=2500; Minimal length (in nt) of input contig to be included in output
	-max	Optional	default=5000
	Note:	The script will split sequence after each 'max' nt; join last part, if remaining seq shorter than 'max', with second-last part
			eg: in default settings, a sequence of 14 kb will be split into a 5 kb and a 9 kb fragment if window_size = 5 kb.

=head1 AUTHORS

	# Anders Andersson, 2007 (anders.andersson@scilifelab.se)
	
=head2 Modified by

	# Itai Sharon, Nov/2010
	# Sunit Jain, 2010 (sunitj [AT] umich [DOT] edu)
		
=head1 CITATION

	If you use this script please cite:
	Dick, G.J., A. Andersson, B.J. Baker, S.S. Simmons, B.C. Thomas, A.P. Yelton, and J.F. Banfield (2009).
	Community-wide analysis of microbial genome sequence signatures. Genome Biology, 10: R85.
	
=head1 LICENSE & COPYRIGHT

	Copyright (C) 2007 Anders Andersson (anders.andersson@scilifelab.se)
	This is free software; see the COPYRIGHT file accompanying this script
	for details. There is NO warranty; not even for MERCHANTABILITY or 
	FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;

my $version=$0." v2.0.4 October, 2013";
my $sfile; #fasta file, may include X:s and N:s
my $annotationfile; #full contig name in left, annotation in right, column. headers (whatever) on first line 
my $min_length = 2500; #Minimal length (in nt) of input contig to be included in output
my $window_size = 5000; #split sequence after each window_size nt, 
                     #join last part, if shorter than window_size, 
                     #with second-last part (a sequence of 
                     #14 kb will be split into a 5 kb and a
                     #9 kb fragment if window_size = 5 kb)
my $ext="fasta";
my $kmer_size 		= 4;

GetOptions(
	'f|fasta=s'=> \$sfile,
	'a|ann=s'=>\$annotationfile,
	'min:i'=>\$min_length,
	'max:i'=>\$window_size,
	'k|kmer:i'=>\$kmer_size,
	'v|version'=> sub{&licensing; exit},
	'h|help'=> sub{system('perldoc', $0); exit;},
	'e|ext:s'=>\$ext,
);

&licensing;

if ((! $sfile) || (! $annotationfile)){print "[ERROR $0] Missing required input.\nFor help using the script, type 'perl $0 -help'\n"; exit;}

print "\n#################### ".$kmer_size."-mer Frequencies #####################\n";
print "Minimum length (in bases) of input contig to be included in output:\n";
print "$min_length\n";
print "Window Size:\n$window_size\n";
my $seqfile = basename($sfile, ".".$ext);
#print "This is the fasta file:$seqfile\n";
#!!! Program will automatically create outfiles called (whatever) infile.lrn and infile.names (and overwrite existing) !!!

my %allowed = ();	# Allowed k-mers
my @mers = ();
my @names = ();
my @tetras = ();
my @seq_list = ();

my $lrnfile = "Tetra_".$seqfile."_$min_length\.lrn";
my $namesfile = "Tetra_".$seqfile."_$min_length\.names";
my $classfile = "Tetra_".$seqfile."_$min_length\.cls";
my $reffile= "Tetra_".$seqfile."_".$min_length."_".$window_size."_split.fasta";
my $n=0;
#$window_size = 5000;

####### main #############
&make_list_of_possible_tetramers('', $kmer_size);
&calc_tetra_freqs;
&make_lrn_file;
&make_names_file;
&make_class_file;
&make_seq_file;
&getRowColESOM;

print  "\nAll Done!!\n";
exit 0;

##### sub routines #######
# This will create the list of all possible k-mers, for any reasonable k-mer size.
sub make_list_of_possible_tetramers {
	my ($mer, $k) = @_;

	if($k == 0) {
		my $rc_mer = make_revcomp($mer);
		if (!defined $allowed{$rc_mer}) {
			push (@mers, $mer);
			$allowed{$mer}++;
		}
		return;
	}

	my @bases = ("A", "T", "C", "G");
	foreach my $na (@bases) {
		make_list_of_possible_tetramers("$mer$na", $k-1);
	}
}
################################################################################################################
sub calc_tetra_freqs {
	print  "calculating ".$kmer_size."-mer frequencies ";
	my $total_index = 0;
	my ($id, $seq) = (undef, "");

	open (INFILE, $sfile) || die "can't open $sfile!";
	my $counter = 0;
	while (<INFILE>) {
    	chomp;
    	if ($_ =~ />(\S+)/) {
    	$counter++;
		print  '.' if($counter % 100000 == 0);
		my $next_id = $1;
		get_tetra_freqs($id, $seq) if (length($seq) >= $min_length);
		($id, $seq) = ($next_id, '');
		
    	} else {
        		$seq .= $_;
    	}
	}

	# Last sequence
    	if (length($seq) >= $min_length) {
        	get_tetra_freqs($id, $seq);
    	}
    	close (INFILE);
	print  "ok, $counter data points\n";
}

################################################################################################################
my $total_index = 0;
sub get_tetra_freqs {
	my ($id, $seq) = @_;

	$seq = uc($seq);

	# filter out short sequences between N's and X's 
	# as well as between these and beginning and end of sequence
	my @lowqual = ();
	push(@lowqual, 0); #to get start position
	foreach my $i (1 ..  (length($seq) - $kmer_size - 1)) {
    	my $base = substr($seq, $i, 1);
	# The following covers all ambiguous characters allowed in the FASTA format (check http://www.boekhoff.info/?pid=data&dat=fasta-codes)
    	if ($base =~ /KMRYSWBVHDNX/) {
        		push(@lowqual, $i);
    	}
	}
	push(@lowqual, length($seq)); #to get end position
	my $filtered_seq = $seq;
	foreach my $i (1 .. $#lowqual) {
    	my $length = $lowqual[$i] - $lowqual[$i - 1] + 1;
    	if ($length < 50) {
    	    	for (my $j = $lowqual[$i - 1]; $j < $lowqual[$i]; $j++) {
    	        	substr($filtered_seq, $j, 1) = "Z";
    	    	}
    	}
	}
	$seq = $filtered_seq;

	#split sequence into subsequences
	my @sub_seq = ();
	if (length($seq) < 2*$window_size) {
    	@sub_seq = ($seq);
	} else {
    	for(my $i=0; $i<length($seq); $i = $i + $window_size) {
        		my $subseq = substr($seq, $i, $window_size);
        		push(@sub_seq, $subseq);
    	}
    	if (length($sub_seq[-1]) < $window_size) {
        		$sub_seq[-2] = $sub_seq[-2].$sub_seq[-1];
        		pop (@sub_seq);
    	}
	}

	#calculate and print freqs for each subsequence
	my $sub_index = 0;
	foreach $seq (@sub_seq) {
    	$sub_index++;
    	my %this_mers = ();
    	my $sum = 0;
    	foreach my $i (0 .. (length($seq)-1)) {
    		my $mer = substr($seq, $i, $kmer_size);
    		if (defined $allowed{$mer}) {
        		$this_mers{ $mer }++;
        		$sum++;
    		} else {
        		my $rc_mer = &make_revcomp($mer);
        		if (defined $allowed{$rc_mer}) {
            			$this_mers{ $rc_mer }++;
            			$sum++;
        		}
    		}
    	}
		my $tetra = "";
		$total_index++;
		my $name = "$total_index\t$id"."_"."$sub_index\t$id";
		push(@names, $name);
		push(@seq_list, $seq);
		foreach my $mer (@mers) {
	    	if (defined $this_mers{$mer}) {
	        		my $counts = $this_mers{$mer}/$sum;
	            	$tetra = $tetra."\t".$counts;
	    	} else {
	        		$tetra = $tetra."\t0";
	    	}
		}            
		push(@tetras, $tetra);
	}
}
################################################################################################################
sub make_revcomp {
	my $rc = $_[0];
	$rc =~ tr/ACGT/TGCA/;
	return reverse($rc);
}
################################################################################################################
sub make_lrn_file {
	print  "printing lrn file $lrnfile ... ";
	open (OUT, ">$lrnfile") || "can't create outfile $lrnfile";
	my $number_rows = @names;
	my $number_cols = @mers + 1;
	print OUT "% $number_rows\n";
	print OUT "% $number_cols\n";
	print OUT "% 9";
	foreach my $mer (@mers) {
    	print OUT "\t1";

	}
	print OUT "\n";
	print OUT "% Key";
	foreach my $mer (@mers) {
    	print OUT "\t$mer";
	}
	print OUT "\n";
	my $key = 0;
	foreach my $tetra (@tetras) {
    	$key++;
    	print OUT "$key$tetra\n";
	}
	close (OUT);
	print  "ok\n";
}
################################################################################################################
sub make_names_file {
 	print  "printing names file $namesfile ... ";
	my $number_rows = @names;
	open (OUT, ">$namesfile");
	print OUT "% $number_rows\n";
	foreach my $name (@names) {
		print OUT "$name\n";
    	}
    	close (OUT);
	print  "ok\n";
}
################################################################################################################
sub make_class_file {
    	print  "printing class file $classfile ... ";
	my %class = ();
	my $line = 0;

	open (INFILE, $annotationfile) || die "can't open $annotationfile";
	while (<INFILE>) {
    	$line++;
		chomp;
		# contig	annotation
		my @fields = split(/\t/, $_);
		$fields[0] =~ s/\s+$//;
		$class{$fields[0]} = $fields[2];
	}
   	close (INFILE);

	open (OUT, ">$classfile");
   	my $number_rows = @names;
	print OUT "% $number_rows\n";
	foreach my $item (@names) {
    	my @fields = split(/\t/, $item);
    	print OUT "$fields[0]\t$class{$fields[2]}\n";
	}
	close (OUT);
	print  "ok\n";
}
################################################################################################################
sub make_seq_file {
    print "printing seq file: $reffile ... ";
    my $number_rows = @names;
	my $getseq= @seq_list;
	open (OUT, ">$reffile");
#    print OUT "% $number_rows\n";

	foreach my $item (@names) {
        my @fields = split(/\t/, $item);
        print OUT ">$fields[1]\n";
		print OUT "$seq_list[$n]\n";
		$n++;
	}
    close (OUT);
	print  "ok\n";
}

################################################################################################################
sub getRowColESOM{
# "The product of rows and columns, i.e. the number of neurons should be at least 1K neurons. The ratio of rows and columns should be significantly different from unity." -- http://databionic-esom.sourceforge.net/user.html#Preprocessing
	my $acceptedSeq=@names;

	my $mapSpace= $acceptedSeq * 5.5;
	my $rows= int(sqrt($mapSpace/2) + 0.5);
	my $cols= 2 * $rows;

	print "\nTry the following values for ESOM Training:\n>Rows:\t$rows\n>Cols:\t$cols\n";
	print "These values are just meant as suggestions, feel free to try your own.\n";
}

sub licensing{
	print $version."\n\n";
	
	print "Copyright (C) 2007 Anders Andersson (anders.andersson\@scilifelab.se)\nThis is free software; see the COPYRIGHT file accompanying this script for details. There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n";

	print "Please cite:\nDick, G.J., A. Andersson, B.J. Baker, S.S. Simmons, B.C. Thomas, A.P. Yelton, and J.F. Banfield (2009).  Community-wide analysis of microbial genome sequence signatures. Genome Biology, 10: R85.\n\n";
}
