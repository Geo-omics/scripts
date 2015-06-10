#!/usr/bin/perl -w

=head1 NAME

getLineageNT - classify BLAST hits by taxonomy (Nucleotides only).

=head2 USAGE

perl getLineage.pl [-i tab_file] [-i second_BLAST_file] [-e evalue_cutoff]
                      [-t dir_where_TAXONOMY_files_are] [-g gi2taxid] 
                      [-z PATH_TO_zcat] [-v]


=head2 DESCRIPTION

Will print out the taxonomic distribution (at the kingdom level) for a
set of hits against the NR database.  This script assumes you've done
a search against the protein database, you'll have to make minor
changes in the gi_taxid part to point to the gi_taxid_nuc.dump file.

This expects BLAST files in tabbed -m9 or -m8 format.  Output with -m
8 or use blast2table.pl to convert (or fastam9_to_table.PLS if using
FASTA).

  Input values:
   -t/--taxonomy  directory where the taxonomy .dmp files are (from NCBI)
   -g/--gi        Location of gi_taxid_prot.dmp (or gi_taxid_nucl.dmp if 
                  the search was against a NT db)
   -i/--in        The name of the tab delimited -m8/-m9 output files to 
                  process.

    -e/--evalue   Provide an E-value cutoff for hits to be considered
    -z/--zcat     Path to the 'zcat' executable, can also be 'gunzip -c'
                  if no zcat on your system.
   Flags
    -v/--verbose  To turn on verbose messages
    -h/--help     Display this helpful information

This is intended to be useful starting script, but users may want to
customize the output and parameters.  Note that I am summarizing the
kingdoms here and Eukaryota not falling into Metazoa, Viridiplantae,
or Fungi gets grouped into the general superkingdom Eukaryota. for
simplicity.  There are comments in the code directing you to where
changes can be made if you wanted to display hits by phylum for
example.  Note that you must wipe out the cache file 'gi2class' that
is created in your directory after making these changes.

=head2 AUTHOR

Jason Stajich jason_at_bioperl_dot_org

=head2 EDITED BY

Sunit Jain sunitj_at_umich_dot_edu

=cut

use strict;
use Bio::DB::Taxonomy;
use DB_File;
use Env;
use File::Spec;
use vars qw($SEP);
my $DEBUG = 0;
use Getopt::Long;
$SEP = '_';

# Set relative path to dmp files.
my $home=$ENV{"HOME"};
my $up=File::Spec->updir();
my $d1="COMMON";
my $d2="extractedDB";
my $d3="dump";
my @dir=($home, $up, $d1, $d2, $d3);


# Default Values.
my $evalue_filter = 1e-4;
my @files;
my $zcat = 'zcat'; # or gunzip -c 
my $prefix = File::Spec->catdir(@dir);
my $gi2taxidfile = File::Spec->catfile(@dir, "gi_taxid_prot.dmp");

GetOptions(
	   'v|verbose|debug' => \$DEBUG,
	   'z|zcat:s'    => \$zcat,
	   'i|in:s'      => \@files,
	   'e|evalue:f'  => \$evalue_filter,
	   't|taxonomy:s'=> \$prefix,
	   'g|gi|gi2taxid:s' => \$gi2taxidfile,
	   'h|help'      => sub { system('perldoc', $0);
				  exit },
	   );

# insure idxSP location is created
mkdir(File::Spec->catfile($prefix,'idxSP')) 
    unless -d File::Spec->catfile($prefix,'idxSP');

# these files came from ftp://ftp.ncbi.nih.gov/pub/taxonomy
my $taxdb = Bio::DB::Taxonomy->new
    (-source => 'flatfile',
     -directory => File::Spec->catfile($prefix, 'idxSP'), 
     -nodesfile => File::Spec->catfile($prefix,'nodes.dmp'),
     -namesfile => File::Spec->catfile($prefix,'names.dmp')
     );
my %query;

my $range = 10000000000;
my $randomNumber = int(rand($range));
my $gi2className='spGi2class'.$randomNumber;

my (%taxid4gi,%gi2node);
my $dbh = tie(%gi2node, 'DB_File', $gi2className);

#open (TEST, $prefix."gi_taxid_nucl.dmp") || die $!;

my $giidxfile = File::Spec->catfile($prefix,'idxSP','gi2taxid');
my $done = -e $giidxfile;
my $dbh2 = tie(%taxid4gi, 'DB_File', $giidxfile);

if( ! $done ) {
    my $fh;
    # this file came from ftp://ftp.ncbi.nih.gov/pub/taxonomy
    # I'm interested in protein hits therefor _prot file.
    if( $gi2taxidfile =~ /\.gz$/ ) {
		open($fh, "$zcat $gi2taxidfile |" ) || die "$zcat $gi2taxidfile: $!";
    } else {
		open($fh, $gi2taxidfile ) || die "[err] ".$gi2taxidfile."\n".$!;
    }
    my $i= 0;
    while(<$fh>) {
    next if ($_=~/^\#/);
		my ($gi, $taxid) = split(/\t/,$_);
		chomp($gi);
		chomp($taxid);
		$taxid4gi{$gi} = $taxid;
		$i++;
		if( $DEBUG && ($i % 100000 == 0) ) {
		    warn "$i\n";
		}
    }
    $dbh2->sync;
}

#open (NOTAXA, ">noTaxIDfound_Sp".$randomNumber.".log");
for my $file ( @files ) {
	open (OUT, ">l_$file.txt");
    warn("$file\n");
    my $gz;
    if( $file =~ /\.gz$/) {
	$gz = 1;
    }
    my ($spname) = split(/\./,$file); 
    my ($i, $fh);
    if( $gz ) {
	open($fh, "$zcat $file |")  || die "$zcat $file: $!";
    } else {
	open($fh, $file) || die "$file: $!";
    }
    while(<$fh>) {
    next if ($_=~/^\#/);
	my ($qname,$hname,$pid,$qaln,$mismatch,$gaps,
	    $qstart,$qend,$hstart,$hend,
	    $evalue,$bits,$score) = split(/\t/,$_);	
	#next if( $evalue > $evalue );
	if( ! exists $query{$spname}->{$qname} ) {
	    $query{$spname}->{$qname} = {};
	}

	if( $hname =~ /gi\|(\d+)/) {		
	    my $gi = $1;	    
	    if( ! $gi2node{$gi} ){ # see if we cached the results from before
		my $taxid = $taxid4gi{$gi};
		if( ! $taxid ) {
#		    warn("no taxid for $gi\n");
		    print NOTAXA $gi."\n";
		    print OUT $gi."\tnotFound\($gi\)\n";
		    next;
		}
		my $node = $taxdb->get_Taxonomy_Node($taxid);
		if( ! $node ) {
		    warn("cannot find node for gi=$gi ($hname) (taxid=$taxid)\n");
		    next;
		}
		my $parent = $taxdb->get_Taxonomy_Node($node->parent_id);

		# THIS IS WHERE THE KINGDOM DECISION IS MADE
		# DON'T FORGET TO WIPE OUT YOUR CACHE FILE
		# spGi2class after you make changes here
		while( defined $parent && $parent->node_name ne 'root' ) { 
		    # this is walking up the taxonomy hierarchy
		    # can be a little slow, but works...
#		    warn( "\t",$parent->rank, " ", $parent->node_name, "\n");
		    # deal with Eubacteria, Archea separate from 
		    # Metazoa, Fungi, Viriplantae differently
		    # (everything else Eukaryotic goes in Eukaryota)
		    if($parent->rank eq 'family') {
			# caching in ... 
				($gi2node{$gi}) = $parent->node_name;
				print OUT $gi."\t".$parent->node_name."\t".$parent->rank."\n";
				last;
		    } elsif($parent->rank eq 'order') {
			# caching in ... 
				($gi2node{$gi}) = $parent->node_name;
				print OUT $gi."\t".$parent->node_name."\t".$parent->rank."\n";
				last;
		    }  elsif($parent->rank eq 'class') {
			# caching in ... 
				($gi2node{$gi}) = $parent->node_name;
				print OUT $gi."\t".$parent->node_name."\t".$parent->rank."\n";
				last;
		    }  elsif($parent->rank eq 'phylum'){
			# caching in ... 
				($gi2node{$gi}) = $parent->node_name;
				$gi2node{$gi} =~ s/ \<(bacteria|archaea)\>//g;
				print OUT $gi."\t".$parent->node_name."\t".$parent->rank."\n";
				last;
		    }
		    elsif($parent->rank eq 'kingdom') {
			# caching in ... 
				($gi2node{$gi}) = $parent->node_name;
				$gi2node{$gi} =~ s/ \<(bacteria|archaea)\>//g;
				print OUT $gi."\t".$parent->node_name."\t".$parent->rank."\n";
				last;
		    }
		   	elsif($parent->rank eq 'no rank') {
			# caching in ... 
				($gi2node{$gi}) = $parent->node_name;
				$gi2node{$gi} =~ s/ \<(bacteria|archaea)\>//g;
				print OUT $gi."\t".$parent->node_name."\t".$parent->rank."\n";
				last;
		    }
   		   	elsif($parent->rank eq '') {
			# caching in ... 
				($gi2node{$gi}) = $parent->node_name;
				$gi2node{$gi} =~ s/ \<(bacteria|archaea)\>//g;
				print OUT $gi."\t".$parent->node_name."\t".$parent->rank."\n";
				last;
		    }
		    $parent = $taxdb->get_Taxonomy_Node($parent->parent_id);
		}		
		$dbh->sync;
	    }
	    my ($kingdom) = $gi2node{$gi};
#		warn("$gi2node{$gi}\n");
	    unless( defined $kingdom && length($kingdom) ) {
#		    warn("no kingdom for $hname\n");
	    } else {
		$query{$spname}->{$qname}->{$kingdom}++;		
	    }	
	} else {
	    warn("no GI in $hname\n");
	}
    }
    last if ( $DEBUG && $i++ > 10000);
    close OUT;
}
close NOTAXA;

# print out the taxonomic distribution
while( my ($sp,$d) = each %query ) {
	open (LG, ">".$sp.".log");
    my $total = scalar keys %$d;
    print LG "$sp total=$total\n";
    my %seen;
    for my $v ( values %$d ) {
	my $tag = join(",",sort keys %$v );
	$seen{$tag}++;
    }
    for my $t ( sort { $seen{$a} <=> $seen{$b} } keys %seen ) {
	printf LG " %-20s\t%d\t%.2f%%\n",
	$t,$seen{$t}, 100 * $seen{$t} / $total;
    }
    print LG "\n\n";

}

system('rm -f '.$gi2className);
#system('rm -rf idxSP');
