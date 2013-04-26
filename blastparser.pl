#!/usr/bin/perl -w
use strict;

use Bio::SearchIO;

# usage: blastparser.pl <blastfile>
# by sheri

my $blast_report = $ARGV[0];
my ($query_length, @results);

my $in = new Bio::SearchIO(-format => 'blasttable', 
                           -file   => $blast_report);
my $outfile = $blast_report . ".parsed.csv";
open(OUT, ">$outfile");
print OUT "query_name\thit_name\tquery_length\thit_length\tpct_ident \n";

while( my $result = $in->next_result ) {
  while( my $hit = $result->next_hit ) {
   while( my $hsp = $hit->next_hsp ) {
	$query_length = $hsp->length('query');
    if( $hsp->length('total') >= 0.5*$query_length && $hit->name ne $result->query_name && $hsp->length('total')>=50) { # print out hit if hit length is greater than half query gene length and percent identity exceeds 96
     if ( $hsp->percent_identity >= 96 ) {
      push(@results, [ $result->query_name,  $hit->name, $query_length, $hsp->length('total'), $hsp->percent_identity]);
     }
    }
   }  
  }
}
for my $i (0 .. scalar(@results)-1) {
  print OUT $results[$i][0], "\t", $results[$i][1], "\t", $results[$i][2], "\t", $results[$i][3], "\t", $results[$i][4], "\n" unless ($results[$i][1] eq $results[$i+1][1]);
}
close OUT;
