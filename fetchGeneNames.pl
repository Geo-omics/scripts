#!/usr/local/bin/perl

use encoding 'utf8';

use XML::LibXML;
use XML::LibXML::NodeList;
use LWP::Simple;
use URI::Escape;

# Search Parameters
my $term = "diabetes";
my $max_results = "3";

# String Constants
my $eutils="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
my $myemail = uri_escape ("biocat\@umich.edu");
my $mytool = uri_escape ("genebot");

# Initialize Fetch Counter
my $fetchcount = 0;

# Output Header
print "\n************ RESULTS FROM QUERY \"$term\" ************\n";

# Construct Disease Query
my $q_ids = $eutils."esearch.fcgi?db=gene&term=".uri_escape ($term).
        "&email=$myemail&tool=$mytool&retmode=xml".
        "&retstart=$start&retmax=$max_results";

# Query NCBI - returns results in XML format
my $xml_ids = my_get ($q_ids);
# print "xml_ids is $xml_ids \n\n";

# Parse XML to obtain ids
my @id_list = parse_ids_from_xml ($xml_ids);
# print "id_list is @id_list \n";

# For each id, query NCBI to retrieve the
# corresponding aliases
foreach (@id_list) {

  # Construct Specific Gene Query
  my $q_summary = $eutils."esummary.fcgi?db=gene&id=".$_.
        "&email=$myemail&tool=$mytool&retmode=xml";

  my $xml_summary = my_get ($q_summary);
  # print "xml_summary is $xml_summary \n\n";

  my ($geneName, $aliases) = parse_aliases_from_xml ($xml_summary);

#  my $aliases = parse_aliases_from_xml ($xml_summary);
  print "GENE NAME: $geneName \n";
  print "ALIASES: $aliases \n\n";

}

#########################################################
#                  HELPER FUNCTIONS                     #
#########################################################

# Make a request to the NCBI server
sub my_get {
  my ($u) = @_;
  for (my $i=0; $i<5; $i++) {
    $fetchcount++;
    # print "fetchcount is: $fetchcount \n";
    if ($fetchcount > 3) { sleep 3; $fetchcount = 0; }
    my $raw = get ($u);
    if (defined $raw && $raw ne "") { return $raw; }
    sleep 3;
  }
  return "";
}

# Helper function to retrieve ids from XML
sub parse_ids_from_xml {

  my ($x_ids) = @_;

  # Set up XML Parser
  my $parser = XML::LibXML->new ();
  my $doc = $parser->parse_string ($x_ids);

  # Get Search Results
  my $parent = $doc->findnodes ('eSearchResult');

  # Traverse Search Results to Acquire IDs
  foreach my $node ($parent->get_nodelist) {

    my $results = $node->findnodes('IdList/Id');
    foreach my $id ($results->get_nodelist) {
      my $n = $id->textContent;
      push @ids, $n;
    }
    last;
  }

  # Return array of IDs
  return @ids;
}

# Helper function to retrieve aliases from XML
sub parse_aliases_from_xml {

  my ($x_summary) = @_;

  # Set up XML Parser
  my $parser = XML::LibXML->new ();
  my $doc = $parser->parse_string ($x_summary);

  # Look for all Item Nodes
  my $item_nodes = $doc->findnodes('eSummaryResult/DocSum/Item');

  # Initialize loop variables
  $n = "";
  $a = "";

  # Find Item Nodes with Name attribute OtherAliases
  foreach my $node ($item_nodes->get_nodelist) {
    my $item_name = $node->getAttribute("Name");

    # Get Gene Name
    if ($item_name eq "Name") {
      $n = $node->textContent;
    }

    # Get Aliases
    if ($item_name eq "OtherAliases") {
      $a = $node->textContent;
    }
  }

  return ($n, $a);
}
