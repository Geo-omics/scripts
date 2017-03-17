#!/user/bin/perl
use strict;
use Bio::SearchIO; 
use Getopt::Long;

=head1 Description

	parse BlastXML results and print out Query, Subj, %id, evalue, bitScore. To get more contact me.

=head2 Usage

	perl parseBlastXML.pl -b <blastXML file>

=head3 Optional

	-o output file name; default processID.tsv,
	-s minimum bit score; default 0

=head4 NOTE:

	To use this script Bioperl MUST be installed.
	For our lab, here is how you can access Bioperl for this script:
		1) Make sure you're on Cayman.
		2) type "/opt/package/perl5/5.12.2/bin/perl" instead of just "perl" while running the script. example:
			/opt/package/perl5/5.12.2/bin/perl parseBlastXML.pl -b <blastXML file>

=head2 Author

	Sunit Jain, July 2011
	sunitj [AT] umich [DOT] edu

=cut

my $blastXml;
my $tab=$$.".tsv";
my $bs=0;
GetOptions(
	'b:s'=>\$blastXml,
	'o:s'=>\$tab,
	's:f'=>\$bs,
);

open (OUT, ">".$tab);
print OUT "#File:$blastXml\tMin_BitScore:$bs\n";
print OUT "#Query\tHit\t\%ID\tEvalue\tBitScore\n";
chomp($blastXml);

my $in = new Bio::SearchIO(-format => 'blastxml', -file=> $blastXml);
while( my $result = $in->next_result ) {
  ## $result is a Bio::Search::Result::ResultI compliant object
	while( my $hit = $result->next_hit ) {
		## $hit is a Bio::Search::Hit::HitI compliant object
		while( my $hsp = $hit->next_hsp ) {
		## $hsp is a Bio::Search::HSP::HSPI compliant object
			if($hsp->score >= $bs){
				print OUT $result->query_name,
				"\t", $hit->name,
				"\t", $hsp->percent_identity,
				"\t", $hsp->evalue,
				"\t", $hsp->score,"\n";
			}
		}
	}
}
