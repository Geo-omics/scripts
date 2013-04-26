#!/user/bin/perl
=head1 DESCRIPTION

	getCol.pl - extract the desired column number from a tab delimited file.
	
=head2 USAGE

	perl getCol.pl -i <tab delimited file>
	
=head3 Optional

	-c cloumn number to extract; default col 1
		NOTE: to count from the last column start counting from '-1' eg: -1, -2, -3 and so on.
	-o output file name; default processID.list
	-u print unique values only; default print everything.
	-t print unique values with their tallies

=head3 Example

	For UNIQUE values in the column:
		perl getCol.pl -i fileWithManyColumns.txt -o fileWithOneColumnOfUniques.txt -c 2 -u

	For ALL values in the column, in order of appearence.
		perl getCol.pl -i fileWithManyColumns.txt -o fileWithOneColumn.txt -c 2

=head2 AUTHOR

	SUNIT JAIN
	sunitj [AT] umich [DOT] edu	
	
=cut

use strict;
use Getopt::Long;

my $in;
my $col=1;
my $out = $$.".list";
my $u;
my $t;
GetOptions(
	'i|in:s'=>\$in,
	'c|col:i'=>\$col,
	'o|out:s'=>\$out,
	'u|unique'=>\$u,
	't|tally'=>\$t,
	'h|help'=>sub{ system("perldoc", $0); exit;},
);

if (! $in){ system("perldoc", $0); exit;}
my $c;
if ($col > 0){
	$c= $col - 1;
}
else{
	$c=$col;
}
my %seen;
open(FILE, $in) || die "[error] $ARGV[0]:$!\n";
open (OUT, ">".$out);
print OUT "\#Column #".$col."\n";

print "printing unique values with tallies.\n" if ($t);
print "printing unique values only.\n" if ($u && ! $t);
print "printing Column $col in order of appearence.\n" if (! $u && ! $t);

while (my $line=<FILE>){
	next if $line=~ m/^#/;
	chomp $line;
	$line=~ s/\r//g;
	next unless $line;
	my(@row)=split(/\t/, $line);
	print OUT $row[$c]."\n" if (! $u && ! $t);
	$seen{$row[$c]}++;
}

if ($u || $t){
	while (my ($k, $v)=each (%seen)){
		if ($t){
			print OUT $k."\t".$v."\n";
		}
		else{
			print OUT $k."\n";
		}
	}
}
else{
	exit;
}
