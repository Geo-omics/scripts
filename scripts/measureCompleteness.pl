#!/usr/bin/perl

=head1 DESCRIPTION

	measureCompleteness.pl -- Given a fasta file and IMG consolidated data look for the 36 essential housekeeping genes (Cicareli et al 2006)
				NOTE1: Just matches COG numbers.
				NOTE2: Users can provide their own list as well. See "myCOGS" option for details.

=head1 USAGE

	perl measureCompleteness.pl -list bin.list -tsv IMG_consolodated.data -p bin -detail

=head2 Options

	-list	-l	<CHAR>	List of headers in your bin matching the contig or IMG_Name column of consolidated data
	-tsv		<CHAR>	Consolidated IMG data generated from "consolidateJGIdata.pl"
	-prefix	-p	<CHAR>	Output file; Contains --> presence/absence, copy numbers and contig names. (extension= .out)
	
	-img_names	<BOOLEAN>	use IMG assigned names to match headers.
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head2 Additional Features

	-binned		<BOOLEAN>	If consolidated file belongs to the list of scaffolds that you're interested in ONLY, use this instead of "-list" flag.
	-detail	-d	<BOOLEAN>	Optional Output File; print line from consolidated data containing the match. (extension= .tsv)
	-myCOGS	-m	<CHAR>	user defined list of cog numbers. One entry per line; only first 2 columns will be read.
			FORMAT: Column1=COG#### <TAB> COG_Description

=head1 Author

	Sunit Jain, (Thu Mar  6 13:41:57 EST 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Spec;

my %COGS;
my($list, $tsv, $imgName, $col, $more, $out, $myCOGS, $binned);
my $prefix=$$;
my $help;
my $version="measureCompleteness.pl\tv0.2.0";
GetOptions(
	'l|list:s'=>\$list,
	'binned'=>\$binned,
	'p|o|prefix:s'=>\$prefix,
	'tsv:s'=>\$tsv,
	'img_names'=>\$imgName,
	'd|detail'=>\$more,
	'm|myCOGS:s'=>\$myCOGS,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

## Sanity Check ##
die &help if (! $tsv);

###################
## Load Defaults ##
###################
if($more){
	$more=$prefix.".tsv";
}
$out=$prefix.".out";
if($imgName){ $col=1; } else { $col=2; }

## Load COGS ##
if ($myCOGS){
	open(COGS, "<".$myCOGS)|| die $!;
	while(my $line=<COGS>){
		next if ($line=~ /^#/);
		chomp $line;
		next unless $line;
		
		my($cog, @d)=split(/\t/, $line);
		my $desc=join(" ", @d);
		$desc="Description not provided" if (scalar(@d) == 0);
		print "COG:$cog\tDESC:$desc\n";
		$COGS{uc($cog)}=$desc;
		exit;
	}
	close COGS;
}
else{
	&loadDefaultCOGS;
}

## Load Contig names of Interest ##
my %index;

unless($binned){
	open(LIST, "<".$list ) || die $!;
	while(my $line=<LIST>){
		next if ($line=~ /^#/);
		chomp $line;
		next unless $line;
		
		$line=~ s/^>//;
		$index{lc($line)}++;
	}
	close LIST;
}
print "# Searching for [ ".scalar(keys %COGS)." ] housekeeping genes in [ ".scalar(keys %index)." ] contigs\n";

## Lookup data in the consolidated file ##
my (%found, @details, $header);
open(TSV, "<".$tsv ) || die $!;
while(my $line=<TSV>){
	if ($line=~ /^#/){
		$header=$line;
	}
	chomp $line;
	next unless $line;
	
	my(@data)=split(/\t/, $line);
	my $name=$data[$col];
	my $cog=uc($data[17]);

	if(!$binned){
		next unless $index{lc($name)};
	}
	next unless $COGS{$cog};
	
	push (@{$found{$cog}}, $name);
	
	if($more){
		push(@details, $line);
	}
}
close TSV;

## Let 'em have it ##
if (scalar(keys %found) > 0){
	print "# Found [ ".scalar(keys %found)." ] out of a possible [ ".scalar(keys %COGS)." ] housekeeping genes, at least once.\n";
	my $measure=sprintf("%.2f", (scalar(keys %found)/scalar(keys %COGS))*100);
	print "# Estimated completeness [ ".$measure." ]\n";
	&writeOutput;
}
else{
	print "# Could not find a single housekeeping gene!\n";
}

##################
## Sub-Routines ##
##################
sub writeOutput{
	open(OUT, ">".$out) || die $!;
	print OUT "# COG_Num\tCOG_Description\tCopy_Number\tContig_Name(s) -->\n";
	foreach my $cog(keys %found){
		print  OUT $cog."\t".$COGS{$cog}."\t".scalar(@{$found{$cog}});
		foreach my $contig(@{$found{$cog}}){
			print  OUT "\t".$contig
		}
		print  OUT "\n";
	}
	close OUT;
	
	if($more){
		open(MORE, ">".$more) || die $!;
		print MORE $header;
		foreach my $line(@details){
			print MORE $line."\n";
		}
		close MORE;
	}
}

sub help{
	system("perldoc $0 \| cat");
	exit 1;
}

### DO NOT MODIFY ###
sub loadDefaultCOGS{
%COGS=(
"COG0080"=>"L11",
"COG0081"=>"L1",
"COG0087"=>"L3",
"COG0091"=>"L22",
"COG0093"=>"L14",
"COG0094"=>"L5",
"COG0097"=>"L6P/L9E",
"COG0102"=>"L13",
"COG0197"=>"L16/L10E",
"COG0200"=>"L15",
"COG0256"=>"L18",
"COG0048"=>"S12",
"COG0049"=>"S7",
"COG0052"=>"S2",
"COG0092"=>"S3",
"COG0096"=>"S8",
"COG0098"=>"S5",
"COG0099"=>"S13",
"COG0100"=>"S11",
"COG0103"=>"S9",
"COG0184"=>"S15P/S13E",
"COG0186"=>"S17",
"COG0522"=>"S4",
"COG0016"=>"Phenylalanyl-tRNA synthethase alpha subunit",
"COG0018"=>"Arginyl-tRNA synthetase",
"COG0060"=>"Isoleucyl-tRNA synthetase",
"COG0124"=>"Histidyl-tRNA synthetase",
"COG0143"=>"Methionyl-tRNA synthetase",
"COG0172"=>"Seryl-tRNA synthetase",
"COG0201"=>"Preprotein translocase subunit SecY",
"COG0495"=>"Leucyl-tRNA synthetase",
"COG0525"=>"Valyl-tRNA synthetase",
"COG0202"=>"DNA-directed RNA polymerase, alpha subunit/40 kD subunit",
"COG0085"=>"DNA-directed RNA polymerase, beta subunit/140 kD subunit",
"COG0012"=>"Predicted GTPase",
"COG0533"=>"Metal-dependent protease"
);
}
