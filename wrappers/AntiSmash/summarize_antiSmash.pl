#!/usr/bin/perl

=head1 DESCRIPTION

	summarize_antiSmash.pl -- Do this.

=head1 USAGE

	perl summarize_antiSmash.pl -dir project_folder

=head2 Options

	-dir	<CHAR>	path to the project folder with all the Antismash outputs
	-prefix	<CHAR>	output file prefix
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Thu Nov 21 16:17:53 EST 2013)
	sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Find;

my $help;
my $version="summarize_antiSmash.pl\tv0.0.2b";
my($DIR, $prefix);
GetOptions(
	'dir:s'=>\$DIR,
	'p|prefix:s'=>\$prefix,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

# Get all File names in the given directory
if ($DIR!~ m/\/$/){
	$DIR=$DIR."/";
}

if(!$prefix){
	$prefix= $$;
}

# Get output names.
my $GeneClustSummary=$prefix.".gclust.summary";
my $GeneClustSeqs=$prefix.".gclust.names";
my $smcogSummary=$prefix.".smcog.summary";
my $smcogSeqs=$prefix.".smcog.names";

my @FILES=<$DIR*>;
die "[ERROR] Can't find \"$DIR\"\nPlease check that the path exist or that you have sufficient privilages.\n" if (scalar(@FILES)==0);

my @content;
find(\&wanted, $DIR);

print "# Files found:\t";
my $numberOfFiles=0;
my(@clusterBlastFiles, $Overview_geneCluster, $Overview_smcogs, @nrpspks);
foreach my $file(@content){
	if($file=~ m#clusterblast/cluster#i){
		push(@clusterBlastFiles,$file);
		$numberOfFiles++;
#		print $file."\n";
	}
	elsif($file=~ m#Overview.geneclusters#i){
		$Overview_geneCluster=$file;
		$numberOfFiles++;
#		print $file."\n";
	}
	elsif($file=~ m#Overview.smcogs.txt#i){
		$Overview_smcogs=$file;
		$numberOfFiles++;
#		print $file."\n";
	}
	elsif($file=~ m#nrpspks_predictions_txt/(.*fasta)#){
		push(@nrpspks, $file);
		$numberOfFiles++;
#		print $file."\n";
	}
}
print $numberOfFiles."\n";

print "Reading Gene Cluster Overview file...\t";
open(OGC, "<".$Overview_geneCluster) || die $!;
my (%GCS, %num_name);
my $num="";
while(my $line=<OGC>){
	chomp $line;
	next unless $line;
	if ($line=~ /^>/){
		$line=~ s/^>//;
		$num=$line;
		next;
	}

#	my($folder, $gcInfo)=split("\n", $line);
	my ($tag, $name, $designation, @data)=split(/\s+/, $line);
	push(@{$GCS{$designation}}, $name);
	$num_name{$num}=$name;
}
close OGC;
print "done\n";

print "Summarizing Gene Cluster Overview file...\t";
open(GCS, ">".$GeneClustSummary) || die $!;
open(GCSeq, ">".$GeneClustSeqs) || die $!;
foreach my $d(keys %GCS){
	print GCS $d."\t".scalar(@{$GCS{$d}})."\n";
	print GCSeq $d;
	foreach my $n(@{$GCS{$d}}){
		print GCSeq "\t".$n;
	}
	print GCSeq "\n";
}
close GCS;
close GCSeq;
print "done\n";

print "Reading smCOG Overview file...\t";
open(OSM, "<".$Overview_smcogs)||die $!;
my $currentSeqNum=0;
my %SMCOGS;
while(my $line=<OSM>){
	chomp $line;
	next unless $line;
	
	if ($line=~ /^#/){
		$line=~ s/^#//;
		$currentSeqNum=$line;
	}
	else{
		if ($line=~ "^>>"){
			$line=<OSM>;
			$line=<OSM>;
			$line=<OSM>;
			my($smcog, @data)=split(/\t/, $line);
			push(@{$SMCOGS{$smcog}}, $currentSeqNum);
		}
	}
}
close OSM;
print "done\n";

print "Summarizing smCOG Overview file...\t";
open(SS, ">".$smcogSummary) || die $!;
open(SSeq, ">".$smcogSeqs) || die $!;
foreach my $s(keys %SMCOGS){
	my $smcog_id=$s;
	$smcog_id=~ s/:/\t/;
	print SS $smcog_id."\t".scalar(@{$SMCOGS{$s}})."\n";
	print SSeq $smcog_id;
	foreach my $n(@{$SMCOGS{$s}}){
		print SSeq "\t".$num_name{$n};
	}
	print SSeq "\n";
}
close SS;
close SSeq;
print "done\n";

exit 0;
sub wanted{
	push @content, $File::Find::name;
	return;
}

__END__
# summarize cluster Blast data
my (%CBs, %descKey);
print "# ".scalar(@clusterBlastFiles)." Cluster Blast files found\n";
foreach my $cb(@clusterBlastFiles){
	chomp $cb;
	next unless $cb;
	
	my $CBFH;
	open($CBFH, "<".$cb)|| die $!;
	while(my $line=<$CBFH>){
		chomp $line;
		if($line=~ /^Type/){

			my($header, $description)=split(/\s+/, $line);
			$CBs{$description}{$cb}++;
		}
	}
	close $CBFH;
}

open(CB, ">".$clusterSummary)||die $!;
print "Cluster Blast Summary:\n";
print CB "Types";
foreach my $f(@clusterBlastFiles){
	print CB "\t".$f;
}
print CB "\tTotal\n";

foreach my $desc(keys %CBs){
	print CB $desc;
	my $total;
	foreach my $file(@clusterBlastFiles){
		if($CBs{$desc}{$file}){
			print CB "\t".$CBs{$desc}{$file};
			$total+=$CBs{$desc}{$file};
		}
		else{
			print CB "\t0";		
		}
	}
	print CB "\t$total\n";
}
close CB;
