#!/usr/bin/perl

=head1 DESCRIPTION

	Consolidate the data obtained from JGI into one tabular file

=head2 Dependencies

	-scripts	<STRING>	location of the script dependencies;	default=/geomicro/data1/COMMON/scripts/SeqTools
	
	extractSubSeq.pl - To extract genes from the contigs
	length+gc.pl - To claculate the length and %GC of the contigs

=head1 USAGE

	perl consolidateJGIdata.pl -DIR path_to_the_JGI_files -OUTDIR output_directory

=head2 Options

	-DIR	-d	<STRING>	path to the files downloaded from JGI;	default=present working directory
	-OUTDIR	-o	<STRING>	Directory with consolidated tab delimited files for each bin;		default=processID
	-genes	-g	<STRING>	produces a fasta file of all the genes in the given file;	default=not produced
	-scripts	<STRING>	location of the extractSubSeq.pl script;	default=/geomicro/data1/COMMON/scripts/SeqTools

	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message. press "q" to exit this screen.

=head3 NOTE:

	In order to get a seperate file for each bin, add a 3rd column in the '*map.txt' file to reflect the bin that the original scaffold belongs to.
	
=head3 Example:

	perl consolidateJGIdata.pl -d ~/JGI_data/ -o output_directory -g output_genes.fasta
	
=head1 Author

	Sunit Jain, (Fri Jun  7 17:53:04 EDT 2013)
	sunitj [AT] umich [DOT] edu

=cut

use strict;
use Getopt::Long;
use File::Spec;
use POSIX ":sys_wait_h";

my $version="0.1.3";
my $DIR="./";
my $outDir=$$;
my $scripts;
my $fasta;
GetOptions(
	'd|DIR:s'=>\$DIR,
	'o|OUTDIR:s'=>\$outDir,
	'g|genes:s'=>\$fasta,
	'scripts:s'=>\$scripts,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system('perldoc', $0); exit;},
);

print "# consolidateJGIdata.pl v$version\n";
# Get all File names in the given directory
unless ($DIR=~ m/\/$/){$DIR=$DIR."/";}
my @FILES=<$DIR*>;
die "[ERROR] Can't find \"$DIR\"\nPlease check that the path exist or that you have sufficient privilages.\n" if (scalar(@FILES)==0);

if (-d $outDir){
	die "[ERROR: $0]Directory '$outDir' already exists!\n";
}
else{
	mkdir($outDir, 0755);
}

my ($cog, $ec, $faa, $fna, $geneProd, $gff, $ko, $contigMap, $pfam, $phyloDist, $config);
foreach my $f(@FILES){
	 $cog=$f if ($f=~ /.*.cog$/i);
	 $ec=$f if ($f=~ /.*.ec$/i);
	 $faa=$f if ($f=~ /.*.faa/);
	 $fna=$f if ($f=~ /.*.fna/);
	 $geneProd=$f if ($f=~ /.*.product*/i);
	 $gff=$f if ($f=~ /.*.gff/);
	 $ko=$f if ($f=~ /.*.ko$/i);
	 $contigMap=$f if ($f=~ /.*.map*/i);
	 $pfam=$f if ($f=~ /.*.pfam*/i);
	 $phyloDist=$f if ($f=~ /.*.phylodist*/i);
}

if ((! $scripts) || (! -d $scripts)){
	if (-d "/geomicro/data1/COMMON/scripts/SeqTools"){
		$scripts="/geomicro/data1/COMMON/scripts/SeqTools";	
	}
	elsif(-e "length+GC.pl"){
		$scripts=`pwd`;
		chomp $scripts;
	}
	else{
		die "[ERROR: $0] Could not locate helper scripts: 'length+GC.pl' and 'extractSubSeq.pl', please provide the location using '-scripts' flag\n";
	}
}
# Run this bit on a seperate thread.
## Do this for all bins seprately
my %PIDs;
my $ess=File::Spec->catfile($scripts, "extractSubSeq.pl");
if ($fasta && -e $ess){
	$fasta=File::Spec->catfile($outDir, $fasta);
	my $pid=run("perl ".$ess." -f ".$fna." -gff ".$gff." -o ".$fasta);
	$PIDs{$pid}++;
}
####

# Continue with the main script.
print "[LGC] Calculating length and GC content of all contigs...\n";
my $lgc=File::Spec->catfile($scripts, "length+GC.pl");
my $tmpLGC=File::Spec->catfile($outDir, $outDir.".lgc");
if(-e $lgc){
	system("perl $lgc -f $fna -gc > $tmpLGC");
}

# Aggregate data from different files.
## Locus Info ##
my %COGS;
if(-e $cog){
print "[COG] File found:\t".$cog."\n";
open(COG, $cog) || die "[COG] $cog :\t$!";
while(my $line=<COG>){
	chomp $line;
	my(@cogData)=split(/\t/, $line);
	$COGS{$cogData[0]}=$cogData[1]."\t".$cogData[2]."\t"; #LocusID =  cog_id <TAB> %id
}
close COG;
}
else{
	warn "[COG] Couldn't find the \`cog\' file\n";
}

my %PFAM;
if(-e $pfam){
print "[PFAM] File found:\t".$pfam."\n";
open(PFAM, $pfam) || die "[PFAM] $pfam :\t$!";
while(my $line=<PFAM>){
	chomp $line;
	my(@pfamData)=split(/\t/, $line);
	$PFAM{$pfamData[0]}=$pfamData[1]."\t"; # LocusID = pfam_id
}
close PFAM;
}
else{
	warn "[PFAM] Couldn't find the \`pfam\' file\n";
}

my %TAXA;
if(-e $phyloDist){
print "[TAXA] File found:\t".$phyloDist."\n";
open(TAXA, $phyloDist); # || die "[PhyloDist] $phyloDist :\t$!";
while(my $line=<TAXA>){
	chomp $line;
	my(@taxaData)=split(/\t/, $line);
	my $locusID=shift @taxaData;
	$TAXA{$locusID}= join("\t", @taxaData);# LocusID = homolog_gene_id <TAB> homolog_taxon_id <TAB> %ID <TAB> Lineage
	$TAXA{$locusID}.="\t";
}
close TAXA;
}
else{
	warn "[TAXA] Couldn't find the \`phylodist\' file\n";
}

my %KO;
if(-e $ko){
print "[KO] File found:\t".$ko."\n";
open(KO, $ko) || die "[KO] $ko :\t$!";
while(my $line=<KO>){
	chomp $line;
	my(@koData)=split(/\t/, $line);
	$KO{$koData[0]}=$koData[2]."\t".$koData[3]."\t"; # LocusID =  ko_term <TAB> %id
}
close KO;
}
else{
	warn "[KO] Couldn't find the \`ko\' file\n";
}

my %EC;
if (-e $ec){
print "[EC] File found:\t".$ec."\n";
open(EC, $ec) || die "[EC] $ec :\t$!";
while(my $line=<EC>){
	chomp $line;
	my(@ecData)=split(/\t/, $line);
	$EC{$ecData[0]}=$ecData[2]."\t"; #LocusID =  EC
}
close EC;
}
else{
	warn "[EC] Couldn't find the \`ec\' file\n";
}

my %PROD;
if (-r $geneProd){
print "[PROD] File found:\t".$geneProd."\n";
open(PROD, $geneProd); # || die "[GENE_PROD] $geneProd :\t$!";
while(my $line=<PROD>){
	chomp $line;
	my(@prodData)=split(/\t/, $line);
	$PROD{$prodData[0]}=$prodData[1]."\t".$prodData[2]."\t"; # LocusID = product <TAB> Source
}
close PROD;
}
else{
	warn "[PROD] Couldn't find the \`gene product\' file\n";
}

my (%LGC);
if(-e $tmpLGC){
#print "[LGC] File found:\t".$tmpLGC."\n";
	open(LGC, $tmpLGC) || die "[LGC] $tmpLGC :\t$!";
	while(my $line=<LGC>){
		chomp $line;
		my @lgcData=split(/\t/, $line);
		$LGC{$lgcData[0]}{"GC"}=$lgcData[1]; # ContigID = %GC
		$LGC{$lgcData[0]}{"Length"}=$lgcData[2];	# ContigID = Length
	}
	close LGC;
	unlink $tmpLGC;
}
else{
	warn "[LGC] Couldn't find the \`$tmpLGC\' file. Check if you have write permissions for this directory\n";
}

## Contig Info ##
my (%contig_name_map, %bins, $totalContigs);
if (-e $contigMap){
print "[MAP] File found:\t".$contigMap."\n";
	open(MAP, $contigMap) || die "[MAP] $contigMap :\t$!\n";
	while(my $line=<MAP>){
		chomp $line;
		my ($original, @data)=split(/\t/, $line);
		$contig_name_map{$data[0]}{"Name"}=$original;
		my $bin;
		if($data[1]){
			$bin=lc($data[1]);
		}
		else{
			$bin="Unclassified";
		}
		
		$contig_name_map{$data[0]}{"Bin"}=$bin;
		$bins{$bin}{"Contigs"}++;
		$bins{$bin}{"GC"}+=$LGC{$data[0]}{"GC"};
		$bins{$bin}{"Bases"}+=$LGC{$data[0]}{"Length"};
		$totalContigs++;
	}
	close MAP;
}
else{
	warn "[MAP] Couldn't find the \`map\' file\n";
}

if (-e $contigMap){
	foreach(keys %bins){
		my $out=File::Spec->catfile($outDir, $_.".tsv");
		open(OUT, ">".$out) || die "[OUT] $out :\t$!\n";
		print OUT "# Bin\tIMG_Contig_Name\tOriginal_Contig_Name\tContig \%GC\tContig Length\tLocus_Tag\tIMG_Gene_ID\tGene_Type\tGene_Start\tGene_Stop\tGene_Length\tHomolog_Gene_ID\tHomolog_Taxon_ID\tLineage \%ID\tLineage\tProduct\tSource\tCOG_ID\tCog \%ID\tPFAM_ID\tKO_Term\tKO \%ID\tEC_Number\n";
		close OUT;
	}
}
else{
	my $out=File::Spec->catfile($outDir, "Unclassified.tsv");
	open(OUT, ">".$out) || die "[OUT] $out :\t$!\n";
	print OUT "# Bin\tIMG_Contig_Name\tOriginal_Contig_Name\tContig \%GC\tContig Length\tLocus_Tag\tIMG_Gene_ID\tGene_Type\tGene_Start\tGene_Stop\tGene_Length\tHomolog_Gene_ID\tHomolog_Taxon_ID\tLineage \%ID\tLineage\tProduct\tSource\tCOG_ID\tCog \%ID\tPFAM_ID\tKO_Term\tKO \%ID\tEC_Number\n";
	close OUT;
}

my (%genes_per_bin, $totalGenes,%printBin, %numCDS);
open(GFF, $gff) || die "[GFF] $gff :\t$!\n";
while(my $line=<GFF>){
	chomp $line;
	my($contigID, $locusID, $geneID, $start, $stop, $type)=parseGFF3($line);
	my($begin, $end)=sort{$a <=> $b}($start, $stop);
	my $len=$end - $begin;
	my $bin=$contig_name_map{$contigID}{"Bin"};
	if ($type=~ /CDS/i){$numCDS{$bin}++}
	my $printThis=$bin."\t".$contigID."\t".$contig_name_map{$contigID}{"Name"}."\t"; # Bin <TAB> IMG Contig Name <TAB> Original Contig Name
	$printThis.=$LGC{$contigID} ? $LGC{$contigID}{"GC"}."\t".$LGC{$contigID}{"Length"}."\t" : "\t\t"; # Contig %GC <TAB> Contig Length
	$printThis.=$locusID."\t".$geneID."\t".$type."\t".$start."\t".$stop."\t"; # Locus_Tag <TAB> IMG_Gene_ID <TAB> Gene_Type <TAB> Gene_Start <TAB> Gene_Stop
	$printThis.=$len."\t"; # Gene_Length
	$printThis.=$TAXA{$locusID} ? $TAXA{$locusID} : "\t\t\t\t"; # homolog_gene_id <TAB> homolog_taxon_id <TAB> %ID <TAB> Lineage
	$printThis.=$PROD{$locusID} ? $PROD{$locusID} : "\t\t"; # product <TAB> Source
	$printThis.=$COGS{$locusID} ? $COGS{$locusID} : "\t\t"; # cog_id <TAB> %id
	$printThis.=$PFAM{$locusID} ? $PFAM{$locusID} : "\t"; # pfam_id
	$printThis.=$KO{$locusID} ? $KO{$locusID} : "\t\t"; # ko_term <TAB> %id
	$printThis.=$EC{$locusID} ? $EC{$locusID} : "\t"; # EC
	
	$printThis=~ s/\t$//;


	$printBin{$bin}.=$printThis."\n";

	$genes_per_bin{$bin}++;
	$totalGenes++;
}
close GFF;
undef %contig_name_map;
undef %COGS;
undef %PFAM;
undef %TAXA;
undef %KO;
undef %EC;
undef %PROD;

foreach my $bin(keys %printBin){
	my $out=File::Spec->catfile($outDir, $bin.".tsv");
	open(OUT, ">>".$out) || die "[OUT] $out :\t$!\n";
	print OUT $printBin{$bin};
	close OUT;
}
undef %printBin;

print "\n# Bin\tNumContigs\tFraction_of_Dataset\tNumGenes\tFraction_of_Dataset\tTotal_Length_of_Bin\tAverage_GC\tNumber_of_CDS\n";

foreach(keys %bins){
	print $_."\t".$bins{$_}{"Contigs"}."\t".($bins{$_}{"Contigs"}/$totalContigs)."\t".$genes_per_bin{$_}."\t".($genes_per_bin{$_}/$totalGenes)."\t".$bins{$_}{"Bases"}."\t".($bins{$_}{"GC"}/$bins{$_}{"Contigs"})."\t".$numCDS{$_}."\n";
}

if ($fasta && -e $ess){
	&REAP;
}

exit;

sub parseGFF3{
# http://gmod.org/wiki/GFF
	my $line=shift;
	my @cols=split(/\t/, $line);
	
	my(@attributes)=split(/\;/, $cols[-1]);
	
	my ($locusID, $geneID);
	my $contigID=$cols[0];
	my $type=$cols[2];
	my $start=$cols[3];
	my $stop=$cols[4];
	foreach my $att(@attributes){
		$locusID= $1 if ($att=~/locus_tag\=(.*)/);
		$geneID= $1 if ($att=~/^ID\=(.*)/);
	}
	if (! $locusID){
		foreach my $att(@attributes){
			if ($att=~/Parent\=(.*)/){
				$locusID=$1."_exon"
			}
			elsif($cols[2]=~/repeat/){
				$locusID=$1 if ($att=~/rpt_type\=(.*)/); # rpt_type=CRISPR;rpt_unit=13023..13055
				$locusID.="[ ".$1." ]" if ($att=~/rpt_unit\=(.*)/);
			}
			else{
				$locusID=$geneID."_".$cols[2];
			}
		}
	}
	return ($contigID, $locusID, $geneID, $start, $stop, $type);
}

sub run{
	my $command=shift;
	my $pid = fork();

	if (!defined($pid)) {
    	die "unable to fork: $!";
	}
	elsif ($pid==0) { # child
		print "Executing:\t$command\n";
		exec($command) || die "unable to exec: [$?]\n$!\n";
		exit(0);
	}
	# parent continues here, pid of child is in $pid
	return($pid);
}

sub REAP{ ## Use this when you want to wait till the process ends before further processing.
	my $numPIDs= scalar(keys %PIDs);

#	print "in REAPER: ".$numPIDs."\n";
	while (scalar(keys %PIDs) > 0){
		my $pid= waitpid(-1, &WNOHANG);
		if ($pid > 0){
#			print "in REAPER:$pid\n";
			print "Process: ".$pid."\tStatus: ".($? == 0 ? "Finished" : "Running")."\nWaiting for ".$numPIDs." more processes...\n";
			if (WIFEXITED($?) && $PIDs{$pid}){
				`echo "Process ID: $pid\tFinished with status $?"`;
#				$numPIDs-- ;
				print "Process: ".$pid."\tStatus: ".($? == 0 ? "Finished" : "Running")."\nWaiting for ".$numPIDs." more processes...\n";
				delete $PIDs{$pid};
			}
		}
		else{
			sleep 10;
		}
	}
	return;
}

