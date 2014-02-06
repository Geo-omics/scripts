#! usr/bin/perl

=head1 NAME

	extractGenomes.pl - extract genomes from NCBI databases (nr/nt) by Taxa curate and concatenate them to form your owm customized database.

=head2 PREPARATION.

	In your terminal window, goto the directory where you've stored your exported file from IMG and type: 'less <your file name>'. 
	This should show you a header line. Start at '1' and count to the column number for 'NCBI_Taxon_ID'. 
	Note the column number. such that, columnA=1, ColumnB=2..and so on.

		eg: less imgCyanoTaxa.tsv (press 'q' to exit)
		
		lets say, my 'NCBI_Taxon_ID' was in column E so my column number is '5'.
	In the same terminal, now type:
		'module load blast/2.2.24'
	
=head2 USAGE

	perl extractGenomes.pl [-i File from IMG] [-o Your Database Name] [-c NCBI Taxon ID Column Number] [-outdir database directory]
				
=head2 OPTIONS

	OPT        Meaning                Expected Values		Default Values.
	-r        redundant Seq File		(y/n)				n
	-d        database			(nr/nt)			nr(proteins)
	-j        concatenate			(y/n)				y
	-f        Individual Genome Files	(y/n)				n
	-c        Column Number		NCBI Taxon Col Num			5
	-outdir   Output Directory		DIR name		<yourUsername>DB
	--help    Help (This page).

=head2 EXAMPLE.

	perl extractDB.pl -i imgCyanoTaxa.tsv -o imgCyanoDB -c 5 -d nt .

=head2 OUTPUT.

	You can find all the individual files in a directory named '<your_user_name>DB'.
	eg. If I were to run this script the folder name would be 'sunitjDB'.

	This folder/directory contains your fasta/faa files with taxon IDs as names and a Concatenated Database
	file  (unless you turned the concatenation off) with name you provided for your database.
	By default, this concatenated file will be formatted.

	Your database is now ready to use.

=head2 AUTHOR

	Sunit Jain: sunitj_at_umich_dot_edu

=cut

#CommandLine
#blastdbcmd -db nr -entry all -outfmt "%g %T" | awk '{if ($2 == <DESIRED_TAX_ID1> || <DESIRED_TAX_ID2> || ...) {print $1 }} ' | blastdbcmd -db nr -entry_batch - -out output_file_name.txt


use strict;
use Getopt::Long;
use File::Spec;
use Benchmark;

# start timer;
my $start=new Benchmark;
my $dateXYZ=`date`;
my @timeXYZ=split(/\s/, $dateXYZ);
my $timestampXYZ=join("_", @timeXYZ);
$timestampXYZ=~ s/:/_/g;
chomp ($timestampXYZ);
## Default Values ##
my $uName = `whoami`;
chomp ($uName);
my $dmp_path="/geomicro/data1/COMMON/publicDB/taxa_dump/";
my $listFile;
my $db='nr';
my $outputDB=$$.'DB.cat';
my $resultDB=$$.'DB';
my $cNum=5;
my $header='y';
my $redun='n';
my $concat='y';
my $faa='n';
my @taxids;

my $version="extractGenomes.pl\tv0.0.9";
GetOptions(
	'i|in:s'=>\$listFile,
	'd|database:s' => \$db,
	'o|out:s'=>\$outputDB,
	'c|col:i'=>\$cNum,
	'header:s'=>\$header,
	'r|redundant:s'=>\$redun,
	'j|join:s'=>\$concat,
	'f|faa:s'=>\$faa,
	'outdir:s'=>\$resultDB,
	't|taxid:s'=>\@taxids,
	'dump:s'=>\$dmp_path,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=> sub {system('perldoc', $0); exit;},
);
print "# $version\n";
if ((!$cNum) && (scalar(@taxids)==0)){ system('perldoc', $0); exit;}


### MAIN ###

my %idList;
if (@taxids){
	undef($cNum);
	chomp(@taxids);
	foreach my $t(@taxids){
		$idList{$t}++;
	}
}

if ($cNum){
	if (!$listFile){ system('perldoc', $0); exit;}
	open (LIST, $listFile) || die "[err] $listFile:\n$!\n";
	my $lineNum=0;
	while (my $line1=<LIST>){
		$lineNum++;
		if (($lineNum == 1) && (lc($header) eq 'y')){ next;}
		else{	
			my(@stuff)=split(/\t/, $line1);
			chomp($stuff[$cNum-1]);
			$idList{$stuff[$cNum-1]}++;
		}
	}
	close LIST;
}

chomp($db);
my ($fType, $ref, $giList, $file, $tmpExt);
if (lc($db) eq "nr"){
	$file="gi_taxid_prot.dmp";
	$ref=File::Spec->catfile($dmp_path, $file);
	$fType='T';
	$tmpExt='fap';
}
elsif (lc($db) eq "nt"){
	$file="gi_taxid_nucl.dmp";
	$ref=File::Spec->catfile($dmp_path, $file);
	$fType='F';
	$tmpExt='fan';
}
else{
	warn "[ERROR]: Invalid Database Type. Enter either \"nr \(for Proteins\)\" OR \"nt\(for Nucleotides\)\".\nYou Entered ".$db."\n";
	print "The Script will now exit.\n";
	exit;
}

if (-e $ref){
	print "Taxa Dump file found: $ref\n"; 
}
else{
	die "Could not locate Taxa Dump file here: $ref\n";
}

my @allGenomes;
mkdir $resultDB unless (-d $resultDB);

foreach my $k(keys %idList){
	print "Extracting Taxon ID: $k\n";
	my $giList=$k."_".$$.".gi";
	open (TAXGI, ">".$giList);
	open (INDEX, $ref) || die "[err] $ref:\n$!\n";

	while(my $line=<INDEX>){
		my($gi, $taxa)=split(/\t/, $line);
		chomp($gi);
		chomp($taxa);
		chomp($k);
		print TAXGI $gi."\n" if ($taxa == $k);
	}
	close INDEX;
	
	my $customDB=$k.".".$tmpExt;
	chomp($customDB);
	system("blastdbcmd -db ".$db." -entry_batch ".$giList." -out ".$customDB);
#	unlink($giList);
		
	print "Checking for Duplicates in the DB...\n";
	print "\tNumber of Sequences before Curation: ";
	system ("grep -c \"\^\>\" ".$customDB);
#	print "\n";
#	my $curated=curate($customDB, $db, $redun, $tmpExt);
	print "Done!\n";
	
	print "Removing temporary files...\n\n";
	system(" mv ".$customDB." ./".$resultDB."/");
	system(" mv reduSeqs_".$customDB.".red ./".$resultDB."/") unless (lc($redun) eq 'n');
	unlink ($customDB);
	push(@allGenomes, $customDB);
}


if (lc($concat) eq 'y'){
	chdir $resultDB;
	my $allG=join(" ", @allGenomes);
	system ('cat '.$allG.' > '. $outputDB);
#	system ('formatdb -i '.$outputDB.' -p '.$fType);
	system ('rm *.'.$tmpExt) if (lc($faa) eq 'n');
}

my $stop=new Benchmark;
my $time=timediff($stop, $start);
print "For ".@allGenomes." Genomes the script took: ". timestr($time, 'all')."\n";
print STDERR "File uncurated.\n";
exit 0;
