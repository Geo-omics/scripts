#!/usr/bin/perl

=head1 DESCRIPTION

	parallel_getGenomesFromTaxa.pl -- Do this.

=head1 USAGE

	perl parallel_getGenomesFromTaxa.pl

=head2 Options


	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Fri Jan 17 09:43:35 EST 2014)
	sunitj [AT] umich [DOT] edu

=head1 License

	This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

	This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use POSIX ":sys_wait_h"; # qw(:signal_h :errno_h :sys_wait_h);
use File::Copy;


my $scripts="/geomicro/data1/COMMON/scripts/SeqTools/";
my $dmp_path="/geomicro/data1/COMMON/publicDB/taxa_dump/";
my $ban="banned.taxa";
my $col=1;
my $procs=30;
my $db="nt";
my ($gi, $list, $outputDB, $resultDir, $giDir);
my $help;
my $version="parallel_getGenomesFromTaxa.pl\tv0.7.1";
GetOptions(
	'ban:s'=>\$ban,
	'gi:s'=>\$gi,
	'giDir:s'=>\$giDir,
	'c|col:i'=>\$col,
	'procs:i'=>\$procs,
	'l|list:s'=>\$list,
	'db|database:s' => \$db,
	'o|out:s'=>\$outputDB,
	'outdir:s'=>\$resultDir,
	'dump:s'=>\$dmp_path,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";
$col--;

if ((!$list) && (! $gi) && (! $giDir)){ system('perldoc', $0); exit;}

chomp($db);
my ($dbType, $ref, $giList, $file, $tmpExt);
if (lc($db) eq "nr"){
	$file="gi_taxid_prot.dmp";
	$ref=File::Spec->catfile($dmp_path, $file);
}
elsif (lc($db) eq "nt"){
	$file="gi_taxid_nucl.dmp";
	$ref=File::Spec->catfile($dmp_path, $file);
}
else{
	warn "[ERROR]: Invalid Database Type. Enter either \"nr \(for Proteins\)\" OR \"nt \(for Nucleotides\)\".\nYou Entered ".$db."\n";
	print "The Script will now exit.\n";
	exit;
}

if (-e $ref){
	print "Taxa Dump file found: $ref\n"; 
}
else{
	die "Could not locate Taxa Dump file here: $ref\n";
}

if (-d $giDir){
	if (! $outputDB){ $giDir.".DB.fasta"; }
	$resultDir=$giDir;
	goto GIDIR
}
elsif(($gi)||($list)){
	if(! $outputDB){
		$outputDB=($gi ? $gi : $list).".DB.fasta";
	}
}

my %bannedTaxa;
if( -e $ban){
	open(BAN, "<".$ban)||die $!;
	while(my $line=<BAN>){
		chomp $line;
		next unless $line;
		next if ($line=~ /^#/);
		
		$bannedTaxa{$line}++;	
	}
	close BAN;	
}

my %idList;
if($gi){
	my %giList;
	open(GI, "<".$gi)|| die $!;
	while (my $line=<GI>){
		chomp $line;
		my(@stuff)=split(/\t/, $line);
		$giList{$stuff[$col]}++;
	}
	close GI;
	my $num=0;
	open(MYREF, "<".$ref) || die $!;
	while(my $line=<MYREF>){
		chomp $line;
		next unless $line;
	
		my ($gi, $taxa)=split(/\t/, $line);
		
		next if ($bannedTaxa{$taxa});
		if ($giList{$gi}){
			$num++;
			$idList{$taxa}++;
			delete $giList{$gi};
		}
		last if (scalar(keys %giList)==0);
	}
	close MYREF;
	print $num." Hits found!\n";
}
elsif($list){
	open (LIST, $list) || die $!;
	while (my $line=<LIST>){
		my(@stuff)=split(/\t/, $line);
		chomp($stuff[$col]);
		$idList{$stuff[$col]}++;
	}
	close LIST;
}
print "Read ".scalar(keys %idList)." unique taxon IDs\n";

if(-d $resultDir){
	warn "[WARNING]".$resultDir." already exists!\n";
	$resultDir.="_".$$;
	print "Output Directory name changed to ".$resultDir."\n";
}
mkdir $resultDir;

print "Aggregating list of GI numbers...\n";
open(REF, "<".$ref)|| die $!;
while(my $line=<REF>){
	chomp $line;
	next unless $line;
	
	my ($gi, $taxa)=split(/\t/, $line);

	next if ($bannedTaxa{$taxa});
	next unless $idList{$taxa};

	my $giList=$ref=File::Spec->catfile($resultDir, $taxa.".gi");
	my $FH;
	open($FH, ">>".$giList)|| die $!;
	print $FH $gi."\n";
	close $FH;
}
close REF;

GIDIR:
print "Getting Fasta files for each taxonomy...\n";
my $i=0;
my $running=0;
my $found=0;
my $errorSeqs=0;
my (%PIDs, %files, $lastStretch);
foreach my $file(glob("$resultDir/*.gi")){
	die "Could not open $file\n" if (! -e $file);
	
	print $file."\n";
	my $tmpFasta=$file;
	$tmpFasta=~ s/\.gi$/\.tmp\.fa/;
	print $tmpFasta."\n";
	
	# Run blastdbcmd on new file
	my $pid= &run("blastdbcmd -db ".lc($db)." -entry_batch $file -out $tmpFasta");
	
	print "[$pid] Getting $tmpFasta...\n";
	$PIDs{$pid}++;
	$files{$pid}=$tmpFasta;
	$running++;

	# Wait for a antismash job to finish before starting another one.
	if ($running >= ($procs-1)){
		&REAP;
	}
	$i++;
}
close FASTA;
$/="\n";

$lastStretch++;
&REAP;

print "Curating new DB...\n";
my $tmpDB=$resultDir."/".($gi ? $gi : $list)."_tmp.DB";
system("cat $resultDir/*.tmp.fa > $tmpDB");
my $script_curateDB=File::Spec->catfile($scripts, "curateDB.pl");
system("perl $script_curateDB -i $tmpDB -o $outputDB -g -d");


sub run{ # create a child process.
	my $command=shift;
	my $pid = fork();

	if (!defined($pid)) {
    	die "unable to fork: $!";
	}
	elsif ($pid==0) { # child
#		print "Executing:\t$command\n"; # The command being executed.
		exec($command) || die "unable to exec: [$?]\n$!\n";
		exit(0);
	}
	# parent continues here, pid of child is in $pid
	return($pid);
}

sub REAP{ ## Use this when you want to wait till a process ends before further processing.
	my $numPIDs= scalar(keys %PIDs);

#	print "in REAPER: ".$numPIDs."\n";
	while ($numPIDs > 0){
		my $pid= waitpid(-1, &WNOHANG);
		if ($pid > 0){
			print "[$pid] Checking ...\n";
			if (WIFEXITED($?) && $PIDs{$pid}){
				if($? > 0){ # There were errors.
					print "[$pid] Crapped out! with status:".WEXITSTATUS($?)."\n";
					delete $PIDs{$pid};
					$numPIDs=scalar(keys %PIDs);
					print "Waiting for ".$numPIDs." more processes...\n";
					$running--;
					#&cleanUp($pid,1);
				}
				else{
					delete $PIDs{$pid};
					$numPIDs=scalar(keys %PIDs);
					print "Waiting for ".$numPIDs." more processes...\n";
					$running--;
					print "[$pid] Completed $files{$pid}\n";
					#&cleanUp($pid);
				}
				return unless $lastStretch;
			}
			else{
				print "[$pid] Waiting...\n";
			}
		}
		else{
			sleep 1;
		}
	}
	return;
}
