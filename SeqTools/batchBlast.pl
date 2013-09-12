#! /usr/bin/perl
=head1 MOTIVATION

	The script began as a wrapper script to run multiple instances of blast jobs when given a list of query files.
	It has since evolved into a wrapper for catalouging all blast jobs run to avoid duplicating results and wasting resources.
	It handles all basic functionalities of blast jobs (both new and old), and more can be added upon request.

=head2 Dependencies

	Blast 2.2.20 (old); OR
	Blast 2.2.22 or above (new)

=head1 USAGE

	Start multiple instances of blast at once:
		batchblast.pl [-p blast type] [-l list of query files]
	OR
	Start a single instance of blast:
		batchblast.pl [-p blast type] [-q single query fasta file]
	OR
	Check if any type of blast was performed on the query file:
		batchBlast.pl [-check] [-q query fasta file]

=head2 Options

	-p OR -blast:	Blast Task (blastn, blastp, blastx, tblastn, tblastx)
	-l OR -list:	List of Query File Names for Blasts; Each name should be in a seperate line.
	-q OR -query:	Single query fasta file.
	-m OR -outfmt:	Output Format. Default: outputs w/ comments.
	-n:		Makes the script use the latest version of blast. Make sure you have loaded the appropriate module.
			Default: uses old blast syntax.
	-a OR -num_threads:	Multi-thread your blast; default: 1
	-e OR -evalue:	Evalue; Default: 1e-3
	-d OR -db:	database; Default: 'nr' OR 'nt', depending on '-p'
	-per	percentage identity; only available with the newer versions of blast (2.2.22 and above).
	-check		:	Check if blasts for this file have been performed earlier.
	-h OR -help:	this page.

=head1 AUTHOR

	Sunit Jain, Aug 2011
	sunitj [AT] umich [DOT] edu

=cut

#######################
## MODULES
#######################
use strict;
use Getopt::Long;
use File::Basename;

#######################
## PARAMETERS
#######################
my ($newBlast,$p,$d,$m,$check);
my $listOfFiles="";
my $a=1;
my $e="1e-5";
my $per=0;
my $w=11;
my $query="";
my $myID=`whoami`;
chomp($myID);
my $mainLogUpdate="test.log";

GetOptions(
	'n|new'=>\$newBlast,
	'p|blast:s'=>\$p,
	'l|list:s'=>\$listOfFiles,
	'i|query:s'=>\$query,
	'a|num_threads:i'=>\$a,
	'e|eval:s'=>\$e,
	'd|db:s'=>\$d,
	'm|outfmt:i'=>\$m,
	'w|word_size'=>\$w,
	'perc_identity:i'=>\$per,
	'check'=>\$check,
	'h|help'=>sub{system("perldoc", $0); exit;},
);

&checkArgs;
my ($commandLine, @fileNames, %procsToFollow, %blastLogTable);
my @table=("Type","Query","Database","E-value","PercentID","WordSize","Output","ProcessorsUsed","Started","Ended","Owner","PWD","BlastVersion","Log","OutFormat","Query_md5");

#######################
## MAIN
#######################
push (@fileNames, $query);
open (LOF, $listOfFiles) if $listOfFiles;
open (LOG, ">>".$mainLogUpdate)|| die "[error] $mainLogUpdate : $! \n";

@fileNames=<LOF> if $listOfFiles;
foreach my $file(@fileNames){
	chomp($file);
	$blastLogTable{"Started"}=localtime();
	$blastLogTable{"Query"}=$file;
#	my ($queryMD5, $fName)=split(/\s+/, `md5sum $file`);
#	$blastLogTable{"Query_md5"}="$queryMD5";
	$newBlast ? &blastPlus($file) : &blastall($file);
}

&getProcIDs;

#######################
## SUB-ROUTINES
#######################
sub help{system('perldoc', $0); exit;}

sub checkArgs{
	&check if ($check);

	if (!$p){ print "[ERROR]Missing Required Options!\n"; &help}
	if (!$listOfFiles && !$query){print "[ERROR]Missing Query Files Options!\n"; &help}

	if (!$m){ $m = $newBlast ? 6 : 8; }
	
	$p= lc($p);

	my %blasts=(
		"blastn"=>"nt",
		"blastp"=>"nr",
		"blastx"=>"nt",
		"tblastn"=>"nr",
		"tblastx"=>"nr",
	);

	die "[error] -p should be one of 'blastn', 'blastp', 'blastx', tblastn', 'tblastx' \nYou entered: $p\n" unless $blasts{$p};

	if (!$d){
		$d = $blasts{$p};
	}

	if (!$newBlast && $per){ warn "Older version of blast(2.2.20) does not support this function\n"; }

	my $pwd=`pwd`;
	chomp($pwd);

	%blastLogTable=(
		"Type"=>"$p",	
		"BlastVersion"=>"",
		"Query"=>"",
		"Database"=>"$d",
		"WordSize"=>"$w",
		"Output"=>"",
		"Log"=>"",
		"E-value"=>"$e",
		"PercentID"=>"$per",
		"ProcessorsUsed"=>"$a",
		"OutFormat"=>"$m",
		"Started"=>"",
		"Ended"=>"",	
		"Owner"=>"$myID",
		"PWD"=>"$pwd",
		"Query_md5"=>"",
	);
}

sub getOutputName{
	my $query=shift;
	my @q=split(/\//, $query);
	my @parts=split(/\./, $q[-1]);
	my @dbParts=split (/\//, $d);
	my $addDB= join("_", @parts[0..($#parts -1)], "vs", $dbParts[-1]); 
	my $out= join(".", $addDB,$p);
	my $log= join(".", $addDB, "log");
	$blastLogTable{"Output"}="$out";
	$blastLogTable{"Log"}="$log";
	return ($out, $log);
}

sub blastPlus{
	my $query=shift;
	my ($out, $log)=getOutputName($query);
	$commandLine= "$p -query $query -out $out -num_threads $a -evalue $e -db $d -outfmt $m";
	system("nohup $commandLine > $log &");
	$blastLogTable{"BlastVersion"}="Blast+";
}

sub blastall{
	my $query=shift;
	my ($out, $log)=getOutputName($query);
#	system("perl countTo100.pl 2");
	$commandLine="blastall -p $p -i $query -o $out -a $a -e $e -d $d -m $m";
	system("nohup blastall $commandLine > $log &");
	$blastLogTable{"BlastVersion"}="Blastall";
}

sub getProcIDs{
	my @allProcIDs;
	@allProcIDs=`ps H -u $myID`;

	foreach my $p(@allProcIDs){
		my($pid)= &parseProcTable($p);

		next unless ($pid);
		$procsToFollow{$pid}++;
	}
	my $num_jobs=keys(%procsToFollow);
	print "Total: ".$num_jobs."\n";
	&keepAnEyeOn($num_jobs);
}

sub keepAnEyeOn{
	my $num_jobs=shift;
	while ($num_jobs){
		foreach my $proc(keys(%procsToFollow)){
			my @p=`ps $proc`;

			my($pid)= &parseProcTable($p[1]);
			if ($pid){
				next;
			}
			elsif(!$pid){
				print $proc." Finished\n";
				$num_jobs--;
				$blastLogTable{"Ended"}=localtime();
			}
			else{
				exit;
			}
		}
		($num_jobs == 0) ? last : sleep(60);
	}
	return;
}

sub updateLog{
	foreach my $c(@table){
		print LOG $blastLogTable{$c}."\t";
	}
	print LOG "\n";
}

sub parseProcTable{
	my $proc=shift;

	$proc=&trim($proc);
	my($pid, $tty, $stat, $time, @command)= split(/\s+/, $proc);
	
	return unless ($pid=~ /\d+/);

	my $c= join(" ", @command);
	return unless ($c=~ /$commandLine/i);

	return ($pid);
}

sub trim{
	my $y=shift;
	$y=~ s/^\s+//;
	$y=~ s/\s$//;
	return $y;
}

sub check{
	die "[error] Query File required\n" if (!$query);
	open(LOGS, "$mainLogUpdate") || die "[error] Cannot access $mainLogUpdate. See if you have permissions for this file.\n";
	my ($qMd5, $fName)=split(/\s+/, `md5sum $query`);
	my $matches=0;
	while(my $line=<LOGS>){
		print $line if ($line=~ /^#/);
		chomp($line);
		$line=~ s/\r//;

		my @cols=split(/\t/, $line);
		if ($cols[-1] eq $qMd5){
			print $line."\n";
			$matches++;
		}
	}
	($matches > 0) ? print "$matches Matches Found!\n" : print "No Matches Found\n";
	exit;
}
