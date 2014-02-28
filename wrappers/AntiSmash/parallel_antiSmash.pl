#!/usr/bin/perl

=head1 DESCRIPTION

	parallel_antiSmash.pl -- Takes a multifasta file and runs antismash, single fasta at a time in an embarrassingly parallel manner.

=head2 Dependencies

	antiSmash v 2.0.2 (http://www.secondarymetabolites.org/download.html)
	*nix environment.
	Tested on RHEL6.

=head1 USAGE

	perl parallel_antiSmash.pl -fasta MultiSequence.fasta

=head2 Options

	-fasta <STRING>	Multi Fasta file.
	-proc	<INTEGER>	Number of processors. default=30
	-ppr	<INTEGER>	Processors per run. default=2
	-len	<INTEGER>	Only run antismash on sequences greater than this value; default=1000
	-pfam	<BOOLEAN>	perform pfam analysis as well
	-pfamdir	<CHAR>	folder to pfamA database provided by Antismash.
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

	Sunit Jain, (Thu Oct 31 10:21:15 EDT 2013)
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

my $help;
my $version="parallel_antiSmash.pl\tv0.1.1";
my ($fasta, $pfam);
my $procs=30;
my $ppr=2;
my $minLen=1000;
my $pfamdir="/opt/packages/AntiSmash/antismash/generic_modules/fullhmmer/";

GetOptions(
	'f|fasta:s'=>\$fasta,
	'p|procs|proc:i'=>\$procs,
	'l|len:i'=>\$minLen,
	'ppr|c:i'=>\$ppr,
	'pfam'=>\$pfam,
	'pfamdir:s'=>\$pfamdir,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

if (! $fasta){
	print STDERR "[FATAL] Missing required input(s)\n";
	system("perldoc $0 \| cat");
	exit;
}

my (%PIDs, %files, $lastStretch);
open(FASTA, "<".$fasta) || die $!;
$/=">";
my $i=0;
my $running=0;
my $found=0;
my $errorSeqs=0;
while(my $line=<FASTA>){
	chomp $line;
	next unless $line;
	
	my ($name, @seqs)=split(/\n/, $line);
	my $seq=join("", @seqs);
	my $seqLen=length($seq);
	if ($seqLen < $minLen){
		warn "[WARNING] Skipping $name\tLength: $seqLen\n";
		next;
	}

	my $FH;
	my $singleFasta=$i.".fasta";
	# Create new file with single sequence
	open($FH, ">".$singleFasta) || die $!;
	print $FH ">".$line."\n";
	close $FH;

	# Run antismash on new file
	my $addPFam;
	if($pfam){
		$addPFam="--full-hmmer --pfamdir $pfamdir ";
		die "# [FATAL] Could not find directory $pfamdir\n# Please check to see if it exists or that you have proper permissions\n" if (! -d $pfamdir);
	}
	my $pid= &run("run_antismash -c $ppr --statusfile $i.status --input-type nucl --smcogs --clusterblast ".($pfam ? $addPFam : "" )."--outputfolder $i $singleFasta");
	print "[$pid] Starting $singleFasta...\n";
	$PIDs{$pid}++;
	$files{$pid}=$singleFasta;
	$running++;

	# Wait for a antismash job to finish before starting another one.
	if ($running >= ($procs/$ppr)){
		&REAP;
	}
	$i++;
}
close FASTA;
$/="\n";

$lastStretch++;
&REAP;

print "# Potential Secondary Metabolite(s) Found: $found\n";
print "# Number of seqiences that generated errors: $errorSeqs\n";

exit 0;

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
					&cleanUp($pid,1);
				}
				else{
					delete $PIDs{$pid};
					$numPIDs=scalar(keys %PIDs);
					print "Waiting for ".$numPIDs." more processes...\n";
					$running--;
					print "[$pid] Completed $files{$pid}\n";
					&cleanUp($pid);
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

sub cleanUp{
	my $pid=shift;
	my $force=shift;
	my $fasta=$files{$pid};
	my $folder=$fasta;
	$folder=~ s#\.fasta##;
	my $status=$folder.".status";	
	my $clusterFile=$folder."/geneclusters.txt";
	my $smcogFile=$folder."/smcogs/smcogs.txt";
	my $resultPresent=0;

	my $fileFound=`find $folder -name "ctg*"`;
	my @results;
	push(@results, $clusterFile, $smcogFile);
	foreach my $r(@results){
		chomp $r;
		next unless $r;

		if (-z $r){
			next;
		}
		else{
			print ">".$r."\n";
			$resultPresent++;
		}
	}

	if ($force){
		system("rm -r $folder");
		system("cat $fasta >> errorSeqs.fasta");
		$errorSeqs++;
	}
	elsif ($resultPresent > 0){
		print ">Found something in $fasta\n";
		$found++;
		
		## Write to Overview Files ##
		`echo ">$folder" >> Overview.geneclusters.txt`;
		`cat $clusterFile >> Overview.geneclusters.txt`;
		
		`echo "#$folder" >> Overview.smcogs.txt`;
		`cat $smcogFile >> Overview.smcogs.txt`;
		
		copy($fasta, $folder) || die "Failed to copy $fasta: $!\n";
		copy($status, $folder) || die "Failed to copy $status: $!\n";
	}
	else{
		system("rm -r $folder");
	}

	unlink $fasta;
	unlink $status;
	delete $files{$pid};
	return;
}
