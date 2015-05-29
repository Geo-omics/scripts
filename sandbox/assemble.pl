#! /usr/bin/perl

=head1 Description

	assemble.pl was created to be THE one stop shop for your illumina assembly needs. Once you're satified with the quality of your reads, put 'em here and the script will create the assemblies for you. The pipeline also includes support for Oases (for meta/transcriptomics) and MetaVelvet (for meta-genomics). Support for MetaIDBA and AMOS (minimus) is in the works.

=head2 Dependencies

=head3 Software

	Velvet 1.1.07 or above
	Oases 0.2.01 or above
	MetaVelvet or 1.0.01 or above

=head3 Scripts - all the required scripts can be found in the 'SeqTools' directory

	-scripts: change the default directory for script dependencies: default: /geomicro/data1/COMMON/scripts/SeqTools
	
	interleave.pl - to interleave the forward and reverse fastq/fasta files
	limit2length.pl	- to generate general statistics of the final assembly
	calcN50.pl - to generate the N50 and L50 statistics of the final assembly
	findStretchesOfNs.pl - to find the stretches of Ns longer than k-mer
	
=head2 Usage

	Case 1: Assemble files as singletons
	perl assemble.pl -singles file1.fastq file2.fastq file3.fastq -k # (kmer value)

	Case 2: Assemble files as paired; this will interleave the fwd and rev files and then assemble.
	perl assemble.pl -fwd fwd.fq -rev rev.fq -k # (kmer value) -i # (insert size)

	Case 3: Assemble files as paired; when you already have an interleaved file
	perl assemble.pl -paired interleaved.fq -k # (kmer value) -i # (insert size)

	Case 4: Mixed
	perl assemble.pl -paired interleaved.fq -singles file1.fastq file2.fastq file3.fastq -k # (kmer value) -i # (insert size)
	OR
	perl assemble.pl -fwd fwd.fq -rev rev.fq -singles file1.fastq file2.fastq file3.fastq -k # (kmer value) -i # (insert size)

=head2 Options

=head3 Required:
	
	fwd		:	Forward Sequence file
	rev		:	Reverse Sequence file
	paired	:	When you have an interleaved file, only accepts one file at a time, can be used with singles.
	singles	:	When you have multiple files and you wish to force velvet into treating them individually.
				You may give it as many files as long as the names are seperated by a space.
	k or kmer	:	K-mer aka Hash size upto 99; We usually use k-mers of 61, 75, 91 but feel free to play around.
	i or insert	:	Insert length; Required if using paired ended or interleaved data. 
			Look at your bioanalyzer results, let's say you see a peak at 391, that's your insert size.
	-sd		:	insert length standard deviation; default=13

=head3 Optional:

	outdir	:	The directory that contains the output
	p or prefix	:	The prefix for your interleaved output
	interval	:	The script also tracks the CPU and Memory consumption for the assemblies; default=10 (seconds)
	log		:	change the name of the log file for the script
	contig_size	:	Minimum Contig Length in your output; default=200

=head3 Boolean Flags

	v		:	Version
	fasta	:	If your sequence files are in the fasta format.
	debug	:	When I need to debug the script. It doesnt actually run the assemblies, just prints out the commands it would have passed;

=head2 Modifying Assembly Type

	The following options require that the corresponding modules be loaded before executing the script.

=head3 Available:

	-trans	:	for meta/transcriptomic reads.(requires Oases)
	-metav	:	for metavelvet (metagenomic reads only).

=head3 ToDo:

	Assembler Support
	-idba	:	uses IDBA-UD for assembly.

=head1 Comments/Accolades/Brickbats/Beers:
	
	Sunit Jain, 2012
	sunitj AT umich DOT edu
	
	last updated: June 2013

=cut



#######################
## MODULES
#######################
use strict;
use Getopt::Long;
use File::Basename;
use File::Spec;
use POSIX ":sys_wait_h"; # qw(:signal_h :errno_h :sys_wait_h);

#######################
## PARAMETERS
#######################
my($intlv, $fwd, $rev, @singles, $KMER, $INS, $OUTDIR, $transcripts, $trim, $derep, $fasta, $DEBUG, $metaV, $amos, $LOG, $prefix);
my $INS_SD= 13;
my $version= "assemble.pl\tv0.0.10";
my $interval=10;
my $scripts;
my $minLen=1999;
my $min_contig_len=200;
GetOptions(
	'paired:s'=> \$intlv,
	'singles:s{,}'=>\@singles,
	'fwd=s'=> \$fwd,
	'rev=s'=> \$rev,
	'k|kmer=i'=> \$KMER,
	'i|insert=i'=>\$INS,
	'sd|insert_sd:i'=>\$INS_SD,
	'contig_size:i'=>\$minLen,
	'minAssemblyLen:i'=>\$min_contig_len,
	'outdir:s'=>\$OUTDIR,
	'trans'=>\$transcripts,
	'trim:s'=>\$trim,
	'derep:s'=>\$derep,
	'fasta'=>\$fasta,
	'p|prefix|o|out:s'=>\$prefix,
	'debug'=>\$DEBUG,
	'metav'=>\$metaV,
	'amos'=>\$amos,
	'interval'=>\$interval,
	'log:s'=>\$LOG,
	'scripts:s'=>\$scripts,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub {system('perldoc', $0); exit;},
);
#######################
## CHECKS
#######################
print "\# $version\n";
## Check if velvet module loaded ##
my @tmp=`velveth 2>&1`; # velvet/1.1.07-MAX99-OPENMP
&helpLoadingModules if ((scalar(@tmp)) < 2);

my $module_error;
if ($amos){
	my @tmp=`bank-transact 2>&1`; # AMOS/3.1.0
	$module_error++ if ((scalar(@tmp)) < 2);
}
if ($metaV){
	my @tmp=`meta-velvetg 2>&1`; # MetaVelvet/1.0.01
	$module_error++ if ((scalar(@tmp)) < 2);
}
if ($transcripts){
	my @tmp=`oases 2>&1`; # oases/0.2.01
	$module_error++  if ((scalar(@tmp)) < 2);
}

&helpLoadingModules if ($module_error >= 1);

if ($amos || $metaV){
	print "WARNING: The Meta-Velvet and AMOS portion of the Script is still under development. Only the processes leading up to this point will be completed.\n";
}

if (! $KMER){	&help;	}
if (! $fwd && ! @singles && !$intlv){	&help;	}

die "[ERROR: $0] Invalid k-mer! Please enter an odd integer between 11 and 99\n" if (($KMER > 99) or ($KMER < 11));

my $sCount="";
my $seqType= $fasta ? "fasta" : "fastq";
my $beginsWith= $fasta ? "^\>" : "^\@";
my ($intPair, $single, $usage);

my $pair++ if ($fwd && $rev);

if ($fwd && ! $rev){
	die "[ERROR: $0] $!: $fwd\n" unless (-e $fwd);
	warn "No reverse strand found. Treating the sequences as singles\n";
	push (@singles, $fwd);
	$usage="singles";
}
elsif ($pair){
	$usage="paired";
	die "[ERROR: $0] $fwd not found\n" unless (-e $fwd);
	die "[ERROR: $0] $rev not found\n" unless (-e $rev);
}
elsif($intlv){
	$usage="paired";
	$intPair=$intlv;
	die "[ERROR: $0] $intlv not found\n" unless (-e $intlv);
}
elsif(! $fwd && ! $rev && @singles){
	$usage="singles";
}

####
if ((! $scripts) || (! -d $scripts)){
	if (-d "/geomicro/data1/COMMON/scripts/SeqTools"){
		$scripts="/geomicro/data1/COMMON/scripts/SeqTools";	
	}
	elsif((-e "interleave.pl") && ($pair)){
		$scripts=`pwd`;
		chomp $scripts;
	}
	else{
		die "[ERROR: $0] Could not locate helper scripts: 'interleave.pl', 'limit2length.pl', 'findStretchesOfNs.pl' and 'calcN50.pl', please provide the location using '-scripts' flag\n";
	}
}
###

#######################
## GLOBAL
#######################
my $email= &identifyUser;
print "A summary file will be emailed to: $email when the job completes\n";

my ($dir, $suf);
if (! $prefix){
	$prefix="assembly_".$usage;
}
($prefix, $dir, $suf)=fileparse($prefix);

my (%PIDs, %useCase);
$OUTDIR="assembly_".$usage."_".$KMER;
die "$OUTDIR already exists! Choose another output directory\n" if (-d $OUTDIR);
$LOG="$OUTDIR.log";
warn "[WARNING: $0] $LOG will be overwritten!\n" if (-e $LOG);
unlink $LOG if (-e $LOG);
#######################
## MAIN
#######################

open(LOG, ">".$LOG)|| die "[ERROR: $0] $!: $LOG\n";

die "[ERROR $0] Insert length required if using paired ended data\n" if (($usage eq "paired") && ! $INS);
print "Interleaving, this may take a while:\n" if ($pair);
&interleave if ($pair);

&addOptions;

&assemble;		

&REAP;

&getStats;
close LOG;
exit 0;

#######################
## SUB-ROUTINES
#######################
sub addOptions{
	foreach my $s(@singles){
		#system("echo $s\n >>$LOG");
		my $fileSize= -s $s;
		next if $fileSize == 0;
		$single.="-".$seqType." -short".$sCount." ".$s." ";
		$sCount++;
	}
	%useCase=(
		'paired'=>"-".$seqType." -shortPaired ".$intPair." ".$single,
		'singles'=>$single,
	);
}

sub assemble{
	if ($DEBUG){
		if (! $pair && ! $intlv && ($sCount==0)){die "[ERROR: $0] Check Singleton files\n"; }

		print "velveth $OUTDIR $KMER $useCase{$usage} >> $LOG\n";
		if ($usage eq "paired"){
			print "velvetg $OUTDIR -exp_cov auto -ins_length $INS -ins_length_sd $INS_SD -read_trkg yes -amos_file yes -min_contig_lgth $min_contig_len -unused_reads yes >> $LOG\n";
		}
		elsif($usage eq "singles"){
			print "velvetg $OUTDIR -exp_cov auto -read_trkg yes -amos_file yes -min_contig_lgth $min_contig_len -unused_reads yes >> $LOG\n";
		}
		if ($metaV){
			if ($usage eq "paired"){
				print "meta-velvetg $OUTDIR -ins_length $INS -amos_file yes -scaffolding yes -min_contig_lgth $min_contig_len >> $LOG\n";
			}
			elsif($usage eq "singles"){
				print "meta-velvetg $OUTDIR -amos_file yes -scaffolding yes -min_contig_lgth $min_contig_len >> $LOG\n";
			}
		}
		if ($transcripts){
			print "oases $OUTDIR  -amos_file yes -alignments yes >> $LOG\n";
		}
	}
	else{
		if (! $pair && ! $intlv && ($sCount==0)){die "[ERROR: $0] Check Singleton files\n"; }
		system("echo **************************** VELVETH ************************** >> $LOG");
		my $pid=&run("perl ".$scripts."usageStats.pl -i $interval -o usageStats_K$KMER.tsv");
		system("velveth $OUTDIR $KMER $useCase{$usage} >> $LOG");
		system("echo **************************** VELVETG ************************** >> $LOG");
		if ($usage eq "paired"){
			system("velvetg $OUTDIR -exp_cov auto -ins_length $INS -ins_length_sd $INS_SD -read_trkg yes -amos_file yes -min_contig_lgth $min_contig_len -unused_reads yes >> $LOG");
		}
		elsif($usage eq "singles"){
			system("velvetg $OUTDIR -exp_cov auto -read_trkg yes -amos_file yes -min_contig_lgth $min_contig_len -unused_reads yes >> $LOG");
		}
		$PIDs{$pid}++;
		if ($metaV){
			system("echo **************************** MetaVelvet ************************** >> $LOG");
			if ($usage eq "paired"){
				system("meta-velvetg $OUTDIR -ins_length $INS -ins_length_sd $INS_SD -amos_file yes -min_contig_lgth $min_contig_len -scaffolding yes >> $LOG");
			}
			elsif($usage eq "singles"){
				system("meta-velvetg $OUTDIR -amos_file yes -scaffolding yes -min_contig_lgth $min_contig_len >> $LOG");
			}
		}
		if ($transcripts){
			system("echo ***************************** OASES *************************** >> $LOG");
			system("oases $OUTDIR -amos_file yes -alignments yes >> $LOG");
		}
		if ($amos){
			warn "WARNING: The Meta-Velvet and AMOS portion of the Script is still under development. All the processes up to this point have been finished.\n";
			warn "This script will now exit\n";
			exit;
		}
	}
}


sub interleave{
	my $script=File::Spec->catfile( $scripts, "interleave.pl");
	print "Creating Output Directory: $OUTDIR\n";
	mkdir $OUTDIR || die "[ERROR $0] $!: $OUTDIR\n";
	my $out=File::Spec->catfile( $OUTDIR, $prefix );
	if (-e $script){
		print "Interleaving files:\t$fwd \; and\n\t$rev\n";
		if ($fasta){
			system("perl $script -fwd $fwd -rev $rev -prefix $out");
			print "perl $script -fwd $fwd -rev $rev -prefix $out\n";
		}
		else{
			system("perl $script -fwd $fwd -rev $rev -prefix $out -fastq");
			print "perl $script -fwd $fwd -rev $rev -prefix $out -fastq\n";
		}
		$intPair=$out."_int.".$seqType;
		push (@singles, $out."_sfwd.".$seqType, $out."_srev.".$seqType);
		$sCount++;
		return;
	}
	else{
		die "[ERROR: $0] ".$script.": Not Found!\n";
	}
}

sub getStats{
	my ($f, $o);
	if ($transcripts){
		$f="transcripts.fa";
		$o="trans";
	}
	elsif($metaV){
		$f="meta-velvetg.contigs.fa";
		$o="metav_contigs";
	}
	else{
		$f="contigs.fa";
		$o="contigs";
	}

	my $contigs=File::Spec->catfile( $OUTDIR, $f );
	my $out=File::Spec->catfile( $OUTDIR, "Log" );
	my $minLenFile=File::Spec->catfile( $OUTDIR, $o."_gt".$minLen.".fasta" );
	my $cwd=`pwd`;
	chomp $cwd;
	open(LOG, ">>".$out) || die "[ERROR: $0] Unable to write to Log file: $!\n";

	print LOG "\n# Pipeline version used:\tassemble.pl version:\t$version\n";
	print LOG "# All relevant files can be found in:\t$cwd\n";
	print LOG "# Assembly Directory:\t$OUTDIR\n";

	my $lim2len=File::Spec->catfile( $scripts, "limit2Length.pl");
	if (-e $lim2len){
		print LOG "\n# Helpful Stats for your sequences:\n";
		print "Calculating stats...\n";
		system("perl $lim2len -f $contigs -l $minLen -o $minLenFile >> $out");
		print LOG "\n";
	}

	my $calcN50=File::Spec->catfile( $scripts, "calcN50.pl");
	if (-e $calcN50){
		print LOG "# N50 stats before setting the $minLen limit:\n";
		system("perl $calcN50 $contigs >> $out");
		print LOG "\n";
		
		print LOG "# N50 stats after setting the $minLen limit:\n";
		system("perl $calcN50 $minLenFile >> $out");
		print LOG "\n";
	}
	
	my $searchNs=File::Spec->catfile( $scripts, "findStretchesOfNs.pl");
	my $searchNsOut=File::Spec->catfile( $OUTDIR, "stretchesOfN_k".$KMER.".out");
	if (-e $searchNs){
		print "# Finding long stretches of Ns...\n";
		print LOG `perl $searchNs -f $contigs -o $searchNsOut`;
	}
	system("mail -s 'Job: $OUTDIR Completed!' $email < $out");
}

sub identifyUser{
	my $uname=`whoami`;
	chomp $uname;
	my @groups=split(/\s+/,`groups`);
	my $found=0;
	foreach my$g(@groups){
		chomp $g;
		if ($g eq "gmb"){
			$email= $uname."\@umich.edu";
			$found++;
		}
	}
	unless($found){
		$email= "sunitj\@umich.edu";
	}

	return $email;
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

	print "in REAPER: ".$numPIDs."\n";
	while (scalar(keys %PIDs) > 0){
		my $pid= waitpid(-1, &WNOHANG);
		if ($pid > 0){
			print "in REAPER:$pid\n";
			if (WIFEXITED($?) && $PIDs{$pid}){
				`echo "Process ID: $pid\tFinished with status $?"`;
#				$numPIDs-- ;
				print "Process: ".$pid."\tStatus: ".$?."\nWaiting for ".scalar(keys %PIDs)." more processes...\n";
				delete $PIDs{$pid};
			}
		}
		else{
			sleep 10;
		}
	}
	return;
}

sub helpLoadingModules{
	print STDERR "## Required modules not found. Please load/install the following:\n";
	print STDERR "## REQUIRED: Velvet version 1.1.07-MAX99-OPENMP or higher\n";
	print STDERR "## For Metatranscriptomic Assembly: Oases version 0.2.01 or higher [OPTIONAL]\n";
	print STDERR "## For Metagenomic Assembly: MetaVelvet version 1.0.01 or higher [OPTIONAL]\n";
	print STDERR "## For combining multiple assemblies: AMOS version 3.1.0 or higher [OPTIONAL]\n";
	exit 1;
}
sub help{
	system('perldoc', $0);
	exit 1;
}
