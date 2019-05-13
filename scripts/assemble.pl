#! /usr/bin/perl

# Copyright 2013, 2019 Regents of The University of Michigan.

# This file is part of geo-omics-scripts.

# Geo-omics-scripts is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.

# Geo-omics-scripts is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with Geo-omics-scripts.  If not, see <https://www.gnu.org/licenses/>.


=head1 NAME

assemble.pl - deprecated one stop shop for your illumina assembly needs


=head1 SYNOPSIS

Case 1: Assemble files as singletons

perl B<assemble.pl> -singles file1.fastq file2.fastq file3.fastq -k # (kmer value)

Case 2: Assemble files as paired; this will interleave the fwd and rev files and then assemble.

perl B<assemble.pl> -fwd fwd.fq -rev rev.fq -k # (kmer value) -i # (insert size)

Case 3: Assemble files as paired; when you already have an interleaved file

perl B<assemble.pl> -paired interleaved.fq -k # (kmer value) -i # (insert size)

Case 4: Mixed

perl B<assemble.pl> -paired interleaved.fq -singles file1.fastq file2.fastq file3.fastq -k # (kmer value) -i # (insert size)

OR

perl B<assemble.pl> -fwd fwd.fq -rev rev.fq -singles file1.fastq file2.fastq file3.fastq -k # (kmer value) -i # (insert size)


=head1 DESCRIPTION

Was created to be THE one stop shop for your illumina assembly needs. Once you're satified with the quality of your reads, put 'em here and the script will create the assemblies for you. The pipeline also includes support for Oases (for meta/transcriptomics) and MetaVelvet (for meta-genomics). Support for MetaIDBA and AMOS (minimus) is in the works.


=head1 OPTIONS

=head2 Required:

=over 8

=item B<-fwd>

Forward Sequence file

=item B<-rev>

Reverse Sequence file

=item B<-paired>

When you have an interleaved file, only accepts one file at a time, can be used with singles.

=item B<-singles>

When you have multiple files and you wish to force velvet into treating them individually.  You may give it as many files as long as the names are seperated by a space.

=item B<-k>, B<-kmer>

K-mer aka Hash size upto 99; We usually use k-mers of 61, 75, 91 but feel free to play around.

=item B<-i>, B<-insert>

Insert length; Required if using paired ended or interleaved data.  Look at your bioanalyzer results, let's say you see a peak at 391, that's your insert size.

=item B<-sd>

insert length standard deviation; default=13

=back

=head2 Optional

=over 8

=item B<-outdir>

The directory that contains the output

=item B<-p>, B<-prefix>

The prefix for your interleaved output

=item B<-interval>

The script also tracks the CPU and Memory consumption for the assemblies; default=10 (seconds)

=item B<-log>

change the name of the log file for the script

=item B<-contig_size>

Minimum Contig Length in your output; default=200

=back

=head2 Boolean Flags

=over 8

=item B<-v>

Version

=item B<-fasta>

If your sequence files are in the fasta format.

=item B<-debug>

When I need to debug the script. It doesnt actually run the assemblies, just prints out the commands it would have passed;

=back

=head2 Modifying Assembly Type

The following options require that the corresponding modules be loaded before executing the script.

=head3 Available

=over 8

=item B<-trans>

for meta/transcriptomic reads.(requires Oases)

=item B<-metav>

for metavelvet (metagenomic reads only).

=back

=head3 ToDo:

=head4 Assembler Support

=over 8

=item B<-idba>

uses IDBA-UD for assembly.

=back


=head1 DEPENDENCIES

=head2 Software

 Velvet 1.1.07 or above
 Oases 0.2.01 or above
 MetaVelvet or 1.0.01 or above

=head2 Scripts

 interleave - to interleave the forward and reverse fastq/fasta files
 calcN50 - to generate the N50 and L50 statistics of the final assembly
 findStretchesOfNs - to find the stretches of Ns longer than k-mer


=head1 Comments/Accolades/Brickbats/Beers:

 Sunit Jain, 2012
 sunitj AT umich DOT edu
 last updated: June 2013


=head1 SEE ALSO

L<omics(1)>, L<illumina-reads-processing(7)>

=head2 Other local resources

=over

=item [1]

L<HTML documentation|file:///usr/share/doc/geo-omics-scripts/html/index.html>

=item [2]

L<Omics workflow documentation [PDF]|file:///usr/share/doc/geo-omics-scripts/Geomicro-Illumina-Reads-Processing-Pipeline.pdf>

=back

=head2 Web

=over

=item [3]

L<Workflow documentation [PDF]|https://drive.google.com/open?id=0BxFSivK8RfJed05wamtrbEVUeE0>

=item [4]

L<Website|http://www.earth.lsa.umich.edu/geomicrobiology/>

=item [5]

L<Github repository|https://github.com/Geo-omics/scripts>

=back

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
		my $pid=&run("usageStats -i $interval -o usageStats_K$KMER.tsv");
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
	print "Creating Output Directory: $OUTDIR\n";
	mkdir $OUTDIR || die "[ERROR $0] $!: $OUTDIR\n";
	my $out=File::Spec->catfile( $OUTDIR, $prefix );
	print "Interleaving files:\t$fwd \; and\n\t$rev\n";
	if ($fasta){
		print "interleave  -fwd $fwd -rev $rev -prefix $out\n";
		system("interleave -fwd $fwd -rev $rev -prefix $out") == 0 or die "Failed running interleave script";
	}
	else{
		print "interleave  -fwd $fwd -rev $rev -prefix $out -fastq\n";
		system("interleave  -fwd $fwd -rev $rev -prefix $out -fastq") or die "Failed running interleave script";
	}
	$intPair=$out."_int.".$seqType;
	push (@singles, $out."_sfwd.".$seqType, $out."_srev.".$seqType);
	$sCount++;
	return;
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

	print LOG "\n# Helpful Stats for your sequences:\n";
	print "Calculating stats...\n";
	system("limit2Length -f $contigs -l $minLen -o $minLenFile >> $out") or die "Failed running limit2Length script";
	print LOG "\n";

	print LOG "# N50 stats before setting the $minLen limit:\n";
	system("calcN50 $contigs >> $out") == 0 or die "Failed running calcN50";
	print LOG "\n";
	
	print LOG "# N50 stats after setting the $minLen limit:\n";
	system("calcN50 $minLenFile >> $out") or die "Failed running calcN50";
	print LOG "\n";
	
	my $searchNsOut=File::Spec->catfile( $OUTDIR, "stretchesOfN_k".$KMER.".out");
	print "# Finding long stretches of Ns...\n";
	print LOG `findStretchesOfNs -f $contigs -o $searchNsOut`;
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
