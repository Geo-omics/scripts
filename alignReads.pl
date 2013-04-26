#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Cwd;
use File::Basename;
use Carp;

use Getopt::Long qw(:config no_ignore_case bundling);


$ENV{LC_ALL} = 'C';

my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  --left and --right <string>   (if paired reads)
#     or
#  --single <string>             (if unpaired reads)
#
#  Required inputs:
#   
#  --target <string>            multi-fasta file containing the target sequences (should be named {refName}.fa )
#
#  --seqType <string>           fa | fq    (fastA or fastQ format)
#
#  --aligner <string>           BLAT, bowtie, or bwa  (for spliced alignments, run BLAT or run tophat separately instead of this)
#
# Optional:
#
#  --SS_lib_type <string>       strand-specific library type:  single: F or R  paired: FR or RF
#                                   examples:  single RNA-Ligation method:  F
#                                              single dUTP method: R
#                                              paired dUTP method: RF
#
#  -o <string>                  output directory (default \${aligner}_out)
#
# 
#  --num_top_hits <int>         (default: 20)
#
#
#  If paired mode:
#
#     --max_dist_between_pairs             default (2000) 
#
#  ## BLAT-related pipeline options:
#
#  -I    maximum intron length  (default: 10000);  (only used in BLAT mode)
#  -P <int>                       min percent identity based on full sequence length  (default: 95)
#  --trim_short_terminal_segments     (trim off short terminal alignment segments that are mostly noise. Default: 10) 
# 
####################################################################################################################


_EOUSAGE_

	;


my $help_flag;
my $target_db; 
my $left_file;
my $right_file;
my $single_file;
my $max_intron = 10000;
my $output_directory = "";
my $trim_short_terminal_segment_length = 10;
my $min_per_ID = 95;
my $num_top_hits = 20;
my $max_dist_between_pairs = 2000;
my $seqType;
my $SS_lib_type;
my $aligner;

&GetOptions ( 'h' => \$help_flag,
			  
			  ## required inputs
			  'left=s' => \$left_file,
			  'right=s' => \$right_file,
			  
			  'single=s' => \$single_file,
			  

			  "target=s" => \$target_db,
			  "seqType=s" => \$seqType,
			  
              "aligner=s" => \$aligner,


			  ## Optional:
			  "SS_lib_type=s" => \$SS_lib_type,
			  
			  "I=i" => \$max_intron,
			  
			  'o=s' => \$output_directory,
			  
			  'trim_short_terminal_segments=i' => \$trim_short_terminal_segment_length,

			  'P=i' => \$min_per_ID,
			  
			  'num_top_hits=i' => \$num_top_hits,
			  			  
			  'max_dist_between_pairs=i' => \$max_dist_between_pairs,
			  
	);






if ($help_flag) { die $usage; }

unless ($target_db && -s $target_db) { 
    die $usage;
}

unless ($seqType && $seqType =~ /^(fq|fa)$/) {
    die $usage;
}

unless ($aligner && $aligner =~ /^BLAT|bowtie|bwa$/) {
    die $usage;
}

unless ( ($single_file && -s $single_file)
		 || 
		 ($left_file && -s $left_file
		  && $right_file && -s $right_file)) {
	die $usage;
}


unless ($output_directory) {
    $output_directory = "${aligner}_out";
}

if ($SS_lib_type && $SS_lib_type !~ /^(F|R|FR|RF)$/) {
	die "Error, SS_lib_type must be one of the following: (F, R, FR, RF)  ";
}



## check for required programs
{
    
    my @required_progs = qw(samtools);
    
    if ($aligner eq "BLAT") {
        push (@required_progs, qw (blat psl2sam.pl samtools) );
    }
    else {
        push (@required_progs, qw (bowtie-build bowtie) );
    }
    
	foreach my $prog (@required_progs) {
		my $path = `which $prog`;
		unless ($path =~ /^\//) {
			die "Error, path to required $prog cannot be found";
		}
	}
}

main: {


	my $start_dir = cwd();

	$single_file = &build_full_paths($single_file, $start_dir) if $single_file;
	$left_file = &build_full_paths($left_file, $start_dir) if $left_file;
	$right_file = &build_full_paths($right_file, $start_dir) if $right_file;
	$target_db = &build_full_paths($target_db, $start_dir);
	
	my $work_dir;
	if ($output_directory =~ /^\//) {
		$work_dir = $output_directory;
	}
	else {
		$work_dir = "$start_dir/$output_directory";
	}
	
	&process_cmd("mkdir -p $work_dir") unless (-d $work_dir);
	
	my $util_dir = "$FindBin::Bin/../util";
	
	my @entries;
	
	if ($single_file) {
		push (@entries, ["single_dir", "single_fa", $single_file]);
	}
	else {
		push (@entries, 
			  ["left_dir", "left_fa", $left_file],
			  ["right_dir", "right_fa", $right_file]);
	}
    
    chdir $work_dir or die "Error, cannot cd to $work_dir";
	
	unless (-s "target.fa") {
		
		# prep the target here for converting sam to bam later.
		
	
        my $cmd = "ln -s $target_db target.fa";
		&process_cmd($cmd);
	
	
        unless (-s "target.fa.fai") {
            my $cmd = "samtools faidx target.fa";
            &process_cmd($cmd);
        }
        
        ## prep reference if not BLAT
        if ($aligner eq "bowtie") {
            my $cmd = "bowtie-build target.fa target";
            &process_cmd($cmd);
        }
        elsif ($aligner eq "bwa") {
            my $cmd = "bwa index target.fa";
            &process_cmd($cmd);
        }
    }
    
    my $num_hits = ($single_file) ? $num_top_hits : 2 * $num_top_hits;
    
    my @process_monitor_files;
    
    
	foreach my $fa_info (@entries) {
		
		## always resume work in the work_dir
		chdir $work_dir or die "Error, cannot cd to $work_dir";
		        
		my ($target_dir, $target_fa, $trans_file) = @$fa_info;
        
        my $process_monitor_file = &build_full_paths("$target_dir/process_info.txt", $work_dir);
        push (@process_monitor_files, $process_monitor_file);
                
        
        my $pid = fork();
        if ($pid) {
            # parent process
            next;
        }
        
        # child does the work:
        eval {
            unless (-d $target_dir) {
                mkdir $target_dir or die "Error, cannot mkdir $target_dir";
            }
            
            # prep target dir with symlinks to components needed.
            if ($seqType eq "fq" && $aligner eq "BLAT") {
                my $cmd = "$util_dir/fastQ_to_fastA.pl -I $trans_file > $target_dir/$target_fa";
                &process_cmd($cmd) unless (-e "$target_dir/$target_fa");
            }
            else {
                ## already fasta format or running in bowtie mode, which can use fq format
                my $cmd = "ln -s $trans_file $target_dir/$target_fa";
                &process_cmd($cmd) unless (-e "$target_dir/$target_fa");
            }
            
            my $cmd = "ln -s ../target.fa* .";
            &process_cmd($cmd) unless (-e "target.fa");
            
            
            ## Work in Target_dir
            chdir ($target_dir) or die "Error, cannot cd to $target_dir";
            
            
            if ($aligner eq "BLAT") {
                ## run BLAT alignment pipeline
                
                $cmd = "ln -s ../target.fa .";
                &process_cmd($cmd) unless (-e "target.fa");
                
                ## Prep sequences>
                if ($seqType eq "fq") {
                    $cmd = "$util_dir/fastQ_to_tab.pl -I $trans_file > $target_fa.tab";
                    &process_cmd($cmd) unless (-s "$target_fa.tab");
                }
                else {
                    $cmd = "$util_dir/fasta_to_tab.pl < $trans_file > $target_fa.tab";
                    &process_cmd($cmd) unless (-s "$target_fa.tab");
                }
                
                $cmd = "sort -T . -S 2G -k 1,1 -k 3,3 $target_fa.tab > $target_fa.sorted.tab";
                &process_cmd($cmd) unless (-s "$target_fa.sorted.tab");
                
                ## run blat
                $cmd = "$util_dir/run_BLAT_shortReads.pl target.fa $target_fa $max_intron $target_fa.psl";
                &process_cmd($cmd) unless (-s "$target_fa.psl");
                
                ## convert to sam
                $cmd = "psl2sam.pl -q 0 -r 0 $target_fa.psl > $target_fa.psl.sam";
                &process_cmd($cmd) unless (-s "$target_fa.psl.sam");
                
                ## sort by name
                $cmd = "sort -T . -S 2G -k 1,1 -k 3,3 $target_fa.psl.sam > $target_fa.psl.nameSorted.sam";
                &process_cmd($cmd) unless (-s "$target_fa.psl.nameSorted.sam");
                
                ## add sequences to sam file.
                $cmd = "$util_dir/blat_sam_add_reads2.pl $target_fa.psl.nameSorted.sam $target_fa.sorted.tab > $target_fa.psl.nameSorted.wReads.sam";
                &process_cmd($cmd) unless (-s "$target_fa.psl.nameSorted.wReads.sam");
                
                ## capture top hits
                $cmd = "$util_dir/top_blat_sam_extractor.pl $target_fa.psl.nameSorted.wReads.sam $num_hits $min_per_ID > $target_fa.nameSorted.sam"; # capture twice as many hits in single mode, distill based on pairs later.
                &process_cmd($cmd) unless (-s "$target_fa.nameSorted.sam");
                
                
            }
            elsif ($aligner eq "bowtie")  {
                ## Run Bowtie Alignment pipeline
                
                my $cmd = "ln -s ../target.* .";
                &process_cmd($cmd) unless (-e "target.fa");
                
                
                my $format = ($seqType eq "fq") ? "-q" : "-f";
                
                $cmd = "bowtie -q -S --sam-nohead $format -v 2 -k $num_hits target $target_fa > $target_fa.pre.sam";
                &process_cmd($cmd) unless (-e "$target_fa.pre.sam");
                
                ## remove unaligned reads
                $cmd = "$util_dir/SAM_filter_out_unmapped_reads.pl $target_fa.pre.sam > $target_fa.sam";
                &process_cmd($cmd) unless (-e "$target_fa.sam");
                
                # name-sort 
                $cmd = "sort -T . -S 2G -k 1,1 -k 3,3 $target_fa.sam > $target_fa.nameSorted.sam";
                &process_cmd($cmd) unless (-e "$target_fa.nameSorted.sam");
            }
            elsif ($aligner eq "bwa") {
                
                my $cmd = "ln -s ../target.* .";
                &process_cmd($cmd) unless (-e "target.fa");
                
                $cmd = "bwa aln target.fa $target_fa > $target_fa.sai";
                &process_cmd($cmd) unless (-e "$target_fa.sai");
                
                $cmd = "bwa samse -n $num_hits target.fa $target_fa.sai $target_fa > $target_fa.pre.sam";
                &process_cmd($cmd) unless (-e "$target_fa.pre.sam");
                
                ## remove unaligned reads
                $cmd = "$util_dir/SAM_filter_out_unmapped_reads.pl $target_fa.pre.sam > $target_fa.sam";
                &process_cmd($cmd) unless (-e "$target_fa.sam");
                
                # name-sort 
                $cmd = "sort -T . -S 2G -k 1,1 -k 3,3 $target_fa.sam > $target_fa.nameSorted.sam";
                &process_cmd($cmd) unless (-e "$target_fa.nameSorted.sam");
            }
            
            
        
        };

        # write monitor result file
        open (my $ofh, ">$process_monitor_file") or die "Error, cannot write to $process_monitor_file";
        
        if ($@) {
            print $ofh "ERROR\t$@\n";
        }
        else {
            print $ofh "SUCCESS\n";
        }
        close $ofh;
        
        exit(0); # child exits. Parent never reaches this statement.
    }
    
    # wait for children to stop running alignments.
    while (wait() > 0) {
        print "-child alignment process completed.\n";
    }
    
    
	chdir $work_dir or die "Error, cannot cd to $work_dir";
	
    ## ensure that each alignment process completed successfully.
    
    my $failure = 0;
    foreach my $monitor_file (@process_monitor_files) {
        if (! -s $monitor_file) {
            print STDERR "Error, cannot find process monitor file: $monitor_file";
            $failure = 1;
            next;
        }
        open (my $fh, $monitor_file);
        my $status_info = <$fh>;
        chomp $status_info;
        my ($status, $msg) = split(/\t/, $status_info);
        if ($status eq "SUCCESS") {
            # no op
        }
        else {
            print STDERR $status_info;
            $failure = 1;
        }
    }
    
    if ($failure) {
        die "Error, alignment pipeline failed due to errors encountered above.";
    }
    else {
        print "\n## Alignment steps succeeded.\n\n";
    }
    
    
    ## NONE OF THE COMMANDS BELOW ARE SKIPPED, INSTEAD OUTPUTS ARE OVER-WRITTEN
    

	## merge into single sam file, setting flags properly

	my $outfile_basename = basename($output_directory);
	
	if ($single_file) {
		
		my $cmd = "sort -T . -S 2G -k 3,3 -k 4,4n single_dir/single_fa.nameSorted.sam > single_dir/single.coordSorted.sam";
		&process_cmd($cmd);
        
        $cmd = "ln -sf single_dir/single.coordSorted.sam $outfile_basename.pre.coordSorted.sam";
        &process_cmd($cmd);
                
	}

	else {
		## paired mode:
	
        ## Now, capture just the top number of hits taking into account read pairing info.
		
		my $cmd = "$util_dir/merge_left_right_nameSorted_SAMs.pl --left_sam left_dir/left_fa.nameSorted.sam --right_sam right_dir/right_fa.nameSorted.sam -C $num_top_hits -D $max_dist_between_pairs > combined.nameSorted.sam";
		&process_cmd($cmd);
		
		## sort by coordinate.
		
		$cmd = "sort -T . -S 2G -k 3,3 -k 4,4n combined.nameSorted.sam > $outfile_basename.pre.coordSorted.sam";
		&process_cmd($cmd);
	}

	if ($aligner eq "BLAT") {
        

        # report splice junctions and remove short terminal exons that are more likely noise.
        my $cmd = "$FindBin::Bin/../Inchworm/bin/cigar_tweaker $outfile_basename.pre.coordSorted.sam target.fa $trim_short_terminal_segment_length | sort -T . -S 2G -k 3,3 -k 4,4n >  $outfile_basename.coordSorted.spliceAdjust.sam";
        &process_cmd($cmd);
    }
    else {
        my $cmd = "ln -sf $outfile_basename.pre.coordSorted.sam $outfile_basename.coordSorted.spliceAdjust.sam";
        &process_cmd($cmd);
    }
	
	# add transcribed orientation info:
	if ($SS_lib_type) {
		my $cmd = "$util_dir/SAM_set_transcribed_orient_info.pl $outfile_basename.coordSorted.spliceAdjust.sam $SS_lib_type > $outfile_basename.coordSorted.sam";
		&process_cmd($cmd);
	}
	else {
		# not strand-specific, keep as is and don't disrupt current flow (so use expected output name)
		my $cmd = "ln -sf  $work_dir/$outfile_basename.coordSorted.spliceAdjust.sam $outfile_basename.coordSorted.sam";
		&process_cmd($cmd);
	}
	
	# convert to bam format
	my $cmd = "samtools view -bt target.fa.fai -S $outfile_basename.coordSorted.sam > $outfile_basename.coordSorted.bam";
	&process_cmd($cmd);
	
		
    ## provide name-sorted SAM
    $cmd = "sort -T . -S 2G -k 1,1 -k 3,3 $outfile_basename.coordSorted.sam > $outfile_basename.nameSorted.sam";
    &process_cmd($cmd);
	
    $cmd = "samtools view -bt target.fa.fai $outfile_basename.nameSorted.sam > $outfile_basename.nameSorted.bam";
    &process_cmd($cmd);


	$cmd = "samtools index $outfile_basename.coordSorted.bam";
	&process_cmd($cmd);
    

	if ($SS_lib_type) {
		## strand-specific
		## separate the sam based on strand, and create separate bam files.  (for convenience sake)
		
		$cmd = "$util_dir/SAM_strand_separator.pl $outfile_basename.coordSorted.sam $SS_lib_type";
		&process_cmd($cmd);

        $cmd = "$util_dir/SAM_strand_separator.pl $outfile_basename.nameSorted.sam $SS_lib_type";
        &process_cmd($cmd);
        
		foreach my $sam_file ("$outfile_basename.coordSorted.sam.+.sam", "$outfile_basename.coordSorted.sam.-.sam",
                              "$outfile_basename.nameSorted.sam.+.sam", "$outfile_basename.nameSorted.sam.-.sam") {
			
			if (-s $sam_file) {

				my $bam_file = $sam_file;
				$bam_file =~ s/\.sam$/\.bam/;
				
				$cmd = "samtools view -bt target.fa.fai $sam_file > $bam_file";
				&process_cmd($cmd);
				
				$cmd = "samtools index $bam_file";
				&process_cmd($cmd) if ($bam_file =~ /coordSorted/);
			}
		}
	}
	
	
	exit(0);
}


####
sub process_cmd {
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";

	my $ret = system($cmd);

	if ($ret) {
	    confess "Error, cmd: $cmd died with ret $ret";
	}

	return;
}


####
sub build_full_paths {
	my ($path, $start_dir) = @_;
	
	if ($path && $path !~ /^\//) {
		$path = "$start_dir/$path";
	}

	return($path);
}
