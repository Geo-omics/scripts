#! /usr/bin/perl

=head1 ESOM Wrapper version 0.1.1

=head1 DESCRIPTION

	Give me:
		- Path to folders that has the fasta files.

	I will:
		- create an annotation file and a concatenated fasta for ESOM binning;
		- run the tetramer frequency script on the files;

=head1 USAGE

	perl esomWrapper.pl -path Folder_Path -ext extension_of_files

=head2 Options

	-p or path	:	path to folder containing fasta files; use "." (dot, without the quotes) for current folder.
	-e or ext	:	file extension to look for in folder; default= fasta
	-prefix	:		prefix filename for annotation and concatenated file; default=esom
	-DIR or dir	:	name of the output directory; default= ESOM
	-min	:	Optional	default=2500; Minimal length (in nt) of input contig to be included in output
	-max	:	Optional	default=5000
	Note:	The script will split sequence after each 'max' nt; join last part, if remaining seq shorter than 'max', with second-last part
			eg: in default settings, a sequence of 14 kb will be split into a 5 kb and a 9 kb fragment if window_size = 5 kb.
	-h	this page.

=head3 Example 1: Required Options

	perl esomWrapper.pl -path .

=head3 Example 2: Other Options

	perl esomWrapper.pl -path . -ext fa -dir MyESOM -prefix esomOutput -min 2000 -max 5000

=head1 Suggestions/Corrections/Feedback/Beer

	Sunit Jain, sunitj@umich.edu
	January 2013

=cut

use strict;
use Getopt::Long;
use File::Spec;
use File::Basename;

### Wrapped Scripts ###
my $tetramerScript="/geomicro/data1/COMMON/scripts/tetramer_freqs_esom.pl";
my $codonModScript="/geomicro/data1/COMMON/scripts/esomCodonMod.pl";
###

my $version="0.1.2";
my $path; # Folder path
my $ext="fasta";
my $prefix="esom";
my $outDir="ESOM";
my $kmer = 4;
my $noMod;
my $min_length = 2500; #Minimal length (in nt) of input contig to be included in output
my $window_size = 5000; #split sequence after each window_size nt, 
                     #join last part, if shorter than window_size, 
                     #with second-last part (a sequence of 
                     #14 kb will be split into a 5 kb and a
                     #9 kb fragment if window_size = 5 kb)

GetOptions(
	'tetramer:s'=>\$tetramerScript,
	'p|path:s'=>\$path,
	'e|ext:s'=>\$ext,
	'k|kmer:i'=>\$kmer,
	'prefix:s'=>\$prefix,
	'DIR|dir:s'=>\$outDir,
	'min:i'=>\$min_length,
	'max:i'=>\$window_size,
	'no_mod'=>\$noMod,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

die "[ERROR: $0] Folder Path Required! See $0 -h for help on the usage" if !$path;

my $annotationFile=$prefix.".ann";
my $concatenatedFasta=$prefix.".".$ext;
my $logFile=$prefix.".log";
if (-e $outDir){
	die "[ERROR: $0]$outDir already exists!\n";
}
else{
	mkdir($outDir, 0755);
}

my $ann=File::Spec->catfile( $outDir, $annotationFile);
my $catFasta=File::Spec->catfile( $outDir, $concatenatedFasta);
my $log=File::Spec->catfile( $outDir, $logFile);

my @files=<$path/*.$ext>;

#$|++;

open(FASTA, ">".$catFasta) || die $!;
open(ANN, ">".$ann) || die $!;
open(LOG, ">".$log) || die $!;

my $class=0;
my $filesProcessed=0;
print "## ESOM Wrapper version: $version ##\n";
print "# FileName\tNumber of Sequences found\n";
print LOG "# FileName\tClass Assigned\tNumber of Sequences\n";
print ANN "# Contig\tAnnotation\tClass\n";
my %cls;
foreach my $file(@files){
	
	my $countSeqs=	parseFasta($file);
	my $fileName = basename($file, ".".$ext);
	$fileName=~ s/\s+/\_/g;
	print $fileName."\t".$countSeqs."\n";
	print LOG $fileName."\t".$class."\t".$countSeqs."\n";
	
	$cls{$class}=$fileName;
	$class++;
	$filesProcessed++;
}
close(IN);
close(FASTA);
print "\n# Files processed:\t $filesProcessed\n";
print LOG "\n# Files processed:\t $filesProcessed\n";
close (LOG);

print "# Calculating Tetramer Frequencies...\n";
chdir $outDir || die $!;
system("perl $tetramerScript -f $concatenatedFasta -a $annotationFile -min $min_length -max $window_size -ext $ext -kmer $kmer >> $prefix.log");

unless($noMod){
	print "# Applying Codon Modification...\n";
	my $lrnfile ="Tetra_".$prefix."_".$min_length.".lrn";
	my $outLrnfile = "Tetra_".$prefix."_".$min_length.".mod.lrn";
#	print "perl $codonModScript -lrn $lrnfile -o $outLrnfile >> $prefix.log\n";
	system("perl $codonModScript -lrn $lrnfile -o $outLrnfile >> $prefix.log");
}

print "# Adding class names and colors to the cls file\n";
my $clsFile= "Tetra_".$prefix."_".$min_length.".cls";
my $tmpCls="tmp.cls";
open(CLS, $clsFile)|| die $!;
open(TMP, ">".$tmpCls)|| die $!;
while(my $line=<CLS>){
	if ($.==2){
#		print TMP"\%0\t".$cls{0}."\t255\t255\t255\n";
		for(my $i=0; $i<  $filesProcessed; $i++){
			my $clsColor=randomColors();
			print TMP "\%".$i."\t".$cls{$i}."\t".$clsColor."\n";
		}
	}
	print TMP $line;
}

system("mv $tmpCls $clsFile");

print STDERR "\nAll done! please check the $prefix.log file for errors/warnings before you proceed\n";
print STDERR "Also make sure that your class (.cls) files have values in both the columns\n";
system ("echo \"All done! please check the $prefix.log file for errors/warnings before you proceed\" >> $prefix.log");
system ("echo \"Also make sure that your class (.cls) files have values in both the columns\" >> $prefix.log");
exit 0;

#### Sub-routines ####
sub parseFasta{
	my $fileName=shift;

	open(IN, $fileName) || die $fileName.":".$!."\n";

	my ($prevHeader, $flag);
	$/=">";
	my $countSeqs=0;
	while(my $line=<IN>){
		chomp $line;
		$line=~ s/\r//;
		next unless $line;

		my($header, @sequence)=split(/\n/, $line);
		my $seq= join("", @sequence);
		$seq=~ s/\s+//g;
		if (length($seq)==0){
			$prevHeader=$header;
			$flag=1;
			next;
		}
		elsif($flag==1){
			$header=$prevHeader."_".$header;
			$flag==0;
			$prevHeader="";
		}

		# Beautify
		$header=~ s/\s+/\_/g;
		$header=~ s/\W+/\_/g;
		$header=~ s/\_+/\_/g;
		$header=~ s/\_+$//;
		$header=~ s/^\_+//;

		$countSeqs++;
		print FASTA ">".$header."\n".$seq."\n";
		print ANN $header."\t".$header."\t".$class."\n";
	}
	$/="\n";
	return $countSeqs;
}

sub randomColors {
    my ($r, $g, $b) = map { int rand 256 } 1 .. 3;
    my $color= join("\t", $r, $g, $b);
    return ($color);
}

