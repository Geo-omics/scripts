#!/usr/local/bin/perl

=head1 DESCRIPTION

nsmpReport.pl -- Given AntiSmash and MetaPathways outputs get best hits and flanking regions of the genes of interest

=head1 USAGE

perl nsmpReport.pl -as path/to/AntiSmash/output/overview/directory -mp path/to/metapathways1/output/directory -out output.tsv

=head2 Options

    -antismash -as  <CHAR>  path to Antismash overview directory
    -metapathways -mp   <CHAR>  path to MetaPathways output directory.
    -prefix -p	<CHAR>	output file prefix   
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Tue Jan 20 15:25:16 EST 2015)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

# Core Perl Modules
use strict;
use Getopt::Long;
use File::Spec;
use FileHandle;
use File::Basename;
use File::Path;

# External Modules
use Bio::SeqIO;
use Bio::Graphics; # confirm requirements for this.

my $help;
my $version=fileparse($0)."\tv0.4.1";
my $nullValue="NA";
my ($asDir, $mpDir, $prefix, $overwrite);
my $imageDir;
GetOptions(
    'as|antismash:s'=>\$asDir,
    'mp|metapathways:s'=>\$mpDir,
    'p|prefix:s'=>\$prefix,
    'img|images:s'=>\$imageDir,
	'f|force'=>\$overwrite,
    'v|version'=>sub{print $version."\n"; exit;},
    'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

##############################
######   SANITY CHECK   ######
##############################

## Output file names
my($details, $concise, @outputs);
if (! $prefix){
	$details=$$."_detailed.tsv";
	$concise=$$."_concise.tsv";
	$prefix=$$;
}
else{
	$details=$prefix."_detailed.tsv";
	$concise=$prefix."_concise.tsv";
}

## Image Folder Name
if (! $imageDir){
	$imageDir=$prefix."_images";
}

## Check for pre-existing file/folders
my %preExist;
my $fatal=0;
unless ($overwrite){
	if(-e $details){
		$preExist{$details}++;
	}

	if (-e $concise){
		$preExist{$concise}++;
	}

	if (-d $imageDir){
		$preExist{$imageDir}++;
	}

	foreach (keys %preExist){
		print STDERR "[FATAL] '$_' already exists!\n";
		$fatal++;
	}
	if($fatal > 0){
		print STDERR "[ACTION] Use the '-force' or '-f' flag to overwrite already existing files and/or directories.\n";
		die "[MOCK] Try again...\n";
	}
}
undef(%preExist);

## Create Image Directory
if (! -d $imageDir) {
    mkpath($imageDir) || die $!;
}
elsif($overwrite){
	rmtree($imageDir) || die $!;
	mkpath($imageDir) || die $!;
}

## Set Color Pallete for image files
# Color Scheme for each layer
my @colors = qw(yellow cyan orange blue purple green chartreuse magenta aqua);
##############################
######		MAIN	    ######
##############################

## Scrape AntiSmash Data from EMBL files
my (%asClusters);
opendir(ASDIR, $asDir) or die $!;
while (my $subdir = readdir(ASDIR)) {
    next unless ($subdir=~ /^\d/);
    next unless (-d "$asDir/$subdir");
    my $seqDir=File::Spec->catdir($asDir, $subdir);
    my @embl=glob "$seqDir/*.embl";
	my ($panel, $clusters);
	my $pngFile=File::Spec->catfile($imageDir,"Clusters_".$subdir.".png");
    ($panel, $clusters)=parseEMBL($embl[0]);
	render($panel, $pngFile, $clusters, 'yellow');
}
closedir(ASDIR);

## Get a mapping key to talk between Original Names and MetaPathways Names 
my %nameMap;
my @mapFile=glob "$mpDir/preprocessed/*mapping.txt";
mapNames($mapFile[0]);

## Scrape MetaPathways annotation from GBK files
my (%annotation,%clusterGenes);
my @gbk=glob "$mpDir/genbank/*.gbk";
parseGBK($gbk[0]);

## Get Best hit information from MetaPathways Blast output
my %blast;
my @blastOut=glob "$mpDir/blast_results/*.refseq.blastout";
my $BLAST=FileHandle->new();
open($BLAST,"<",$blastOut[0]) || die $!;
while (my $line=<$BLAST>) {
    chomp $line;
    next if ($line=~ /^#/);
    next unless $line;
    
    my($locusID,$topHit, $percID, @stuff)=split(/\t/,$line);
    next if $blast{$locusID};
    
    my $evalue=$stuff[7];
    my $bitScore=$stuff[8];
    $blast{$locusID}{"Hit"}=$topHit;
    $blast{$locusID}{"PercID"}=$percID;
    $blast{$locusID}{"Evalue"}=$evalue;
    $blast{$locusID}{"Score"}=$bitScore;
}
close $BLAST;

## Bring it all together
my $header= join("\t", qw(Scaffold_Name Scaffold_Length Cluster_Number Cluster_Length Cluster_Start Cluster_Stop Cluster_Product Cluster_Rule Cluster_Monomer Cluster_Notes))."\t";
$header.=join ("\t", qw(Gene_Locus Gene_Start Gene_Stop Gene_Length Gene_Product))."\t";

$header.=join("\t", qw(Best_RefSeq_Hit Refseq_Hit_Percent_ID RefSeq_Evalue RefSeq_Bit_Score));

my $DETAIL=FileHandle->new();
open($DETAIL, ">", $details)|| die $!;
print $DETAIL $header."\n";

my $CONCISE=FileHandle->new();
open($CONCISE, ">", $concise) if ($concise);
print $CONCISE (join("\t", qw(Scaffold_Name Cluster_Number Cluster_Length Cluster_Product Cluster_Rule Cluster_Monomer Cluster_Notes Gene_Locus Gene_Length Gene_Product Best_RefSeq_Hit Refseq_Hit_Percent_ID)))."\n";
foreach my $scafName(keys %asClusters){
    foreach my $clustNum(keys %{$asClusters{$scafName}}){
        my ($clustStart, $clustStop)=sort{$a<=>$b}($asClusters{$scafName}{$clustNum}{"Start"}, $asClusters{$scafName}{$clustNum}{"Stop"});

        foreach my $locusID(keys %{$annotation{$scafName}}){
            my ($output,$cout)=("","");
	    my ($geneStart, $geneStop)=sort{$a<=>$b}($annotation{$scafName}{$locusID}{"START"}, $annotation{$scafName}{$locusID}{"STOP"});
            
            next if($clustStop < $geneStart); 
            next if($geneStop < $clustStart);

	    my $clustLen=$asClusters{$scafName}{$clustNum}{"ClustLen"};

# DETAILED REPORT: AntiSmash Data
	    $output=$scafName."\t".$asClusters{$scafName}{$clustNum}{"Scaffold_Length"}."\t".$clustNum."\t".$clustLen."\t";
	    $output.=$asClusters{$scafName}{$clustNum}{"Start"}."\t".$asClusters{$scafName}{$clustNum}{"Stop"}."\t";
	    $output.=$asClusters{$scafName}{$clustNum}{"Product"}."\t".$asClusters{$scafName}{$clustNum}{"Rule"}."\t";
		$output.=$asClusters{$scafName}{$clustNum}{"Monomer"}."\t".$asClusters{$scafName}{$clustNum}{"Notes"}."\t";

# DETAILED REPORT: MetaPathways Data
	    $output.=$locusID."\t".$annotation{$scafName}{$locusID}{"START"}."\t".$annotation{$scafName}{$locusID}{"STOP"}."\t".$annotation{$scafName}{$locusID}{"LEN"}."\t".$annotation{$scafName}{$locusID}{"PRODUCT"}."\t";

# DETAILED REPORT: RefSeq Blast Output
	    $output.=$blast{$locusID}{"Hit"}."\t".$blast{$locusID}{"PercID"}."\t".$blast{$locusID}{"Evalue"}."\t".$blast{$locusID}{"Score"};
	    print $DETAIL $output."\n";

# CONCISE REPORT: AntiSmash
	    $cout=$scafName."\t".$clustNum."\t".$clustLen."\t".$asClusters{$scafName}{$clustNum}{"Product"}."\t";
		$cout.=$asClusters{$scafName}{$clustNum}{"Rule"}."\t".$asClusters{$scafName}{$clustNum}{"Monomer"}."\t";
		$cout.=$asClusters{$scafName}{$clustNum}{"Notes"}."\t";
# CONCISE REPORT: MetaPathways
	    $cout.=$locusID."\t".$annotation{$scafName}{$locusID}{"LEN"}."\t".$annotation{$scafName}{$locusID}{"PRODUCT"}."\t";
# CONCISE REPORT: RefSeq Blast
	    $cout.=$blast{$locusID}{"Hit"}."\t".$blast{$locusID}{"PercID"};
	    
	    print $CONCISE $cout."\n";
        }
    }
}
close $DETAIL;
close $CONCISE;

exit 0;

##############################
######   SUB-ROUTINES   ######
##############################
sub render{
# Make the png image.
	my $panel=shift;
	my $pngFile=shift;
	my $clusters=shift;
	my $color=shift;
	my $adjusted=shift;
    my $idx    = 0;
	print "[RENDER]\n";
    for my $tag(sort keys %$clusters){
	my $featObj=${$clusters}{$tag};
	# my ($start, $stop, $strand)=\&gene_location();
	$panel->add_track($featObj,
				-offset		=> -500,
			    -glyph       => 'generic',
			    -bgcolor     => $color,
			    -fgcolor     => 'black',
			    -font2color  => 'red',
			    -key         => "Clusters",
			    -bump        => +1,
			    -height      => 12,
			    -label       => \&gene_product,
			    -description => \&gene_locus,
			    -stranded	=> +1,
			   );
    }

    my $PNG=FileHandle->new();
    open($PNG, ">", $pngFile)|| die $!;
    binmode($PNG);
    print $PNG $panel->png;
    close $PNG;
    return;
}

sub scaffold_info{
	my $seqObj=shift;

	# Definition of the sequence
	my $definition=$seqObj->desc;
	# LOCUS Identifier
	my $LocusID=$seqObj->display_id;
	# Contig/Scaffold Length
	my $Length=$seqObj->length;
	return ($definition, $LocusID, $Length);
}

sub gene_location{
	my $feature=shift;
	my ($start, $stop, $strand);
	# print "[LOCATION]\n";
	# Get notes, if present

	$start=$feature->location->start;
	$stop=$feature->location->end;
	$strand=$feature->location->strand;
	return ($start, $stop, $strand);
}

sub gene_locus{
	my $feature = shift;
	my (@notes,$locus);
	# Get Locus Tag
	foreach (qw(locus_tag)) {
	    @notes=eval{$feature->get_tag_values($_)};
        last;            
	}
	if (! $notes[0]){
		return gene_cluster($feature);
	}
	return $notes[0];
}

sub gene_product {
  my $feature = shift;
  my @notes;
  foreach (qw(product gene)) {
    @notes = eval {$feature->get_tag_values($_)};
  }
  return join(" ", @notes);
}
 
sub gene_cluster {
    my $feature = shift;
    my (@notes, $clustNum);
    foreach (qw(note)) {
		@notes = eval{$feature->get_tag_values($_)};
		foreach (@notes){
			chomp $_;
			if($_=~ /Cluster number: (\d+)/){
				$clustNum="$1";
				last;
			}
		}
    }
    if (length $clustNum > 10){
		substr($clustNum,10) = '...';
    }
    elsif(length $clustNum < 1){
		$clustNum="0'";
    }
    
    if(! @notes){
		$clustNum="UNK"
    }
    return "Cluster_".$clustNum;
}

sub generic_description {
  my $feature = shift;
  my $description;
  foreach ($feature->get_all_tags) {
	my @values = $feature->get_tag_values($_);
	$description .= $_ eq 'note' ? "@values" : "$_=@values; ";
  }
  $description =~ s/; $//; # get rid of last
  return $description;
}

sub makeBackbone{
	my $Length=shift;
	my $definition=shift;
	my $color=shift;
	my $start=shift;
	my $pad=shift;
	my $end=$Length;

	if(($pad) && ($start)){
		$end= ($start + $Length) + $pad;
		$start= ($start - $pad) < 1 ? 1 : ($start - $pad);
		$Length=($Length + (2*$pad));
	}
	$start=1 if (! $start);
	my $segment = Bio::Graphics::Feature->new(-start=>$start,-stop=>$end);
	my $panel;
	# Create a whole sequence object
	my $wholeseq = Bio::SeqFeature::Generic->new(
                                             -start        => $start,
                                             -end          => $end,
                                             -display_name => $definition
                                            );
	
	# Create the Image panel
	$panel = Bio::Graphics::Panel->new(
                                      -start    => $start,
									  -end		=> $end,
                                      -key_style => 'between',
                                      -width     => 1800,
                                      -pad_left  => 100,
                                      -pad_right => 100,
                                     );
	
	# Add the arrow track for scale
	$panel->add_track($wholeseq,
				  -start  => $start,
				  -end	  => $end,
                  -glyph  => 'arrow',
                  -bump   => 0,
                  -double => 1,
                  -tick   => 2,
                 );
	
	# Add a genereic blue line to represent the column
	$panel->add_track($wholeseq,
                  -glyph   => 'generic',
                  -bgcolor => $color,
                  -label   => 1,
                 );
	return($panel,$wholeseq);
}

sub parseEMBL{
    my $embl = Bio::SeqIO->new(
			    -format=> 'embl',
			    -file=> shift,
			    );
	my ($panel,%Clusters);
	print "[EMBL]\n";
    while (my $seqObj=$embl->next_seq) {
	# Contig Level Details that can be overwritten at the CDS level
        my($definition, $LocusID, $Length)=scaffold_info($seqObj);

		my $wholeseq;
		($panel, $wholeseq)=makeBackbone($Length,$definition,'blue');
        
	# Features for all cluster in file.
        foreach my $featObj ($seqObj->get_SeqFeatures){
			if($featObj->primary_tag eq "cluster"){
			# Add Cluster to panel
			push (@{$Clusters{"cluster"}}, $featObj);
			
		
			# Cluster Level details that should be specific to the current cluster. Can overwrite Contig level variables as well.
		    my($rule,$product, $start, $stop, $note, $cb, $clustNum,$monomer, $clustLen);
		    
			($start, $stop, undef)=gene_location($featObj);
		    # $start=$featObj->location->start;
		    # $stop=$featObj->location->end;
			$clustLen=abs($start - $stop)+1;
		
			# Get notes, if present
			if ($featObj->has_tag('note')) {
				my @notes=$featObj->get_tag_values("note");
				my @selNotes;
				foreach(@notes){
				if($_=~ /Cluster number: (.*)/){
					$clustNum=$1;
				}
				elsif ($_=~ /Detection rule/) {
					$rule=$_;
				}
				elsif ($_=~ /Monomers prediction: (.*)/) {
					$monomer=$1;
				}
				else{
					push(@selNotes,$_);
				}
				}
				$note=join(";", @selNotes);
			}
			else{
				$note=$nullValue;
			}
		
			# Get notes, if present
			if ($featObj->has_tag('clusterblast')) {
				my @clustBlast=$featObj->get_tag_values("clusterblast");
				$cb=join("\t", @clustBlast);
			}
			else{
				$cb=$nullValue;
			}
		
			$product=gene_product($featObj);
			
	#		print "$definition\t$clustNum\t$Length\n";
			$asClusters{$definition}{$clustNum}{"Start"}=$start;
			$asClusters{$definition}{$clustNum}{"Stop"}=$stop;
			$asClusters{$definition}{$clustNum}{"Product"}=$product;
			$asClusters{$definition}{$clustNum}{"Notes"}=$note;
			$asClusters{$definition}{$clustNum}{"ClustBlast"}=$cb;
			$asClusters{$definition}{$clustNum}{"Monomer"}=$monomer;
			$asClusters{$definition}{$clustNum}{"Rule"}=$rule;
			$asClusters{$definition}{$clustNum}{"ClustLen"}=$clustLen;
			$asClusters{$definition}{$clustNum}{"Scaffold_Length"}=$Length;
			my $clustDef=$definition."_Cluster_".$clustNum;
			my($clustPanel, $clustWholeSeq)=makeBackbone($clustLen,$clustDef, 'yellow');
			$asClusters{$definition}{$clustNum}{"Panel"}=$clustPanel;
			}
		}
    }
	return ($panel, \%Clusters);
}

sub mapNames{
    my $mapFile=shift;
#    print "Reading MP Mapping file\t".$mapFile."\n";
    my $FH=FileHandle->new();
    open($FH, "<", $mapFile) || die $!;
    while(my $line=<$FH>){
        chomp $line;
        next if ($line=~ /^#/);
        next unless $line;
        
        my($mpName, $origName)=split(/\t/, $line);
        $nameMap{$mpName}=$origName;
    }
    close $FH;
    return;
}

sub parseGBK{
	my $genbank = Bio::SeqIO->new(
                                -format=> 'genbank',
                                -file=> shift,
				-verbose=>-1
                                );
	print "[GBK]\n";

    while (my $seqObj=$genbank->next_seq) {
        # Contig Level Details that can be overwritten at the CDS level
        my($definition, $LocusID, $Length)=scaffold_info($seqObj);
		my $contig=$nameMap{$LocusID};

## NOTE ## $contig == $definition in hash $asCluster{$definition}	
        # Features for all CDSs in file.
		my (%clusterGenes, %adjusted);
        	foreach my $featObj ($seqObj->get_SeqFeatures){
			my $adjFeature;
        		if($featObj->primary_tag eq "CDS"){
				my $type="CDS";
                # CDS Level details that should be specific to the current CDS. Can overwrite Contig level variables as well.
                		my($locus,$product, $start, $stop, $strand,$clustNum, $adjStart, $adjStop);
				($start, $stop, $strand,$clustNum, $adjStart, $adjStop)=genes2cluster($contig,$featObj);
				next unless $clustNum;

				$locus=gene_locus($featObj);
				$product=gene_product($featObj);
# Add Cluster to panel
				push (@{$clusterGenes{$clustNum}}, $featObj);

				$annotation{$contig}{$locus}{"START"}=$start;
				$annotation{$contig}{$locus}{"STOP"}=$stop;
				$annotation{$contig}{$locus}{"TYPE"}=$type;
				$annotation{$contig}{$locus}{"LEN"}=($stop-$start)+1;
				$annotation{$contig}{$locus}{"STRAND"}=$strand;
				$annotation{$contig}{$locus}{"PRODUCT"}=$product;

				$adjusted{$featObj}{"START"}=$adjStart;
				$adjusted{$featObj}{"STOP"}=$adjStop;
			}
		}	
		
		foreach (keys %clusterGenes){
			my %Clusters;
			my $panel=$asClusters{$contig}{$_}{"Panel"};
			my $png=$contig."_Cluster_".$_.".png";
			my $pngFile=File::Spec->catfile($imageDir,$png);
			foreach my $feature(@{$clusterGenes{$_}}){
				push(@{$Clusters{$_}},$feature);
			}
			render($panel,$pngFile,\%Clusters, 'green', \%adjusted);
		}
	}
	return;
}

sub genes2cluster{
	my $contig=shift;
	my $featObj=shift;
	my ($clustNum, $adjStart, $adjStop);
	
	my ($start, $stop, $strand)=gene_location($featObj);

	my ($gStart, $gStop)=sort{$a <=> $b}($start,$stop);
	foreach (keys %{$asClusters{$contig}}){
		my $clustStart = $asClusters{$contig}{$_}{"Start"};
		my $clustStop = $asClusters{$contig}{$_}{"Stop"};
		my ($cStart, $cStop)=sort{$a<=>$b}($clustStart, $clustStop);

		next if($cStop < $gStart); 
	        next if($gStop < $cStart); 
		$clustNum=$_;
		$adjStart=$gStart-$cStart;
		$adjStop=$gStop-$cStart;
		#die "Tag 'adjusted_start' already exist in this Feature!\n" if($featObj->has_tag("adjusted_start"));
		#$featObj->add_tag_value("adjusted_start", $adjStart);		#$featObj->add_tag_value("adjusted_stop", $adjStop);
		last;
	}
	
	return($start, $stop, $strand,$clustNum, $adjStart, $adjStop);
}
