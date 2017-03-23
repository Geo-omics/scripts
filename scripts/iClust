#! usr/bin/perl
#USAGE: perl iClust.pl <inputFile> <id1> <id2>
#	INPUT FILE: should contain a list of file names of individual genomes. each name should be in a newline.


use strict;
use Benchmark;

sub uClust{
	my $genome1=shift;
	my $iden1=shift;
	my $id1=$iden1/100;
	my ($fName1, $ext1)=split (/\./, $genome1);
	my $sorted1=$fName1.".sorted";
	my $out1=$fName1."_".$iden1.".count";
	system("uclust --mergesort ".$genome1." --output ".$sorted1." --amino");
	system("uclust --input ".$sorted1." --uc ".$out1." --id ".$id1." --amino");
	my (%c1)=clustCount($out1);
	return (%c1);

}

sub clustCount{
	my $genome2=shift;
	open (CLUST, "$genome2") || die "ERROR: $genome2\n$!";
	my %clusters2;
	while (my $line2=<CLUST>){
		unless ($line2=~/^\#|^C/){
			chomp $line2;
			my @stuff2=split(/\t/, $line2);
			my ($number2, @etc2)=split (/\s/, $stuff2[-1]);
			if ($line2=~/^S/){
				if ($clusters2{$stuff2[1]}){
					push (@{$clusters2{$stuff2[1]}}, $number2);
				}
				elsif (! $clusters2{$stuff2[1]}){
					push (@{$clusters2{$stuff2[1]}}, $number2);
				}
			}	
			
			elsif($line2=~/^H/){
				if ($clusters2{$stuff2[1]}){
					push (@{$clusters2{$stuff2[1]}}, $number2);
				}
				elsif (! $clusters2{$stuff2[1]}){
					push (@{$clusters2{$stuff2[1]}}, $number2);
				}	
			}	
		}
	}
	
	my($fName2, $ext2)=split (/\./, $genome2);
	my $out2=$fName2.".clust";
	open (OUT, ">$out2");
	print OUT "\#CNum\tTotal\tClusteredSeqs\n";

	my $totalGenes2=0;
	foreach my $key2(keys %clusters2){
		my $arraySize2= scalar (@{$clusters2{$key2}});
		my $everything2= join("\t", @{$clusters2{$key2}});
		print OUT $key2."\t".$arraySize2."\t".$everything2."\n";
		$totalGenes2=$totalGenes2+$arraySize2;
	}
	return (%clusters2);	
}
sub prepFasta{
	my $in3=shift;
	my $numGenomes3=shift;
	my %wg3;
	chomp($in3);
	open (FAA, $in3) || die "ERROR: $in3\n".$!;
	$/= ">";
	while (my $a3 = <FAA>) {
		chomp $a3;
		next unless $a3;
		my ($query3, @seqs3) = split (/\n/, $a3);
		my($numbers3, @etc3)=split(/\s/, $query3);
		chomp($numbers3);
		my $nuQuery3=$numbers3."@".$numGenomes3;
		my $seq3 = join ("", @seqs3);
		chomp($nuQuery3, $seq3);
		$wg3{$nuQuery3}=$seq3;
		#print "$nuQuery3\t";
	}
	close FAA;
	$/= "\n";
	
	## print out edited version.
	my($fName3, @ext3)=split(/\./, $in3);
	my $out3=$fName3."_ed.faa";
	open (ED, ">".$out3);
	while (my($k3, $v3)=each(%wg3)){
		chomp($k3, $v3);
		print ED ">".$k3."\n";
		print ED $v3."\n";
	}
	close ED;
	undef %wg3;
	return $out3;
}

sub makeItReadable{
	my $in4=shift;
	my $numGenomes4=shift;
	my %finalClusters4=%{$in4};
	my %output4;
	while (my($k4,$v4)=each(%finalClusters4)){
		my %tracker4;
		@{$output4{$k4}}=(0) x $numGenomes4;
		foreach my $geneName4(@{$v4}){
			my($gene4,$pos4)=split(/\@/, $geneName4);
			$tracker4{$pos4}++;
		}
		while(my($key4,$value4)=each(%tracker4)){
			$output4{$k4}[$key4-1]=$value4;
		}
	}	
	return (%output4);
}

sub magic{
	my $allCombined5=shift;
	my $id5=shift;
	my $go2nxt5=shift;
	my $numGenomes5=shift;
	my $editedList5=shift;
	my @edGenomeList5=@{$editedList5};
	my $allSeq5=shift;
	my %allSequences5=%{$allSeq5};
	my $level5=shift;
	print "inFile:$allCombined5\n";
	my (%clusters5)=uClust($allCombined5, $id5);
	print "\nTotal Clusters:\t".keys(%clusters5)."\n";

	# deleting Temporary files.
	system("rm -f *.sorted");
	system("rm -f *.count");
	system("rm -f *_ed.faa");

	print "\nRecieved Gibberish...\nMaking sense out of the data...\n";
	my %readable5=makeItReadable(\%clusters5, $numGenomes5);
	my ($rfName5, @etc5)=split(/\./, $allCombined5);
	my $result5="result_".$rfName5."G".$numGenomes5."id".$id5.".tsv";
	open(RES, ">$result5");

	print "$numGenomes5 were read\n";
	print "Total Clusters Read:".keys(%readable5)."\n";

	print RES "\#CNum\t";
	foreach my $gName5(@edGenomeList5){
		print RES $gName5."\t"
	}
	print RES "TotalGenes\tLumps\n";
	print "Done!!\n\n";

	my $summary5="summary_".$rfName5."_".$id5.".log";
	open (RSUMM, ">$summary5");
	print "These Clusters were present in all Genomes at $id5\% identity:\n";
	print RSUMM "\#These Clusters were present in all Genomes:\n";
	print RSUMM "\#Refer to the allCombined_$ARGV[1].clust file for more details.\n";
	print RSUMM "\#ClustNum\tGenesInCluster\tGenes\n";

	my %shortlist5;
	my $totalGenes5=my $level5_1=my $allPresent5=my $flag5=0;
	while (my($k5,$v5)=each(%readable5)){
		print RES $k5."\t";
		my ($total5, $presence5);
		foreach my $value5(@{$v5}){
			print RES $value5."\t";
			$total5+= int($value5);
			$presence5++ if ($value5!=0);
			if ($presence5 == $numGenomes5){
				print $k5."\t".$total5."\n";
				my $details5=join ("\t", @{$clusters5{$k5}});
				print RSUMM $k5."\t".$total5."\t".$details5."\n";
				$allPresent5++;
				$totalGenes5+=$total5;
			}
			if (($total5 >= int($numGenomes5)) && $level5 ==0){
				$flag5=1;				
			}
		}
		if ($flag5==1){
			$shortlist5{$k5}=$clusters5{$k5};
			$flag5=0;
			$go2nxt5=1;
			$level5_1++;
		}
		print RES "$total5\t$presence5\n";
	}
	

	close RES;
	print "TotalClusters:\t$allPresent5\n";
	print "TotalGenes in these Clusters:\t$totalGenes5\n";
	print RSUMM "Total:\t$allPresent5\n";
	print RSUMM "TotalGenes in these Clusters:\t$totalGenes5\n";
	close RSUMM;
	print "Current Level:\t$level5\n";
	print "$level5_1 clusters qualified for the next level.\n";
	if ($go2nxt5 == 1){
		if (scalar(keys %shortlist5) >0){
			print "Short5:";
			print keys(%shortlist5)."\t";
			print ".";
			my @iterate5=nextLevel(\%shortlist5, \%allSequences5);
			andAgain(\@iterate5, $numGenomes5, \@edGenomeList5, \%allSequences5);
		}
	}
	return();
}


sub andAgain{
	print "Going at it again\n";
	my $lof6=shift;
	my @listOfFiles6=@{$lof6};
	my $id6=$ARGV[2];
	my $go2nxt6=0;
	my $numGenomes6=shift;
	my $editedList6=shift;
	my @edGenomeList6=@{$editedList6};
	my $allSeq6=shift;
	my %allSequences6=%{$allSeq6};
	my $level6=1;
	foreach my $file6(@listOfFiles6){
		&magic($file6,$id6, $go2nxt6, $numGenomes6, \@edGenomeList6, \%allSequences6, $level6);
	}
	print "Iterative UClust Successfully Completed!!\n";
	return();
}

sub nextLevel{
	print "New Level\n";
	my $sl7=shift;
	my %nLevel7=%{$sl7};
	my $sq7=shift;
	my %allCombined7=%{$sq7};
	my @clusterFiles7;
	while (my($cNum7, $seqNames7)=each(%nLevel7)){
		my $clusterFAA7="Cluster".$cNum7.".faa";
		open (CLEVEL, ">".$clusterFAA7);
		foreach my $sn7(@{$seqNames7}){
			print CLEVEL">$sn7\n$allCombined7{$sn7}\n";	
		}
		close CLEVEL;
		push (@clusterFiles7, $clusterFAA7);
	}
	print "All Cluster Files Created\n";
	return(@clusterFiles7);
}

sub allSeqs{
	my $in8=shift;	

	open (CONTIGS, $in8) || die "Couldn't open $in8\n";
	$/= ">";
	my %sequences8;
	while (my $b8 = <CONTIGS>) {
	    chomp $b8;
	    next unless $b8;
	    my ($name8, @sequence8) = split (/\n/, $b8);
	    my $seq8 = join ("", @sequence8);
	    $sequences8{$name8} = uc $seq8;
	}
	close CONTIGS;
	$/="\n";	
	return (%sequences8);
}

### MAIN ###

#start timer
my $start=new Benchmark;

chomp ($ARGV[0]); # files list
chomp ($ARGV[1]); # id for the First uClust
chomp ($ARGV[2]); # id for the Second uClust

## Check 1 of 2
my $n_args=@ARGV;
if ($n_args != 3) {
	print "You entered $n_args argument(s):\n";
	foreach my $arg(@ARGV){
		print $arg."\t";
	}
	print "\nThis program takes 3 arguments: 1 input file, percent identity for first iteration, percent identity for secind iteration\n"; 
	print "\n USAGE: perl iClust.pl <INPUT> Pid1 Pid2\n";
	print "where,\nINPUT= A file that contains a list  of names of all genomes in fasta format.\n";
	print "Pid1 = id for the First uClust\n";
	print "Pid2 = id for the Second uClust\n";
	print "This script will now exit.\nPlease check the arguments and try again.\n";
    exit;
}
my (@genomeList0);
my $nG0=0;

open (ALL, "$ARGV[0]") || die "ERROR: $ARGV[0]\n".$!;

while (my $genomeName0=<ALL>){
	unless ($genomeName0=~ /^\s+|^\#/){
		chomp $genomeName0;
		$nG0++;
		my $edGenomeName0=prepFasta($genomeName0, $nG0);
		push(@genomeList0, $edGenomeName0);
	}
}
close ALL;

## Check 2 of 2
if ($nG0<2){
	print "This script requires that the number of Genomes entered be at least 2.\n You entered $nG0.\n";
	print "Please check the input file for spaces or special characters before filenames and make sure each file name has been entered in a seperate line\n";
	print "This script will now exit\n";
	exit;
}

## RANDOMIZE ##
my (@edGL0, %noRep0);
for (my $i0=0; $i0<$nG0; $i0++){
	AGAIN:
	my $r0=int(rand($nG0));
	if ($noRep0{$r0}){
		goto AGAIN;
	}
	else{
		$noRep0{$r0}++;
	}
	my $editedGenome0=$genomeList0[$i0]; ## To randomize: replace '$i0' with '$r0'
	push (@edGL0, $editedGenome0);
}
undef %noRep0;

my $concatenate0=join(" ", @edGL0);
my $combined0="allCombined.faa";
print "Concatenating ".@edGL0." Files...\n";
system("cat ".$concatenate0." > ".$combined0);

system("mkdir Results");
system("mkdir Summary");
system("mkdir ClusterFastas");
system("mkdir Clusters");

my %sequences0=allSeqs($combined0);
my $go2nxt0=0; #anything but 1
my $id0=$ARGV[1];
my $level0=0;

magic($combined0,$id0,$go2nxt0,$nG0,\@edGL0,\%sequences0, $level0);

system ("mv result*.tsv ./Results/");
system ("mv *.log ./Summary/");
system ("mv Cluster*.faa ./ClusterFastas/");
system ("mv *.clust ./Clusters/");
# end timer
my $end=new Benchmark;
my $time=timediff($end,$start);
print "The Process took ". timestr($time, 'all')."\n";
