#!/usr/local/bin/perl
use Cwd;
use Switch;
### Version 1.1
### 160508 - INCLUDES the "HARD WAY"!!!! :)
### 190508 - Added prediction of multi-core machines to use PTHREADS version
### 190508 - Split "hard way" into three steps...
### 260608 - Fixed bad comparison to get file/directory - caused problems selecting file/dir #0... 
### 260608 - Moved bootstrap # selection to before BKL tree option - allows you to walk away from analyses when running.

$working_directory = getcwd;

print `clear`, "\n";
&detect_multi_core;
&menu;

sub detect_multi_core {

    $cmd = "grep cores /proc/cpuinfo > out.txt";
    system($cmd);

    open( CORES, "<out.txt" );

    $line = <CORES>;
    $re1  = '.*?';
    $re2  = '(\\d+)';
    $re   = $re1 . $re2;
    if ( $line =~ m/$re/is ) {
        $int1 = $1;

    }
    close(CORES);
    unlink("out.txt");
    if ( $int1 >= 2 ) {
        print "You have $int1 cores.\nUsing the RAxML-PTHREADS version.\n";
        $raxml_ver = "raxmlHPC-PTHREADS -T $int1 ";
    }
    else {
        print
"You do not have a multi-core processor or we cannot identify more than 1 core.\nUsing the standard RAxML version.\n";
        $raxml_ver = "raxmlHPC ";
    }

}

sub menu {
    print "**************************************\n";
    print "* An easy to use 'wrapper' for RAxML *\n";
    print "**************************************\n";
    print "1) Fasta to Phylip Conversion\n\n";
    print "2) Predict Substitution Model\n";
    print "3) Predict Substitution Model - Partitioned\n\n";
    print "4) 'Fast & Easy' Bootstrap Search\n";
    print "5) 'Fast & Easy' BS and ML Search\n\n";
    print "6) Best-Known Likelihood (BKL) Tree\n";
    print "7) Bootstrapping\n";
    print "8) Obtain Confidence Values\n\n";
    print "9) 'Hard & Slow' Options 6 - 8\n\n";
    print "Q) Quit\n\n";

    print ">:";

    chomp( $menu_choice = <STDIN> );

    switch ($menu_choice) {

        case "1" {

            &get_sequence_directory;
            &directory_menu;
            $file_type = "fasta";
            &get_file_names;
            &file_menu;
            $in_fasta_file = $menu_file;
            &fasta_to_phylip;
        }
        case "2" {
            &get_sequence_directory;
            &directory_menu;
            $file_type = "phylip";
            &get_file_names;
            &file_menu;
            $in_phylip_file = $menu_file;

            if ( -e $log_directory && -d $log_directory ) {

                $cmd = "cp $sequence_directory\/$in_phylip_file $log_directory";
                system($cmd);
                $alignmentName = "$log_directory\/$in_phylip_file";
                chdir($log_directory);
                &substitution_model_prediction;
            }
            else {
                mkdir( "$log_directory", 0755 );
                $cmd = "cp $sequence_directory\/$in_phylip_file $log_directory";
                system($cmd);
                $alignmentName = "$log_directory\/$in_phylip_file";
                chdir($log_directory);
                &substitution_model_prediction;
            }
            chdir($working_directory);
        }
        case "3" {
            &get_sequence_directory;
            &directory_menu;
            $file_type = "phylip";
            &get_file_names;
            &file_menu;
            $in_phylip_file = $menu_file;
            $alignmentName  = "$log_directory\/$in_phylip_file";
            $file_type      = "part";
            &get_file_names;
            &file_menu;
            $partition = $menu_file;

            if ( -e $log_directory && -d $log_directory ) {

                $cmd = "cp $sequence_directory\/$in_phylip_file $log_directory";
                system($cmd);

                $cmd = "cp $sequence_directory\/$partition $log_directory";

                system($cmd);
                chdir($log_directory);
                &partition_model_prediction;
            }
            else {
                mkdir( "$log_directory", 0755 );
                $cmd = "cp $sequence_directory\/$in_phylip_file $log_directory";
                system($cmd);

                $cmd = "cp $sequence_directory\/$partition $log_directory";

                system($cmd);
                chdir($log_directory);
                &partition_model_prediction;
            }
            chdir($working_directory);

        }
        case "4" {
            &get_sequence_directory;
            &directory_menu;
            $file_type = "phylip";
            &get_file_names;
            &file_menu;
            $in_phylip_file = $menu_file;
            $alignmentName  = "$sequence_directory\/$menu_file";
            chdir($sequence_directory);
            &bootstrap_search;
            chdir($working_directory);
        }
        case "5" {
            &get_sequence_directory;
            &directory_menu;
            $file_type = "phylip";
            &get_file_names;
            &file_menu;
            $in_phylip_file = $menu_file;
            $alignmentName  = "$sequence_directory\/$menu_file";
            chdir($sequence_directory);
            &bootstrap_ml_search;
            chdir($working_directory);
        }
        case "6" {
            &get_sequence_directory;
            &directory_menu;
            $file_type = "phylip";
            &get_file_names;
            &file_menu;
            $in_phylip_file = $menu_file;
            $alignmentName  = "$sequence_directory\/$menu_file";
            chdir($sequence_directory);
            &bkl_tree;
            chdir($working_directory);
        }
        case "7" {
            &get_sequence_directory;
            &directory_menu;
            $file_type = "phylip";
            &get_file_names;
            &file_menu;
            $in_phylip_file = $menu_file;
            $alignmentName  = "$sequence_directory\/$menu_file";
            chdir($sequence_directory);
            print "\nPlease enter the number of bootstrap runs (default: 100)\n";
    	    print ">:";
            chomp( $num_runs = <STDIN> );
            if ( $num_runs le "0" ) {
           	print "\nUsing default value.\n";
           	$num_runs = "100";
            }
            &bootstrapping;
            chdir($working_directory);

        }
        case "8" {
            &get_sequence_directory;
            &directory_menu;
            $file_type = "phylip";
            &get_file_names;
            &file_menu;
            $in_phylip_file = $menu_file;
            $alignmentName  = "$sequence_directory\/$menu_file";
            chdir($sequence_directory);
            $bkl_directory = "$sequence_directory/bkl";
            $bs_directory  = "$sequence_directory/bs";
            &obtain_confidence_values;
            chdir($working_directory);
        }
        case "9" {
            &get_sequence_directory;
            &directory_menu;
            $file_type = "phylip";
            &get_file_names;
            &file_menu;
            $in_phylip_file = $menu_file;
            $alignmentName  = "$sequence_directory\/$menu_file";
            chdir($sequence_directory);
            print "\nPlease enter the number of bootstrap runs (default: 100)\n";
    	    print ">:";
    	    chomp( $num_runs = <STDIN> );
   	    if ( $num_runs le "0" ) {
     	 	   print "\nUsing default value.\n";
            	   $num_runs = "100";
   	    }
            &bkl_tree;
            &bootstrapping;
            &obtain_confidence_values;
            chdir($working_directory);
        }
        case /Q/i {
            print `clear`, "\n";
            print "Goodbye...\n";
            last;
        }
        else {
            print `clear`, "\n";
            print "Incorrect Selection - Please Try Again...\n\n";
            &menu;
        }
    }

}

sub get_sequence_directory {

    @all_files = glob("*");
    @folders;
    foreach my $item (@all_files) {
        if ( -d $item ) {    #Put all folders into array
            push( @folders, $item );
        }
    }
    $folder_num = @folders;
}

sub directory_menu {

    print "Directory File Menu\n";
    print "*******************\n";
    for ( $i = 0 ; $i < $folder_num ; $i++ ) {
        print "$i) $folders[$i]\n";
    }

    print "Please enter the number for the directory where your sequences are located.\n";
    print ">:";

    chomp( $menu_choice = <STDIN> );
    if (   $menu_choice ge $folder_num
        || $menu_choice lt "0"
        || $menu_choice eq m/[a-z]/ )
    {
        print `clear`, "\n";
        print "\nIncorrect Menu Choice!\n\n";
        $menu_choice = "";
        &directory_menu;

    }
    else {
        $sequence_directory = "$working_directory/$folders[$menu_choice]";
        $log_directory      = "$sequence_directory/logs";
    }
}

sub get_file_names {
    if ( -e $sequence_directory && -d $sequence_directory ) {
        @file_names = <$sequence_directory/*.$file_type>;
        foreach $file (@file_names) {

            $file = reverse $file;
            @fn   = split( /\//, $file );
            $file = reverse $fn[0];
        }
        $file_num = @file_names;
    }
    else {
        print `clear`, "\n";
        print "\nThere is no sequence directory, please create one.\n";
        last;
    }
}

sub file_menu {

    if ( $file_num < 1 ) {

        print `clear`, "\n";
        print
"\nThere are no $file_type sequences in the directory, please make sure you have selected the right directory.\n";
        last;

    }
    else {
        print "$file_type File Menu\n";
        print "******************\n";
        for ( $i = 0 ; $i < $file_num ; $i++ ) {
            print "$i) $file_names[$i]\n";
        }

        print "Please enter the number for your sequence file.\n";
        print ">:";

        chomp( $menu_choice = <STDIN> );
        if (   $menu_choice ge $file_num
            || $menu_choice lt 0
            || $menu_choice eq m/[a-z]/ )
        {
            print `clear`, "\n";
            print "\nIncorrect Menu Choice!\n\n";
            $menu_choice = "";
            &file_menu;

        }
        else {
            $menu_file = $file_names[$menu_choice];
        }
    }

}

## Adapted from code, without permission, from the main RAxML website.
sub substitution_model_prediction {

    #print "\nStart = " . scalar localtime(time) . "\n";
    @AA_Models = (
        "DAYHOFF",   "DCMUT",  "JTT",      "MTREV",  "WAG",      "RTREV",
        "CPREV",     "VT",     "BLOSUM62", "MTMAM",  "DAYHOFFF", "DCMUTF",
        "JTTF",      "MTREVF", "WAGF",     "RTREVF", "CPREVF",   "VTF",
        "BLOSUM62F", "MTMAMF"
    );
    print "\nBuilding initial parsimony tree.\n";
    $cmd =
        $raxml_ver
      . "-y -m PROTCATJTT -s "
      . $alignmentName
      . " -n ST_"
      . $in_phylip_file
      . " \> ST_"
      . $in_phylip_file . "_out";

    #print($cmd);
    system($cmd);
    $numberOfModels = @AA_Models;
    print "Starting to optimise model parameters.\n";
    for ( $i = 0 ; $i < $numberOfModels ; $i++ ) {
        $aa = "PROTGAMMA" . $AA_Models[$i];
        $cmd =
            $raxml_ver
          . "-f e -m "
          . $aa . " -s "
          . $alignmentName
          . " -t RAxML_parsimonyTree.ST_"
          . $in_phylip_file . " -n "
          . $AA_Models[$i] . "_"
          . $in_phylip_file
          . "_EVAL \> "
          . $AA_Models[$i] . "_"
          . $in_phylip_file
          . "_EVAL.out\n";

        #print($cmd);
        system($cmd);
        print "Optimising model $i of 20: $aa\n";
    }
    print "Completed optimising models.\n";
    print "Generating logs";
    for ( $i = 0 ; $i < $numberOfModels ; $i++ ) {
        $logFileName =
          "RAxML_log." . $AA_Models[$i] . "_" . $in_phylip_file . "_EVAL";

        #print $logFileName."\n";
        $lh[$i] = getLH($logFileName);
        print ".";
    }

    $bestLH = -1.0E300;
    $bestI  = -1;

    for ( $i = 0 ; $i < $numberOfModels ; $i++ ) {

        #print "\nModel: ".$AA_Models[$i]." LH: ". $lh[$i]."\n";
        if ( $lh[$i] > $bestLH ) {
            $bestLH = $lh[$i];
            $bestI  = $i;
        }
    }
    open( OUTFILE, ">>$sequence_directory/best_model_$in_phylip_file.txt" );

    print "\n\nBest Model : " . $AA_Models[$bestI] . "\n\n";
    print OUTFILE "\n\nBest Model : " . $AA_Models[$bestI] . "\n\n";
    close(OUTFILE);

    $bestModel = $AA_Models[$bestI];

    #print "\nEnd = " . scalar localtime(time) . "\n";

}

## Adapted from code, without permission, from the main RAxML website.
sub partition_model_prediction {

    print "How many parts to your alignment? (0-9)\n";
    print ">:";
    chomp( $count = <STDIN> );
    print "Splitting up multi-gene alignment\n";
    $cmd =
        $raxml_ver
      . "-f s -m PROTCATJTT -s "
      . $alignmentName . " -q "
      . $partition
      . " -n SPLIT_"
      . $in_phylip_file
      . " \> SPLIT_"
      . $in_phylip_file . "_out";
    $alignment_temp = $alignmentName;
    $phylip_temp    = $in_phylip_file;

    system($cmd);
    for ( $x = 0 ; $x < $count ; $x++ ) {

        #print "\n\n$in_phylip_file . ".GENE." . $x\n\n";
        print "PARTITION: " . $x . "\n";
        $alignmentName  = $alignmentName . ".GENE." . $x;
        $in_phylip_file = $in_phylip_file . ".GENE." . $x;
        &substitution_model_prediction;
        $alignmentName  = $alignment_temp;
        $in_phylip_file = $phylip_temp;
    }

}

## Adapted from code, without permission, from the main RAxML website.
sub getLH {
    my $fileID = $_[0];
    open( CPF, $fileID );
    my @lines = <CPF>;
    close(CPF);
    my $numIT  = @lines;
    my $lastLH = pop(@lines);
    my $k      = index( $lastLH, '-' );
    my $LH     = substr( $lastLH, $k );
    return $LH;
}

## Adapted from code, without permission, from the main RAxML website.
sub getTIME {
    my $fileID = $_[0];
    open( CPF, $fileID );
    my @lines = <CPF>;
    close(CPF);
    my $numIT  = @lines;
    my $lastLH = pop(@lines);
    my $k      = index( $lastLH, '-' );
    my $TIME   = substr( $lastLH, 0, $k - 1 );
    return $TIME;
}

## Based on a script available from the author of RAxML
#  from http://icwww.epfl.ch/~stamatak/index-Dateien/Page443.htm
sub fasta_to_phylip {

    $out_fasta_file =
      " > " . $sequence_directory . "/" . $in_fasta_file . ".phylip";
    open( INFILE, " < " . "$sequence_directory/$in_fasta_file" );

    $taxa = -1;

    while ( $line = <INFILE> ) {
        if ( $line =~ />/ ) {
            $taxa++;
            $name = $line;
            $name =~ s/[\s+|\(|\)|\,|;]//g;
            $name =~ s/,//g;
            $name =~ s/>//g;
            $taxonNames[$taxa] = $name;
        }
        else {
            $seq = $line;
            $seq =~ s/\s+//g;
            $sequences[$taxa] = $sequences[$taxa] . $seq;
        }
    }

    close(INFILE);

    $s  = $taxa + 1;
    $bp = length( $sequences[0] );

    print "\nConverting";
    open( OUTFILE, $out_fasta_file );
    print OUTFILE $s . " " . $bp . "\n";

    for ( $i = 0 ; $i <= $taxa ; $i++ ) {
        print OUTFILE $taxonNames[$i] . " " . $sequences[$i] . "\n";
        print ".";
    }
    print "\nFile Converted$out_fasta_file\n";
}

sub model_selection {

    print "DNA Data or AA Data?\n1) DNA\n2) Amino Acid\n>:";
    chomp( $dna_or_aa = <STDIN> );
    if ( $dna_or_aa eq "1" ) {
        $model = "GTRMIX";
        print "\nEstimate proportion of invariable sites?\n1) Yes\n2) No\n>:";
        chomp( $base_freq = <STDIN> );
        if ( $base_freq eq "1" ) {
            $model = $model . "I";
        }
        $is_model_set = "true";
    }
    elsif ( $dna_or_aa eq "2" ) {
        $model = "PROTMIX";
        print
"\nEstimate proportion of invariable sites? (Model + I)\n1) Yes\n2) No\n>:";
        chomp( $base_freq = <STDIN> );
        if ( $base_freq eq "1" ) {
            $model = $model . "I";
        }
        @aa_models = (
            "DAYHOFF", "DCMUT", "JTT",      "MTREV", "WAG", "RTREV",
            "CPREV",   "VT",    "BLOSUM62", "MTMAM", "GTR"
        );
        $aa_num = @aa_models;
        print "\nWhich Amino Acid Model?\n";
        for ( $i = 0 ; $i < $aa_num ; $i++ ) {
            print "$i) $aa_models[$i]\n";
        }
        print ">:";
        chomp( $input_model = <STDIN> );
        $model = $model . "$aa_models[$input_model]";
        print
          "\nUse empirical base frequencies? (Model + F)\n1) Yes\n2) No\n>:";
        chomp( $base_freq = <STDIN> );
        if ( $base_freq eq "1" ) {
            $model = $model . "F";
        }
        $is_model_set = "true";
    }
    else {
        print `clear`, "\n";
        print "Incorrect selection. Please try again!";
        &model_selection;
    }

}

sub bootstrap_search {

    print "Please enter the number of runs (default: 100)\n";
    print ">:";
    chomp( $num_runs = <STDIN> );
    if ( $num_runs le "0" ) {
        print "\nUsing default value.\n";
        $num_runs = "100";
    }

    #
    print
"Please enter the substituion model (default: GTRMIX+Model or PROTMIX+Model)\n";
    &model_selection;

    $cmd =
        $raxml_ver
      . "-x 12345 -p 12345 -N "
      . $num_runs . " -m "
      . $model . " -s "
      . $alignmentName . " -n "
      . $in_phylip_file . ".out";

    #print "\n\n$cmd\n\n";
    system($cmd);

}

sub bootstrap_ml_search {

    print "Please enter the number of runs (default: 100)\n";
    print ">:";
    chomp( $num_runs = <STDIN> );
    if ( $num_runs le "0" ) {
        print "\nUsing default value.\n";
        $num_runs = "100";
    }

    #
    print
"Please enter the substituion model (default: GTRMIX+Model or PROTMIX+Model)\n";
    &model_selection;

    $cmd =
        $raxml_ver
      . "-f a -x 12345 -p 12345 -N "
      . $num_runs . " -m "
      . $model . " -s "
      . $alignmentName . " -n "
      . $in_phylip_file . ".out";

    #print "\n\n$cmd\n\n";
    system($cmd);

}

sub bkl_tree {
    print "Finding the Best-Known Likelihood (BKL)\n";
    if ( $is_model_set eq "true" ) {
        print "\nModel already selected = $model";
    }
    else {
        &model_selection;
    }
    print "\nPlease enter number of inferences (default/min: 2)\n";
    print ">:";
    chomp( $num_infs = <STDIN> );
    if ( $num_infs <= "2" ) {
        print "\nUsing default value of 2.\n";
        $num_infs = "2";
    }
    $bkl_directory = "$sequence_directory/bkl";

    ## This should be obsolete once I implement a 'run' system as per fdfBLAST...
    if ( -e $bkl_directory && -d $bkl_directory ) {

        print "\nYou appear to have already run this...\n";
    }
    else {
        mkdir( "$bkl_directory", 0755 );
        chdir($bkl_directory);

        $cmd =
            $raxml_ver
          . "-f d -m "
          . $model . " -s "
          . $alignmentName . " -N "
          . $num_infs . " -n "
          . $in_phylip_file . ".out";

        #print "\n\n$cmd\n\n";
        system($cmd);
    }
    chdir($working_directory);

}

sub bootstrapping {
    print "\nBootstrapping...\n";
    if ( $is_model_set eq "true" ) {
        print "\nModel already selected = $model";
    }
    else {
        &model_selection;
    }
    $bs_directory = "$sequence_directory/bs";

    ## This should be obsolete once/if I implement a 'run' system as per fdfBLAST...
    if ( -e $bs_directory && -d $bs_directory ) {

        print "\nYou appear to have already run this...\n";
    }
    else {
        mkdir( "$bs_directory", 0755 );
        chdir($bs_directory);
        $cmd =
            $raxml_ver
          . "-f d -m "
          . $model . " -s "
          . $alignmentName . " -N "
          . $num_runs . " -b "
          . "12345 -n "
          . $in_phylip_file . ".boot";

        #print "\n\n$cmd\n\n";
        system($cmd);
    }
    chdir($working_directory);
}

sub obtain_confidence_values {

    $infofile = "$bkl_directory/RAxML_info.$in_phylip_file.out";
    open( INFOFILE, " < " . "$infofile" );

    foreach $line (<INFOFILE>) {
        chomp($line);
        $txt = $line;
        $re1 = '.*?';       # Non-greedy match on filler
        $re2 = '(\\d+)';    # Integer Number 1
        $re3 = '(:)';       # Single Character 1

        $re = $re1 . $re2 . $re3;
        if ( $txt =~ m/$re/is ) {
            $int1 = $1;
            print "\nBest Scoring Tree found in RUN $int1 \n";
        }
    }

    if ( $is_model_set eq "true" ) {
        print "\nModel already selected = $model";
    }
    else {
        &model_selection;
    }

    ##
    $tree_directory = "$sequence_directory/tree";
    mkdir( "$tree_directory", 0755 );
    chdir($bs_directory);
    $cmd = "cp RAxML_bootstrap.$in_phylip_file.boot $tree_directory";
    system($cmd);
    chdir($bkl_directory);
    $cmd = "cp RAxML_result.$in_phylip_file.out.RUN.$int1 $tree_directory";
    system($cmd);
    chdir($tree_directory);

    print "\nObtaining Confidence Values...\n";
    $cmd =
        $raxml_ver
      . "-f b -m "
      . $model . " -s "
      . $alignmentName . " -z "
      . "RAxML_bootstrap."
      . $in_phylip_file
      . ".boot -t "
      . "RAxML_result."
      . $in_phylip_file
      . ".out.RUN."
      . $int1
      . " -n final."
      . $in_phylip_file . ".tree";

    #print "\n\n$cmd\n\n";
    system($cmd);
    unlink("RAxML_bootstrap.$in_phylip_file.boot");
    unlink("RAxML_result.$in_phylip_file.out.RUN.$int1");

}
