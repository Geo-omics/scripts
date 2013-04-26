#### Anders Andersson 2007 ####

# Infiles (may be changed):
$blast_out_file = "UBAnlCONTIGSvsAMD_combined_scaffs2.fasta";

#-------------- main program --------
&extract_blast_hits; 

#----------- sub routines ----------
sub extract_blast_hits {
  local(@fields);
  local(@sub_fields);
  local($score);
  local($query);
  local($hit);
  local($a_length); 
  local($ratio);
  local($a_id);
  local($a_pos);
  print"query\tq_length\thit\th_length\ta_length\ta_id\ta_pos\ta_score\te_value\n";
  open (INFILE, $blast_out_file) || die ("could not open !");
  $read = 0;
  while (<INFILE>) {
    #print"$_";
    $read++;
    chomp;
    if (/^Query= /) {
        @fields = split(/\s+/);
        $query = $fields[1]." ".$fields[2];
    }
    if (/(\d+ letters)/) {
        ($q_length) = $_ =~ /(\d+)/;
    }
    if (/^>\w/) {
        $hit = $_;
    }
    if (/Length = /) {
        @fields = split(/\s+/);
        $h_length = $fields[3];
    }
    if (/Score = /) {
        @fields = split(/\s+/);
        $a_score = $fields[3];
        $e_value = $fields[8];  
    }
    if (/Identities = /) {
        @fields = split(/,/);
        ($ratio) = $fields[0] =~ /((\d)+\/(\d)+)/;
        @sub_fields = split(/\//,$ratio);
        $a_length = $sub_fields[1];
        $a_id = $sub_fields[0];

        ($ratio) = $fields[1] =~ /((\d)+\/(\d)+)/;
        @sub_fields = split(/\//,$ratio);
        $a_pos = $sub_fields[0];
        
        @fields = split(/\[/,$hit);
        $hit = "[".$fields[-1];
        #print"$hit\n";
        $hits{$hit}++;
        print"$query\t$q_length\t$hit\t$h_length\t$a_length\t$a_id\t$a_pos\t$a_score\t$e_value\n";
    }
  }
  close (INFILE);
  #foreach $hit (keys %hits) {
  #  print"$hit\t$hits{$hit}\n";
  #}
  #$hit = (keys %hits);
  #print"$hit\n";
}


