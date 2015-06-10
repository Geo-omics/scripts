use Bio::DB::EUtilities;
use strict;
 
my @ids;
open (IN, $ARGV[0]) || die "[error] $ARGV[0] : $!\n";
while (my $line=<IN>){
	next if $line=~ m/^#/;
	chomp $line;
	$line=~ s/\r//;
	next unless $line;

	push(@ids, $line);
}
 
my $factory = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                       -email => 'sunitj@umich.edu',
                                       -db    => 'protein',
                                       -id    => \@ids);
 
open (OUT, ">".$ARGV[1]);
while (my $ds = $factory->next_DocSum) {
	my $id=$ds->get_id;
    print OUT $id."\t";
    # flattened mode
    while (my $item = $ds->next_Item('flattened'))  {
        # not all Items have content, so need to check...
		if ($item->get_content){    
	    	my $name= $item->get_name;
			my $content= $item->get_content;
			print OUT $name."\t".$content;
		}
    }
	print OUT "\n";
}

