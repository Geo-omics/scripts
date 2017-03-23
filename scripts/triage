#! /usr/bin/perl

=head1 DESCRIPTION

	triage.pl - Give it forward and reverse fastq files and this will seperate the sequences based on whether a fwd-rev pair can be found or not.

=head2 USAGE

	triage.pl -fwd forward.fastq -rev reverse.fastq [-n for the new illumina pipeline casava 1.8+ that has a space instead of a '/' in sequence headers]

=head2 OUTPUTS

	Created 4 files, with the processID (pid) as the prefix.
	pid.fp.fastq : forward pair found in reverse
	pid.fs.fastq : forward singleton
	pid.rp.fastq : reverse pair found in forward
	pid.rs.fastq : reverse singleton

=head2 Questions/Suggestions/Feedback

	Contact: Sunit Jain (sunitj [AT] umich [DOT] edu)

=cut

# @DGFC08P1:94:C065LACXX:7:1101:1218:2051 1:N:0:
# @DGFC08P1:94:C065LACXX:7:1101:1218:2051 2:N:0:
use strict;
use Getopt::Long;

my $forward;
my $reverse;
my $fwdP=$$.".fp.fastq";
my $fwdS=$$.".fs.fastq";
my $revP=$$.".rp.fastq";
my $revS=$$.".rs.fastq";
my $new;

GetOptions(
	'fwd:s'=>\$forward,
	'rev:s'=>\$reverse,
	'n|new'=>\$new,
	'h|help'=>sub{system('perldoc', $0); exit;},
);

my $delim= $new ? " " : "\/";

my %fwd;
open(FWD, $forward) || die "[$0] ERROR: $!\n";
while(my $line=<FWD>){
	$line= &polish($line);
	next unless $line;

	if ($line=~ /^@\w+/){
#header
		my ($name, $strand)= &splitHeader($line);
		$fwd{$name}{"Strand"}=$strand;
		$fwd{$name}{"Count"}++;
#sequence
		$line=<FWD>;
		$line= &polish($line);
		$fwd{$name}{"Seq"}=$line;
#qual header
		$line=<FWD>;
#qual score
		$line=<FWD>;
		$line= &polish($line);
		$fwd{$name}{"Qual"}=$line;
	}
	else{
		print "Script Borked! Discuss with Sunit.\n";
		exit;
	}
}
close FWD;

#=head1
open(REV, $reverse) || die "[$0] ERROR: $!\n";
open(FWDP, ">".$fwdP) || die "[$0] ERROR: $!\n";
open(REVP, ">".$revP) || die "[$0] ERROR: $!\n";
open(FWDS, ">".$fwdS) || die "[$0] ERROR: $!\n";
open(REVS, ">".$revS) || die "[$0] ERROR: $!\n";

while(my $line=<REV>){
	$line= &polish($line);
	next unless $line;

	if ($line=~ /^@\w+/){
		my $header=$line;
		my($name, $strand)= &splitHeader($line);#split(/\//, $line);

		if($fwd{$name}){

## Forward
			print FWDP "@".$name.$delim.$fwd{$name}{"Strand"}."\n"; # header
			print FWDP $fwd{$name}{"Seq"}."\n"; # sequence
			print FWDP "+\n"; #.$name.$delim.$fwd{$name}{"Strand"}."\n"; # qual header
			print FWDP $fwd{$name}{"Qual"}."\n"; # qual score

## Reverse
			print REVP "@".$name.$delim.$strand."\n"; # header

			$line=<REV>;
			$line= &polish($line);
			print REVP $line."\n"; # sequence

			$line=<REV>;
			$line= &polish($line);
			print REVP "+\n"; #$line."\n"; # qual header

			$line=<REV>;
			$line= &polish($line);
			print REVP $line."\n"; # qual score

			delete($fwd{$name});
		}
		else{
			print REVS $line."\n"; # header

			$line=<REV>;
			$line= &polish($line);
			print REVS $line."\n"; # sequence

			$line=<REV>;
			$line= &polish($line);
			print REVS $line."\n"; # qual header

			$line=<REV>;
			$line= &polish($line);
			print REVS $line."\n"; # qual score
		}
	}	
}
close REV;
close REVS;
close REVP;
close FWDP;
#=cut

foreach my $name(keys %fwd){
	if ($fwd{$name}{"Count"} > 1){
		print $name."\t".$fwd{$name}{"Count"}."\n";
	}
#=head1
	print FWDS "@".$name.$delim.$fwd{$name}{"Strand"}."\n"; # header
	print FWDS $fwd{$name}{"Seq"}."\n"; # sequence
	print FWDS "+\n"; #.$name.$delim.$fwd{$name}{"Strand"}."\n"; # qual header
	print FWDS $fwd{$name}{"Qual"}."\n"; # qual score
#=cut
}
close FWDS;

sub splitHeader{
	my $header=shift;
	$header=~ s/^@//;
	my (@headerParts)=split(/ /, $header);
	if (scalar(@headerParts)==2){
		return($headerParts[0], $headerParts[1]);
	}
	else{
		@headerParts=split(/\//, $header);
		return($headerParts[0], $headerParts[1]);
	}
}


sub polish{
	my $line=shift;
	chomp($line);
	$line=~ s/\r//;
	$line=~ s/^\s+//;
	$line=~ s/\s+$//;
	return $line;
}

exit;
