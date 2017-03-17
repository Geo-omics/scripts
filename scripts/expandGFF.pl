#!/usr/bin/perl

=head1 DESCRIPTION

expandGFF.pl -- Take the GFF file attribute column and split it into multiple columns.

=head1 USAGE

perl expandGFF.pl -gff -out

=head2 Options

    -gff        <FILE>      GFF File Name
    -out        <FILE>      Output File Name
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Thu Aug  6 08:01:10 EDT 2015)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use FileHandle;
use File::Basename;

my $gffFile;
my $outFile;
my $help;
my $version=fileparse($0)."\tv0.0.1b";
GetOptions(
    'gff:s'=>\$gffFile,
    'o|out:s'=>\$outFile,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my @gffHeaders=qw(Scaffold Source Type Start End Score Strand Phase);
my $GFF=FileHandle->new();
my (%data, %uniqCols);
open( $GFF, "<", $gffFile) || die $!;
while(my $line=<$GFF>){
    chomp $line;
    next if ($line=~/^#/);
    next unless $line;
    
    my @columns=split(/\t/, $line);
    my $att=pop(@columns);
    add_columns(\@columns, $.);
    my @attr=split(/;/, $att);
    add_columns(\@attr, $.);
}
close $GFF;

my @sortedColNames=sort{$a<=>$b} keys %uniqCols;

my $OUT=FileHandle->new();
open($OUT, ">", $outFile) or die $!;
my $header;
$header.=$_."\t" foreach(@sortedColNames);
$header=~s/\s+$/\n/;
print $OUT $header;
foreach my $row(keys %data){
    my $line;
    $line.=$data{$row}{$_}."\t" foreach(@sortedColNames);
    $line=~s/\s+$/\n/;
    print $OUT $line;
}
close $OUT;
exit 0;

sub print_columns{
    my $arr=shift;

}
sub add_columns{
    my $arr=shift;
    my $rowNum=shift;
    my $colNum=0;
    foreach my $a(@{$arr}){
        if ($a=~/\=/) {
            my($key,$value)=split(/=/,$a);
            $data{$rowNum}{$key}=$value;
            $uniqCols{$key}++;
        }
        else{
            my $colName=$gffHeaders[$colNum];
            $data{$rowNum}{$colName}=$a;
            $uniqCols{$colName}++;
            $colNum++;
        }
    }
    return;
}


