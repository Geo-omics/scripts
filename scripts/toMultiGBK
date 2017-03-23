#!/usr/local/bin/perl

=head1 DESCRIPTION

toMultiGBK.pl -- split a multi GBK file to multiple single GBK file.

=head1 USAGE

perl toMultiGBK.pl -gbk file.gbk

=head2 Options

    -gbk        <CHAR>  Input multi GBK file
    -dir -d     <CHAR>  Output directory. [default= SPLIT]
	-version -v	<BOOLEAN>	version of the current script
	-help	-h	<BOOLEAN>	This message.

=head1 Author

Sunit Jain, (Mon Jan 19 16:09:25 EST 2015)
sunitj [AT] umich [DOT] edu

=head1 License

This script is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=head1 Disclaimer

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut

use strict;
use Getopt::Long;
use File::Basename;
use FileHandle;

my $help;
my $version="toMultiGBK.pl\tv0.0.6b";
my ($gbk, $directory);
GetOptions(
    'gbk:s'=>\$gbk,
    'd|dir:s'=>\$directory,
	'v|version'=>sub{print $version."\n"; exit;},
	'h|help'=>sub{system("perldoc $0 \| cat"); exit;},
);
print "\# $version\n";

my $counter=1;

unless ($directory) { $directory="SPLIT"; }
unless(-e $directory or mkdir $directory) {
        die "Unable to create $directory\n";
    }
my $filename=fileparse($gbk, ".gbk");
$/="//\n";

open(GBK, "<".$gbk)|| die $!;
while(my $line=<GBK>){
    my $FH=FileHandle->new;
    open($FH, ">".$directory."/".$filename."_part".$counter.".gbk") || die $!;
    print $FH $line;
    close $FH,
    $counter++;
}
close GBK;





__END__
#open($FH, ">".$directory."/".$filename."_part".$counter.".gbk") || die $!;
open(GBK, "<".$gbk)|| die $!;
while(my $line=<GBK>){
    open($FH, ">".$directory."/".$filename."_part".$counter.".gbk") || die $!;
    unless ($line=~/^\/\//){
        print $FH $line;
    }
    else{
        print $FH $line;
        close $FH;
        $counter++;
        $FH=FileHandle->new;
    }
    
    die "[FATAL] Too many files written!!\n" if ($counter == 10);
}
close GBK;