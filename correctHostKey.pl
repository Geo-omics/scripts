#! /usr/bin/perl

use strict;
#use File::HomeDir;
 
my $home=`echo \$HOME`;
chomp($home);

my $known_host = $home . "/.ssh/known_hosts";
my $known_bkp = $home ."/known_host.bkp.txt";

print "#Please enter the line number that has the 'Offending Key' Or press 'enter', if you don't know: ";
my $offendingKey=<STDIN>;
chomp $offendingKey;

my $errorHelp="simply try to SSH into the server and in the error message look for the 'Offending key'. If you don't see the offending key copy the error message and mail it to sunitj [AT] umich [DOT] edu\n";

## Just Being Thorough ##
die "#[ERROR: known_host file]$!\n" unless (-e $known_host);

die "#Offending key required.\n To get it, ".$errorHelp if (! $offendingKey);

my $maxOffKey=`wc -l $known_host | cut -d " " -f 1`;
chomp ($maxOffKey);
die "#Invalid Offending key.\n Your Offending key must be between '1-$maxOffKey'. To get it, ".$errorHelp if ($offendingKey > $maxOffKey);

die "#Invalid Offending key.\n To get it, ".$errorHelp if ($offendingKey==0);

warn "\n#Please confirm the line number or this will SERIOUSLY mess up your server settings!\n";
print "#Are you sure this \( $offendingKey \) is the 'Offending Key' line number?[Y/N]: ";
my $sure=<STDIN>;
chomp ($sure);
print "\n";

$sure= "y" if lc($sure) eq "yes";

lc($sure) eq "y" ? print "#Ok then, here goes. I warned you!\n" : die "#To confirm, ".$errorHelp;

##

`cp $known_host $known_bkp`;

if (-e $known_bkp){
	print "#In case something goes wrong. The back of $known_host file is here $known_bkp\n";
	print "rm -f $known_host\n";
	`rm -f $known_host`;

	open(BKP, $known_bkp) || die "[ERROR: BKP] $0: $!\n";
	open(KH, ">".$known_host)|| die "[ERROR: KH] $0: $!\n";
	while(my $line=<BKP>){
		print KH $line unless $.==$offendingKey;
	}
	close KH;
	close BKP;

	print "#All Done! Try logging into the server again.\n";
}
else{
	print STDERR "#Script Borked! Get Sunit sunitj [AT] umich [DOT] edu.\n";
}
