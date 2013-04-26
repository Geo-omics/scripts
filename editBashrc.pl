#! /usr/bin/perl

use strict;
#use File::HomeDir;

my $home=`echo \$HOME`;
chomp($home);

my $bashProfile = $home . "/.bash_profile";
my $bpBackup = $home ."/bashProfile.bkp.txt";
my $add2Path="/geomicro/data1/COMMON/bin";

system("cp $bashProfile $bpBackup");

if (-e $bpBackup){
	system("rm $bashProfile");

	open(NBP, ">".$bashProfile)|| die "[ERROR:$0] NBP: $!\n";
	open(BP, $bpBackup)|| die "[ERROR:$0] BP: $!\n";
	my $edited=0;
	while(my $line=<BP>){
		if (($line=~ /^PATH\s*\=/) || ($line=~ /export PATH\s*\=/)){
			my ($var, $currentPath)=split(/\=/, $line);
			print NBP $var."=".$add2Path.":".$currentPath;
			print $var."=".$add2Path.":".$currentPath;
			$edited++;
		}
		else{
			print NBP $line;
		}
	}
	close BP;
#	print $edited."\n";

	if ($edited==0){
		print NBP "# User specific environment and startup programs\n";
		print NBP "export PATH=$add2Path:\$PATH\n"
	}
	close NBP;
	close BP;

	print "You're all set! To have the new settings take effect, please log out of the current server and back in.\n";
}
else{
	die "[ERROR] Could not create backup file. The script can't proceed further without creating a backup. Get Sunit.\n";
}


