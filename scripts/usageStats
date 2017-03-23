#! /usr/bin/perl

=head1 DESCRIPTION

	usageStats -- notes the CPU and Memory usage of a process at each defined time interval. Presents it in a usable format.

=head2 USAGE

	perl usageStats.pl

=head3 Options

	-t	0 or more search terms;		default: velveth velvetg mVelvetPipe_paired.sh meta-velvetg
	-i time interval in seconds;	default: 10
	-o output file; 		default: processID.log
	-e [flag] email the results file once the jobs have finished.

=head1 Suggestions/Comments/Feedback/Beer

	Sunit Jain (sunitj [AT] umich [DOT] edu)
	May 2012	

=cut

use strict;
use Getopt::Long;

my $thisScript=lc($0);
my @searchTerms;
my $intv=10;
my $out=$$.".log";
my $mail;
GetOptions(
	't|search:s{,}'=>\@searchTerms, # " :{,} " 0 or more terms for search; " ={,} " 1 or more terms for search
	'i|interval:i'=>\$intv,
	'o|out:s'=>\$out,
	'e|email'=>\$mail,
	'h'=>sub{system('perldoc', $0); exit;},
);

my $start=`date`;
@searchTerms= qw(velveth velvetg meta-velvetg oases) if scalar(@searchTerms)==0;

open (OUT, ">".$out) || die "[ERROR: $0] $!\n";
print OUT "# Start Time:\t$start\n";

my (@master, $end);
while (1){
	my $m= &getSnapshot;
	my @gotInfo=@{$m};
	if (scalar(@gotInfo) >0){
		push (@master, @gotInfo);
		sleep $intv;
	}
	else{
		$end=`date`;
		last;
	}
}
chomp ($start, $end);

print OUT "# End Time:\t$end\n";
print OUT "# Monitor Interval:\t $intv\n";
print OUT "# PID\tCommand\t\%Mem\t\%CPU\tUser\n";
foreach my $m(@master){
	my($user, $pid, $percCPU, $percMem, $vsz, $rss, $tty, $status, $start, $cpuTime, @command)= split(/\s+/, $m);
	my $comm=join(" ", @command);
	print OUT "$pid\t$comm\t$percMem\t$percCPU\t$user\n";
}
close OUT;

if ($mail){
	my $uname=`whoami`;
	chomp $uname;
	my $email=$uname."\@umich.edu";

	system("mail -s 'Job Completed' $email < $out");
}

sub getSnapshot{
	my @thisInstant;
	foreach my $s(@searchTerms){
		next if (lc($s) eq "perl"); # else, will cause the script to stay in an infinite loop;
		next if (lc($s) eq $thisScript);
		chomp $s;
#		print $s."\n";
		my @arr;
		if ($s=~ /^\d+$/){
			my @arr=`ps u -p $s`;
			my $shifted=shift @arr;
			print scalar(@arr)."\n";
		}
		else{
			@arr=`ps u -C $s`;
			shift @arr;
		}
		push (@thisInstant, @arr);
	}
	return (\@thisInstant);
}

