#!/usr/bin/perl -w

use strict;

my $delay= $ARGV[0] ? $ARGV[0] : 2 ;

my $count;
for (1..100){
	$count++;
	sleep($delay);
}
