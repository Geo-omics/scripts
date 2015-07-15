#! /bin/sh
if [ -z $1 ]; then echo "command failed: Give a file name with list"; exit; fi
if [ ! -s $1 ]; then echo "$1 does not exist"; exit; fi
for i in $(grep "^" $1); do
    if [ ! -s $i ]; then echo "$i is empty; skipping..."; continue; fi
    OUT="$i.out"
    # replace the following line with the desired command and $i as input and $OUT as output
	echo "IN: $i, OUT:$OUT"
done
