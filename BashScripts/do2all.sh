#function do2all() {
if [ -z $1 ]; then echo "command failed: Give a directroy name with files"; return 1; fi
if [ ! -d $1 ]; then echo "$1 does not exist"; return  1; fi
for i in $1/*; do
        OUT="$i.out"
        if [ ! -s $i ]; then echo "$i is empty; skipping..."; continue; fi
		echo "perl name.pl $i $OUT"
done
#}
