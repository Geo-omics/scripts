#!/bin/bash

set -e

echo "GeneClusters: $(grep -v "^>" Overview.geneclusters.txt | wc -l)"
echo "smcogs: $(grep -c "^>>" Overview.smcogs.txt)"

ls */structures/* | cut -f 1 -d "/" | sort -u > structures.list
echo "Structures: `wc -l structures.list`"
