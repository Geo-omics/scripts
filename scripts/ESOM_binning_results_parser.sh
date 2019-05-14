cat *.conf > esom_scaffolds2bin.tsv
sed '/^#/ d' esom_scaffolds2bin.tsv > esom_scaffolds2bin.cleaned.tsv
awk 'BEGIN{OFS="\t"}{print $2,$1}' esom_scaffolds2bin.cleaned.tsv > esom_scaffolds2bin.tsv
awk 'BEGIN{OFS="\t"}{$2="Bin_"$2; print}' esom_scaffolds2bin.tsv > ESOM_binning_results.txt
perl -pe 's/(?<=\d)_(?=\d)/./g' ESOM_binning_results.txt > ESOM_binning_results.txt.fixed
sed 's/k141\./k141_/g' ESOM_binning_results.txt.fixed > ESOM_binning_results.txt
rm esom_scaffolds2bin.cleaned.tsv
rm esom_scaffolds2bin.tsv
rm ESOM_binning_results.txt.fixed
