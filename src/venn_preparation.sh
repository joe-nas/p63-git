#! /bin/bash
# create folder shortcuts
UNIONFOLDER="../results/tables/venn_diagrams/h3k4me1_h3k4me3_h3k27ac"
ENCODEFOLDER="../data/hg19/encode/peaks/"
mkdir -p $UNIONFOLDER; # make destination folder
# finding files, concatenate, sort, merge/collapse, identifyer
find $ENCODEFOLDER -iname *h3k4me1* -o -iname *h3k4me3* -o -iname *h3k27ac* | xargs cat | bedtools sort -i stdin |bedtools merge -i stdin| awk -F t '{print $1$2$3"\tpeak_"NR}' > $UNIONFOLDER/union_h3k4me1_h3k4me3_h3k27ac.bed;
# overlap files with union_h3k4me1_h3k4me3_h3k27ac
UNIONFILE="$UNIONFOLDER/union_h3k4me1_h3k4me3_h3k27ac.bed" # shortcut to unionfile
# intersect unionfile with histone files
for i in $(find $ENCODEFOLDER -iname *h3k4me1* -o -iname *h3k4me3* -o -iname *h3k27ac*); do 
  HISTONE=$(basename $i);
  bedtools intersect -wa -a $UNIONFILE -b $i | cut -f 4 >> $UNIONFOLDER/intersect_union_$HISTONE;
done;
