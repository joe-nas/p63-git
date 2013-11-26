#! /bin/bash

# 1000 iterations of: 1. p63 shuffle, 2. intersect/count overlap with histone #marks, 3. append count to file.

mkdir -p ../results/tables/p63_exp_vs_ovs_overlap;
#for i in $(seq 1 1000);
#  do for j in $(ls ../data/hg19/encode/peaks/);
#    do echo "iteration $i file $j";
#    bedtools shuffle -chrom -seed $i -i ../data/hg19/p63all_sort_merge_hg.19.bed -g ../data/hg19/hg19.genome | bedtools intersect -a stdin -b ../data/hg19/encode/peaks/$j | wc -l >> ../results/tables/p63_exp_vs_ovs_overlap/p63_exp_$j
#  done;
#done;

# parameter j determines the number of cores used.
/home/jfalck/bin/parallel -j20 --gnu --keep-order --eta 'bedtools shuffle -chrom -seed {1} -i ../data/hg19/p63all_sort_merge_hg.19.bed -g ../data/hg19/hg19.genome | bedtools intersect -a stdin -b ../data/hg19/encode/peaks/{2} | wc -l | xargs echo exp  >> ../results/tables/p63_exp_vs_ovs_overlap/hg19_p63_obs_vs_exp_{2}' ::: $(seq 1 1000) ::: $(ls ../data/hg19/encode/peaks/);

# append the observed overlap to the file
for i in $(ls ../data/hg19/encode/peaks/);
  do bedtools intersect -a ../data/hg19/p63all_sort_merge_hg.19.bed -b ../data/hg19/encode/peaks/$i | wc -l |xargs echo obs >> ../results/tables/p63_exp_vs_ovs_overlap/hg19_p63_obs_vs_exp_$i
done;
