#!/bin/sh

# directory of this script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#DATASET=$DIR/../../data/p3_clean_A-14-Ileum_S27.fna
#DATASET_SORTED=$DIR/../../data/p3_clean_A-14-Ileum_S27_sorted.fna
DATASET=$DIR/../../data/RDP_10000.fna
DATASET_SORTED=$DIR/../../data/RDP_10000_sort-decr.fna

NUM=`seq 0.0 0.01 1.0`

# echo 0.1, ..., 1.0
for i in $NUM; do echo -n $i ' '; done

echo

# echo number of klusters for klust
for i in $NUM; do
  KLUST_COUNT=$(../../src/klust $DATASET klust_cts klust_cls \
    -t $i | grep '# of clusters' | sed 's/[^0-9]*//g')
  echo -n $KLUST_COUNT ' '
done

echo

# echo number of klusters for usearch
for i in $NUM; do
  # USEARCH prints copyright/license info to stdout and results to stderr (!?)
  USEARCH_COUNT=$(usearch -cluster_smallmem $DATASET_SORTED \
    -id $i 2>&1| grep 'Clusters ' | sed 's/[^0-9]*//g')
  echo -n $USEARCH_COUNT ' '
done

echo
