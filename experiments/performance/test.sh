#!/bin/bash

# directory of this script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

DATASET=$DIR/../../data/SILVA_119_SSURef_tax_silva.fasta
KLUST=$DIR/../../src/klust

for k in `seq 4 8`; do
  for t in `seq 0.75 0.05 0.95`; do
    $KLUST $DATASET -k $k -t $t -c 100000 2>&1 | tee "klust_${k}_${t}"

    $KLUST $DATASET -k $k -t $t -c 100000 -i 2>&1 | tee "klust_${k}_${t}_sort_incr"

    $KLUST $DATASET -k $k -t $t -c 100000 -d 2>&1 | tee "klust_${k}_${t}_sort_decr"
  done
done

for k in `seq 5 6`; do
  for t in `seq 0.8 0.05 0.9`; do
    $KLUST $DATASET -o "SILVA_FULL_${k}_${t}_cts" -u "SILVA_FULL_${k}_${t}_cls" \
    	 -k $k -t $t  2>&1 | tee "SILVA_FULL_klust_${k}_${t}"

    $KLUST $DATASET -o "SILVA_FULL_${k}_${t}_cts_incr" -u "SILVA_FULL_${k}_${t}_cls_incr" \
    	 -k $k -t $t -i 2>&1 | tee "SILVA_FULL_klust_${k}_${t}_sort_incr"
  done
done