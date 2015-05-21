#!/bin/bash

# directory of this script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

DATASET=$DIR/../../data/SILVA_119_SSURef_tax_silva.fasta

KLUST=$DIR/../../src/klust

for k in `seq 4 8`; do
  for t in `seq 0.75 0.05 0.95`; do
    $KLUST $DATASET "cts_${k}_${t}" "cls_${k}_${t}" -k $k -t $t \
      2>&1 | tee "klust_${k}_${t}"
  done
done
