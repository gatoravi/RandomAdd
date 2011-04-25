#!/bin/sh

for i in `ls *.txt`
do
  ./nw_order $i > $i
done
