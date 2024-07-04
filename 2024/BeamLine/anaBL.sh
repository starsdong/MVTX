#!/bin/bash

for RUN in `cat $1.list`;
do
  echo "Processing RUN # $RUN"
  root -l -q 'fitBL.C('${RUN}')'
  sleep 2
done
