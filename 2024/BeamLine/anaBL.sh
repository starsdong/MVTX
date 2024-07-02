#!/bin/bash

for RUN in `cat run.list`;
do
  echo "Processing RUN # $RUN"
  root -l -q 'fitBL.C('${RUN}')'
  sleep 2
done
