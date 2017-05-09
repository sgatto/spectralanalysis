#! /bin/bash

for comp in {0..4..1};
  for time in {0..400..50};
    do
      python analysis.py -tstart $time -cstart $comp
    done
  done
