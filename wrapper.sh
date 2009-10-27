#!/bin/bash
region=('Region1' 'Region2' 'Region3')

cd [subjectdirectories]

for K in [subjectlist]
  do
  j=0
  for i in [ROIlist]
  do
    [Toolboxpath]/resting_connectivity.sh -S [subjectdirectories] -s $K -V ${i} -R ${region[$j]} --files [fileidentifier]
    j=$[${j}+1]
  done
done
