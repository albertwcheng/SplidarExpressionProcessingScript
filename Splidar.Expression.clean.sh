#!/bin/bash


rm */*.expression.txt
rm */*.merged.*

rm *.00

for i in `ls -l | grep ^d | awk '($NF~/\./){printf("%s\n",$NF);}'`; do
rm -R $i
done