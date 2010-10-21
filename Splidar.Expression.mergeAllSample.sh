#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 settingFile #FDRCutOffForExpression
fi

source $1

for i in ${samples[@]} ; do ###
	echo "trying to merge sample $i"
	Splidar.Expression.mergeSample.sh $i $FDRCutOffForExpression
done
