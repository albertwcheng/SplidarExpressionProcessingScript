#!/bin/bash

if [ $# -lt 1 ]; then
	echo `basename $0` settingFile #FDRCutOff(0.05) FoldCutOff(2)
	exit
fi

settingFile=$1

#add myself to path
source addMyPath.sh
myPath=`absdirname.py $0`
export PATH=${PATH}:${myPath}/AudicClaverieStat
export PYTHONPATH=${PYTHONPATH}:${myPath}

#combine count files from Splidar into one file per sample and then use Poisson background poisson model to process data
Splidar.Expression.mergeAllSample.sh $settingFile

#now normalize using TMM procedures and make a normalized file for all data using the first sample (samples[0] ) as the reference for normalization.
Splidar.Expression.normalizeAllSamples.sh $settingFile
