#!/bin/bash



if [ $# -lt 2 ]; then
	echo "usage:" $0 sampleName FDRCutOffForExpression
	exit
fi

sampleName=$1
FDRCutOffForExpression=$2
folder=$sampleName

echo "merging counts in $sampleName"

mergedName=$folder/$sampleName.merged.count

if [ -e $mergedName ]; then
	echo "$mergedName existed. remove"
	rm $mergedName
fi

cat $folder/$sampleName.*.count > $mergedName
echo "done merging $sampleName"


mergedNetName=$mergedName.net
slimName=$mergedNetName.slim
#now remove all NGB (No GBlock) rows
awk '$4!~/N/ && $5>0' $mergedName > $mergedNetName

#get main read info
echo -e "Chr\tStrand\tGeneName\t$sampleName.UniqueReads\t$sampleName.UniquPos" > $slimName
cut -f1-5 $mergedNetName >> $slimName

#filename,startRow1,prefixHeader,ReadCol1,PosCol1,totalReadsMapped,statEI
Splidar.Expression.calculateExpressionPvalueFDR.py $slimName 2 "$sampleName." 4 5 $folder/$sampleName.statEI $FDRCutOffForExpression > $folder/$sampleName.expression.txt



