#!/bin/bash


if [ $# -lt 1 ]; then
	echo $0 settingFile 
fi

source $1


normalizer=${samples[0]}
TAB=`echo -e "\t"`


rm *.00

echo "normalizing normalizer $normalizer"
Splidar.Expression.normalizeSample.sh $normalizer $normalizer $FoldCutOff


splitlines.py $normalizer.$normalizer/$normalizer.RPKM-$normalizer.RPKM.AC.txt 1 $normalizer.phead.00,$normalizer.pcontent.00
cut -d"$TAB" -f1-3 $normalizer.pcontent.00 | sort -k1,1 > $normalizer.content.00
cut -d"$TAB" -f1-3 $normalizer.phead.00 > $normalizer.head.00
cat $normalizer.head.00 $normalizer.content.00 > $normalizer.values.00


counter=0



for((normalizeeI=1;normalizeeI<${#samples[@]};normalizeeI++)); do ####
	normalizee=${samples[$normalizeeI]}
	echo "normalizing $normalizee on $normalizer"


	Splidar.Expression.normalizeSample.sh $normalizee $normalizer $FoldCutOff


	splitlines.py $normalizee.$normalizer/$normalizee.RPKM-$normalizer.RPKM.AC.txt 1 $normalizee.head.00,$normalizee.$normalizer.pcontent.00
	pwd
	cut -d"$TAB" -f1-5 $normalizee.$normalizer.pcontent.00 | sort -k1,1 > ordered.00
	
	cut -d"$TAB" -f1-3 ordered.00 > $normalizee.$normalizer.content.00
	cut -d"$TAB" -f4-5 ordered.00 > tomerge.content.00

	cut -d"$TAB" -f4-5 $normalizee.head.00 > tomerge.head.00
	cat tomerge.head.00 tomerge.content.00 > tomerge.00

	diff -s -q $normalizer.content.00 $normalizee.$normalizer.content.00

	counter=`expr $counter + 1`
	if [ $counter -lt 2 ]; then	
		paste $normalizer.values.00 tomerge.00 > merged.00
	else
		mv merged.00 prev.00
		paste prev.00 tomerge.00 > merged.00
	fi
done

mv merged.00 merged.expression.txt

getNamesDesc.py $genome.config merged.expression.txt 1 2 > merged.expression.anno.txt 2> merged.expression.anno.err









