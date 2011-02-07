#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 normalizee normalizer foldCutOff
	exit
fi

normalizee=$1
normalizer=$2
foldCutOff=$3

normalizeeFile=$normalizee/$normalizee.expression.txt
normalizerFile=$normalizer/$normalizer.expression.txt

source $normalizee/$normalizee.statEI
totalRReads=$totalReadsMapped
echo $normalizee has $totalRReads total reads

source $normalizer/$normalizer.statEI
totalGReads=$totalReadsMapped
echo $normalizer has $totalGReads total reads


TAB=`echo -e "\t"`

nn=$normalizee.$normalizer

mkdir $nn

cut -d"$TAB" -f1-3 $normalizerFile > $nn/$normalizer.info.00
cut -d"$TAB" -f1-3 $normalizeeFile > $nn/$normalizee.info.00

diff -s -q $nn/$normalizer.info.00 $nn/$normalizee.info.00

cut -d"$TAB" -f4-16  $normalizerFile > $nn/$normalizer.00
paste $normalizeeFile $nn/$normalizer.00 > $nn/$normalizee.$normalizer.00

source $normalizee/$normalizee.statEI
totalRReads=$totalReadsMapped
echo $normalizee has $totalRReads total reads
totalURReads=`expr $totalExonicReads + $totalNonExonicReads`
echo $normalizee has $totalURReads total unique reads

source $normalizer/$normalizer.statEI
totalGReads=$totalReadsMapped
echo $normalizer has $totalGReads total reads
totalUGReads=`expr $totalExonicReads + $totalNonExonicReads`
echo $normalizee has $totalUGReads total unique reads

nn=$normalizee.$normalizer



#/net/coldfact/data/awcheng/scripts/plotTMMGeneExpression.py [options] filename GeneIdCol RCol[RPKM] RFlagCol RReadCol totalRReads GCol[RPKM] GFlagCol GReadCol totalGReads #outFilePrefix
#[options]
#--Acutoff x [-1e10]  minimal average log2 value for TMM trimming
#--not-doWeighting    turn off weighting by inverse of variance in TMM normalization
#--sumTrim p [0.05]  double trim the lowest and highest p percentile of sum (A) in TMM normalization
#--logratioTrim p [0.3] double trim the lowest and highest p percentile of log ratio (M) in TMM normalization
#--sep c  [TAB]  set separator symbol
#--minValueAdd x [0.02] add a small value for xlim and ylim low bound in order to show a whole dot for the lowest points
#--alpha x [0.5] specify the alpha for the RG AC plot
#--beta y [0.5] specify the beta for the RG AC plot
#--showNumber   turn on the show number feature
#--sizeDot x [5] set the size of dots in graphs
#--foldCutOff f [2.0] set the fold cutoff
#--FDRCutOff F [0.05] set the FDR cutoff for AC statistics
#--geneNamesToPlotFile filename [] plot gene listsed in the file
#--plotName name [] explicitly specify the plot name
#--startRow r [2] data starts from row r
#--formatFigure ext [png] the format for figure output
#--not-useSmartName turn off the use smart name feature
#--logData logbase [10] log the data by logbase
#cols format:
#        * inserts all columns
#        number: 1-5,3
#        number preceded by '_' means from the end. _1 means last row, _2 second last
#        field name preceded by '.': .exp-5,.pvalue-.fdr
#        regex preceded by '@': @FDR[a-z]
#        %filename%sep or just %filename: open file which contains the fields to include. default separator is [TAB]
#        if last field is a. Then following columns are added as ascending. If appears at end, all are sorted ascending
#        is last field is d. Then following columns are added as descending. If appears at end, all are sorted descending

#housekeeping genes from http://www.cgen.com/supp_info/Housekeeping_genes.html
#getNamesDesc.py --onlyReplace Sym hg18.config human_housekeepinggenes_nm.txt 1 1 | sort | uniq > human_housekeepinggenes_sym


#readsToUse=Reads+BG
readsToUse=UniqueReads


Splidar.Expression.plotTMMGeneExpression.py  --foldCutOff $foldCutOff --plotName "Gene Expression Profile M/E" $nn/$normalizee.$normalizer.00 .GeneName .$normalizee.RPKM .$normalizee.expressedFlag .$normalizee.$readsToUse $totalRReads .$normalizer.RPKM .$normalizer.expressedFlag .$normalizer.$readsToUse $totalGReads $nn/

#plotTMMGeneExpression.py  --foldCutOff $foldCutOff --houseKeepingGenes human_housekeepinggenes_sym  --geneNamesToPlotFile $scriptDir/toPlotGenes.txt --plotName "Gene Expression Profile M/E" $nn/$normalizee.$normalizer.00 .GeneName .$normalizee.RPKM .$normalizee.expressedFlag .$normalizee.$readsToUse $totalRReads .$normalizer.RPKM .$normalizer.expressedFlag .$normalizer.$readsToUse $totalGReads $nn/toPlot_

cd $nn

gcountCol=`colSelect.py "$normalizee.RPKM-$normalizer.RPKM.AC.txt" .GCount`
rcountCol=`colSelect.py "$normalizee.RPKM-$normalizer.RPKM.AC.txt" .RCount`

awk -v gc=$gcountCol -v rc=$rcountCol -v FS="\t" -v OFS="\t" '{if(FNR==1){n=NF; $(n+1)="B1"; $(n+2)="B2"; $(n+3)="B3";}else{$(n+1)=$gc*$rc; $(n+2)=$gc+$rc; if($gc>$rc){$(n+3)=$gc;}else{$(n+3)=$rc;} }print;}'  "$normalizee.RPKM-$normalizer.RPKM.AC.txt" > $normalizee.RPKM-$normalizer.RPKM.AC.wBStat.txt

cd ..



