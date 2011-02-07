#!/usr/bin/env python
######
##version 2.0
##modified 6/9/2009
##albert cheng. the poisson noise model for solexa-based RPKM gene expression estimation
##awcheng@mit.edu
########

import warnings;

warnings.simplefilter("ignore",DeprecationWarning);
from sys import stdout;
from sys import stderr;
from poisson import poisson_cdf;
from sys import argv;
from configobj import ConfigObj

##ExonicUniqueLength is actually not required. because it is cancelled out in calculation of lambda

def calculateExpressionPvalueFDR(filename,startRow1,prefixHeader,ReadCol1,PosCol1,statEI,FDRCutOff,useUniqueReadsAsTotal):
	ReadCol0=ReadCol1-1;
	PosCol0=PosCol1-1;
	
	if useUniqueReadsAsTotal:
		totalReadsMapped=int(statEI["totalExonicReads"])+int(statEI["totalNonExonicReads"])
	else:
		totalReadsMapped=int(statEI["totalReadsMapped"])
	NonExonicUniqueReads=int(statEI["totalNonExonicReads"])
	NonExonicUniqueLength=int(statEI["totalNonExonicPos"])
	print >> stderr,"statEI"
	print >> stderr, statEI
	lineOrderedRecord=[]

	fin=open(filename);
	lino=0;
	pNonExonic=float(NonExonicUniqueReads)/NonExonicUniqueLength;   ##*ExonicUniqueLength

	PvaluesLinesMap=dict();
	
	scalingRPKM=float(1000)/(float(totalReadsMapped)/1000000); ##the so-called a Kb Exon model per million reads mapped~~

	total=0;

	minlogRPKM=100000

	for line in fin:
		line=line.strip();
		lino+=1;
		if lino<startRow1:
			print >> stdout, line+"\t"+prefixHeader+"expressedFlag"+"\t"+prefixHeader+"Reads+BG"+"\t"+prefixHeader+"DPB"+"\t"+prefixHeader+"RPKM"+"\t"+prefixHeader+"lambda_Counts"+"\t"+prefixHeader+"lamdba_DPB"+"\t"+prefixHeader+"lamdba_RPKM"+"\t"+prefixHeader+"var_DPB"+"\t"+prefixHeader+"var_RPKM"+"\t"+prefixHeader+"pvalue"+"\t"+prefixHeader+"FDR";
			continue;
		
		spliton=line.split("\t");
		readCount=int(spliton[ReadCol0]);
		pos=int(spliton[PosCol0]);
		
		#pg=float(pos)/ExonicUniqueLength;
		
		lam=pNonExonic*float(pos) #nr*pg;

		#print >> stderr, "nr,pg,lam,p-value=",nr,pg,lam,
		if(readCount<=1):
			pvalue=1;
		else:		
			pvalue=1-poisson_cdf(readCount-1,lam);
		
		#print >> stderr, pvalue;		
		
		#if readCount==0:
		#	readCountpBG=1 #readCount+lam; #poisson lam=mean=sd
		#else:
		
		readCountpBG=readCount+lam; #poisson lam=mean=sd
		
		flag="";		
		if(readCount==0):
			flag=-1;		
		#elif(readCount<lam):
		#	flag=0;
		else:
			flag=1;


		dpb=float(readCountpBG)/pos;
		RPKM=dpb*scalingRPKM;			
		
		lamDPB=float(lam)/pos;
		lamRPKM=lamDPB*scalingRPKM;

		varDPB=float(lam)/(pos**2);
		varRPKM=varDPB*(scalingRPKM**2);

		#print >> stdout, line+"\t"+str(lam)+"\t"+str(pvalue);	
		
		#a new record:
		thisRecord=[line,flag,str(readCountpBG),str(dpb),str(RPKM),str(lam),str(lamDPB),str(lamRPKM),str(varDPB),str(varRPKM),str(pvalue)]
		
		lineOrderedRecord.append(thisRecord)

		if not PvaluesLinesMap.has_key(pvalue):
			PvaluesLinesMap[pvalue]=[];
		PvaluesLinesMap[pvalue].append(thisRecord);
		total+=1;

	pvalueSorted=sorted(PvaluesLinesMap.keys());

	included=0;

	for pvalue in pvalueSorted:
		PvalueLine=PvaluesLinesMap[pvalue];
		thisSize=len(PvalueLine);
		included+=thisSize;
		FDR=pvalue*total/included;
		for line in PvalueLine:
			line.append(FDR)
			#print >> stdout, line +"\t"+str(FDR);

	fin.close();


	#now output in order of input:
	for line in lineOrderedRecord:
		FDR=line[-1]
		flag=line[1]
		if FDRCutOff<1 and FDR>=FDRCutOff and flag==1:
			flag=0
		line[1]=str(flag)
		line[-1]=str(FDR)
		print >> stdout,"\t".join(line)



if(len(argv)<9):
	print >> stderr, "Argument vector:",argv;
	print >> stderr, "Usage:",argv[0],"filename,startRow1,prefixHeader,ReadCol1,PosCol1,configFile,FDRCutOff[=2:noCutOff],useUniquelyReadsAsTotal[=y|n]";	
	exit();

##read in stat from config file!
configFile=argv[6]
statEI=ConfigObj(configFile);

FDRCutOff=float(argv[7])
useUniqueReadsAsTotal=(argv[8]=='y')
calculateExpressionPvalueFDR(argv[1],int(argv[2]),argv[3],int(argv[4]),int(argv[5]),statEI,FDRCutOff,useUniqueReadsAsTotal);





