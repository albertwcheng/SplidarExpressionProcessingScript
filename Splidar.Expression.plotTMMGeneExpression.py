#! which python

'''



Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

'''
##
#
#
#  TMM normalization and plotting gene expression by Audic-Claverie Statistics
#
#
#
#

##ZEROEXP=0.00000001

from pylab import *
from math import log
from math import fabs
from sys import argv
from sys import stderr
from sys import stdout
from sys import exit
from math import sqrt
from getopt import getopt
from math import ceil,fabs

from numpy import median,mean

from pvalue_module import *

from AudicClaverieStatInterface import *

from TMMNormalization import TMMnormalizationFactor

from albertcommon import *

def loadList(filename):
	print >> stderr,"loading list",filename
	L=[];	
	
	try:
		fin=open(filename);
		for line in fin:
			line=line.strip();
			if(line=="_end_"):
				break;		
			L.append(line);

		fin.close();	
	except IOError:
		print >> stderr, filename,"not openable";
	
	return L;


def drawGeneByName(GeneID,X,Y,geneName,display,color):
	try:
		if(len(geneName)<1):
			return;

		#if(geneName[len(geneName)-1]!=','):
		#	geneName+=",";

		if(display==""):
			display=geneName; #.split(",")[0];

		inx=GeneID.index(geneName);
		plot((X[inx],),(Y[inx],),color+'o');
		text(X[inx],Y[inx],display,{'color' : color});
		print >> stderr, geneName, "[",X[inx],Y[inx],"]";
	except ValueError:
		print >> stderr, geneName, " not found, not texted";


def fillDEFlag(DEFlags,GeneIDs,RFlag,GFlag,logNormRValues,logGValues,ACPvalue,ACFDR,RReads,GReads,FDRCutOff,foldCutOff,logb):
	logfoldCutOff=log(foldCutOff)/logb
	#print >> stderr,"fill"
	for gid,rflag,gflag,lognormrvalue,loggvalue,acpval,acfdr,rval,gval in zip(GeneIDs,RFlag,GFlag,logNormRValues,logGValues,ACPvalue,ACFDR,RReads,GReads):
		#print >> stderr,"filling DEFlag"
		if int(rflag)<1 and int(gflag)<1:
					#Not Expressed
					DEFlags.append("NotExpressed")
					#print >> stderr,"filling DEFlag NE"
					continue
			
		if acfdr<FDRCutOff:
			#differential!
			if rval>=gval:
				#r value is bigger!
				#so it's up
				if fabs(lognormrvalue-loggvalue)>=logfoldCutOff:  #fold and stat pass
					#Up2
					DEFlags.append("Up2")				
				else:
					#Up
					DEFlags.append("Up")

				#print >> stderr,"filling DEFlag Pan Up"
			else:
				if fabs(lognormrvalue-loggvalue)>=logfoldCutOff:  #fold and stat pass
					#Down2
					DEFlags.append("Down2")	

				else:
					#Down
					DEFlags.append("Down")		
				#print >> stderr,"filling DEFlag pan Down"	
		else:
			#Middle
			DEFlags.append("Middle")	
			#print >> stderr,"filling DEFlag Middle"

def fillDEFlagRPKM(DEFlags,GeneIDs,RFlag,GFlag,logNormRValues,logGValues,ACPvalue,ACFDR,RReads,GReads,FDRCutOff,foldCutOff,logb):
	logfoldCutOff=log(foldCutOff)/logb
	#print >> stderr,"fill"
	for gid,rflag,gflag,lognormrvalue,loggvalue,acpval,acfdr,rval,gval in zip(GeneIDs,RFlag,GFlag,logNormRValues,logGValues,ACPvalue,ACFDR,RReads,GReads):
		#print >> stderr,"filling DEFlag"
		if int(rflag)<1 and int(gflag)<1:
					#Not Expressed
					DEFlags.append("NotExpressed")
					#print >> stderr,"filling DEFlag NE"
					continue
			
		if acfdr<FDRCutOff:
			#differential!
			if rval>=gval:
				#r value is bigger!
				#so it's up
				if fabs(lognormrvalue-loggvalue)>=logfoldCutOff:  #fold and stat pass
					#Up2
					DEFlags.append("Up2")				
				else:
					#Up
					DEFlags.append("Up")

				#print >> stderr,"filling DEFlag Pan Up"
			else:
				if fabs(lognormrvalue-loggvalue)>=logfoldCutOff:  #fold and stat pass
					#Down2
					DEFlags.append("Down2")	

				else:
					#Down
					DEFlags.append("Down")		
				#print >> stderr,"filling DEFlag pan Down"	
		else:
			#Middle
			DEFlags.append("Middle")	
			#print >> stderr,"filling DEFlag Middle"

def selectVectorByMatchingFlag(listOfListSelectee,FlagVector,FlagValues,ignoreNAN=True):
	selectedListOfList=[]
	
	for LL in listOfListSelectee:
		selectedListOfList.append([])

	for i in range(0,len(FlagVector)):
		#print >> stderr,FlagVector[i],FlagValues
		
		if FlagVector[i] in FlagValues:
			#selected!
			#print >> stderr,"selected!"
			hasNAN=False

			if ignoreNAN:
				for j in range(0,len(listOfListSelectee)):
					if listOfListSelectee[j][i]=="nan":
						hasNAN=True
						break

				if hasNAN:
					continue
			for j in range(0,len(listOfListSelectee)):
				selectedListOfList[j].append(listOfListSelectee[j][i])
			
	return selectedListOfList
	



def writeOutDataAC(foutData,UpGeneIDAC,UpRFlagAC,UplogRAC,UpGFlagAC,UplogGAC,Flag,logData,UpPvalueAC,UpFDRAC,UpGCounts,UpRCounts):
	for (geneid,rflag,logr,gflag,logg,pval,fdr,ugc,urc) in zip(UpGeneIDAC,UpRFlagAC,UplogRAC,UpGFlagAC,UplogGAC,UpPvalueAC,UpFDRAC,UpGCounts,UpRCounts):
		
		m=logr-logg
		
		#foutData.write("GeneID\t"+GLabel+".Flag\tlog"+GLabel+"\t"+RLabel+".Flag\tlog"+RLabel+"\tA(avglog)\tM(logRatio)\tRatio\tFlag\tAudicClaveriePvalue\tAudicClaverieFDR\n");
		#foutData.write(geneid+"\t"+gflag+"\t"+str(logg)+"\t"+rflag+"\t"+str(logr)+"\t"+str(a)+"\t"+str(m)+"\t"+str(logData**m)+"\t"+Flag+"\t"+str(pval)+"\t"+str(fdr)+"\t"+str(gread)+"\t"+str(rread)+"\t"+str(totalGReads)+"\t"+str(totalRReads)+"\n");
		
		#foutData.write("GeneID\t"+GLabel+".Flag\tlog"+GLabel+"\t"+RLabel+".Flag\tlog"+RLabel+"\tlogRatioRPKM\tRatio\tFlag\tAudicClaveriePvalue\tAudicClaverieFDR\n");
		foutData.write(geneid+"\t"+gflag+"\t"+str(logg)+"\t"+rflag+"\t"+str(logr)+"\t"+str(m)+"\t"+str(logData**m)+"\t"+Flag+"\t"+str(pval)+"\t"+str(fdr)+"\t"+str(ugc)+"\t"+str(urc)+"\n")


def plotTMMGeneExpression(filename,houseKeepingGeneList,geneNamesToPlotFileName,plotName,startRow1,GeneIdCol,RCol,RFlagCol,GCol,GFlagCol,outFilePrefix,formatFigure,useSmartName,logData,RReadCol,GReadCol,totalRReads,totalGReads,FDRCutOff,foldCutOff,sizeDot,showNumber,alpha,beta,minValueAdd,sep,logratioTrim,sumTrim,doWeighting,Acutoff):
	

	lino=0
	
	logb=log(logData)
	
	RReads=[]
	GReads=[]

	RValues=[]
	GValues=[]

	MValues=[]
	AValues=[]
	ACPvalue=[]
	normRValues=[]
	normRReads=[]

	RFlag=[]
	GFlag=[]

	logNormRValues=[]
	logRValues=[]
	logGValues=[]

	GeneIDs=[]
	SNR=[]

	TMM_M=[]
	TMM_A=[]
	TMM_normM=[]
	TMM_normA=[]
	
	#load housekeeping genes
	houseKeepingGenes=loadList(houseKeepingGeneList)
	
	#load gene names to plot
	toPlotGenes=loadList(geneNamesToPlotFileName)

	#start reading data from file
	fin=open(filename)

	
		
	if(useSmartName):
		#use the first line
		line=fin.readline();
		line=line.strip();
		spliton=line.split(sep);
		RLabel=spliton[RCol];
		GLabel=spliton[GCol];
		
		if(plotName==""):
			plotName=RLabel+" vs "+GLabel;

		plotNameLessFancy=RLabel+"-"+GLabel;
		yaxisLabel=RLabel+"/"+GLabel;

		outFileFigurenormMA=outFilePrefix+plotNameLessFancy+".normMA."+formatFigure
		outFileFigureMA=outFilePrefix+plotNameLessFancy+".MA."+formatFigure     #use as prefix
		outFileFigureRGAC=outFilePrefix+plotNameLessFancy+".RGAC."+formatFigure
		outFileDataAC=outFilePrefix+plotNameLessFancy+".AC.txt"

		fin.close();
		fin=open(filename);
	



	
	for line in fin:
		lino+=1;
			
		line=line.strip();
		spliton=line.split(sep);


		if(lino<startRow1):
			
			continue; #not neccessary?						
		else:
			RValue=float(spliton[RCol]);
			GValue=float(spliton[GCol]);

			
			RRead=(float(spliton[RReadCol]))
			GRead=(float(spliton[GReadCol]))

			RReads.append(RRead)
			GReads.append(GRead)

			
			logRValue=log(RValue)/logb;
			logGValue=log(GValue)/logb;	

			GeneIDs.append(spliton[GeneIdCol].split(",")[0]);
			RValues.append(RValue);
			GValues.append(GValue);

			
			logRValues.append(logRValue)
			logGValues.append(logGValue)

			try:
				TMM_M.append((log(float(RRead)/totalRReads)-log(float(GRead)/totalGReads))/logb)
				TMM_A.append((log(float(RRead)/totalRReads)+log(float(GRead)/totalGReads))/logb/2.0)
			except:
				TMM_M.append("nan")
				TMM_A.append("nan")

			if(RFlagCol<0):
				RFlag.append('-');
			else:
				RFlag.append(spliton[RFlagCol]);

			if(GFlagCol<0):
				GFlag.append('-');
			else:			
				GFlag.append(spliton[GFlagCol]);
			


	fin.close();


	#filled in:
	#
	#RValues = RPKM of R (not normalized)
	#GValues = RPKM of G (ref)
	#RRead = R reads (not normalized)
	#GRead = G reads (ref)
	#logRValues = log RPKM of R (not normalized)
	#logGValues = log RPKM of G (ref)
	#GeneIDSs
	

	#to update
	#normTotalRReads

	#to fill in:
	#normRValues
	#lognormRValues
	#normRReads
	#TMM_normM
	#TMM_normA

	Fnorm=TMMnormalizationFactor(RReads,GReads,logratioTrim,sumTrim,doWeighting,Acutoff)
	print >> stderr,"TMM normalization factor=",Fnorm
	#get the normalization constant
	#first of all the effectiveSizes
	totalNormRReads=((totalRReads*Fnorm))
	#now the Rvalues

	for rvalue,rreads,greads in zip(RValues,RReads,GReads):
		norm_rvalue=rvalue*Fnorm
		norm_rreads=(rreads*Fnorm)		

		normRValues.append(norm_rvalue)
		normRReads.append(norm_rreads)
		logNormRValues.append(log(norm_rvalue)/logb)

		try:
			TMM_normM.append((log(float(norm_rreads)/totalRReads)-log(float(greads)/totalGReads))/logb) #specifically not using totalNormRReads
			TMM_normA.append((log(float(norm_rreads)/totalRReads)+log(float(greads)/totalGReads))/logb/2.0) #specifically not using totalNormRReads
		except:
			TMM_normM.append("nan")
			TMM_normA.append("nan")
	ACpvalue=[]
	ACFDR=[]
	

	

	nnnn=0
	NNNN=len(RReads)
	
		
	XYN1N2=[]

	#using normalized R reads and normalized R library size

	#RtoUse=normRReads
	#RTotalToUse=totalNormRReads
	RtoUse=RReads
	RTotalToUse=totalRReads

	for r,g in zip(RtoUse,GReads):
		nnnn+=1
			
		if r>=g: #the bigger one is the y
			y=r
			x=g
			n2=RTotalToUse
			n1=totalGReads
		else:
			y=g
			x=r
			n2=totalGReads
			n1=RTotalToUse

		XYN1N2.append([int(x),int(y),int(n1),int(n2)])
			
		
		if nnnn%1000==0:
			print >> stderr,"calculate Audic Claverie",nnnn,"of",NNNN
			if len(XYN1N2)>0:
				AudicClaverieStatInPlace(XYN1N2)

				for x,y,n1,n2,dummy,pvalue in XYN1N2:
					ACpvalue.append(pvalue)

			XYN1N2=[]

	if len(XYN1N2)>0:
		AudicClaverieStatInPlace(XYN1N2)
		for x,y,n1,n2,dummy,pvalue in XYN1N2:
			ACpvalue.append(pvalue)

	ACFDR=getFDRfromPvalue(ACpvalue)
		
		
		
	

	DEFlags=[]   #"Up2" "Down2" "Up" "Down" "Middle" "NotExpressed"

	#logNormRValues,logGValues are RPKM
	fillDEFlag(DEFlags,GeneIDs,RFlag,GFlag,logNormRValues,logGValues,ACpvalue,ACFDR,RtoUse,GReads,FDRCutOff,foldCutOff,logb)
	
	#print >> stderr,"DEFlaginputlengths=",len(GeneIDs),len(RFlag),len(GFlag),len(logNormRValues),len(logGValues),len(ACpvalue),len(ACFDR),len(RReads),len(GReads)


	#print >> stderr,"DEFlag len",len(DEFlags)

#outFileFigurenormMA=outFilePrefix+plotNameLessFancy+".normMA."+formatFigure
#outFileFigureMA=outFilePrefix+plotNameLessFancy+".MA."+formatFigure     #use as prefix
#outFileFigureRGAC=outFilePrefix+plotNameLessFancy+".RGAC."+formatFigure
#outFileDataAC=outFilePrefix+plotNameLessFancy+".AC.txt"

#MA plot before normalization
	print >> stderr,"plotting pre-normalized MA plot"
	figure(figsize=(8, 8));
	


	UpA,UpM=selectVectorByMatchingFlag([TMM_A,TMM_M],DEFlags,["Up","Up2"]) ##
	MiddleA,MiddleM=selectVectorByMatchingFlag([TMM_A,TMM_M],DEFlags,["Middle"]) ##
	DownA,DownM=selectVectorByMatchingFlag([TMM_A,TMM_M],DEFlags,["Down","Down2"]) ##
	#NEA,NEM=selectVectorByMatchingFlag([TMM_A,TMM_M],DEFlags,["NotExpressed"]) ##

	HKA,HKM=selectVectorByMatchingFlag([TMM_A,TMM_M],GeneIDs,houseKeepingGenes) ##

	if len(UpA)>0:
		scatter(UpA,UpM,marker='o',color='red',s=sizeDot);
	if len(MiddleA)>0:
		scatter(MiddleA,MiddleM,marker='o',color='blue',s=sizeDot);
	
	if len(DownA)>0:
		scatter(DownA,DownM,marker='o',color='green',s=sizeDot);

	#if len(NEA)>0:
	#	scatter(NEA,NEM,marker='o',color='grey',s=sizeDot);

	if len(HKA)>0:
		scatter(HKA,HKM,marker='o',color='grey',s=sizeDot)
		medianM=median(HKM)
		#meanM=mean(HKM)
		axhline(y=medianM, linewidth=1, color='cyan', linestyle="-")
		axhline(y=log(Fnorm)/logb, linewidth=1, color='yellow', linestyle="-.")
		#axhline(y=meanM, linewidth=1, color='gray', linestyle=":")

	axhline(linewidth=1, color='gray', linestyle="-")

	

	for geneName in toPlotGenes:
		geneNameS=geneName.split("=");
		geneID=geneNameS[0];
		
		lspliton=len(geneNameS);
		geneDisplay="";		
		color='k';		
		
		if(lspliton>1 and len(geneNameS[1])>0):
			geneDisplay=geneNameS[1];
		
		if(lspliton>2 and len(geneNameS[2])>0):
			color=geneNameS[2];
		
		
		drawGeneByName(GeneIDs,TMM_A,TMM_M,geneID,geneDisplay,color); ##

	

	title(plotName);
	xlabel("non normalized A" );
	ylabel("non normalized M");

	savefig(outFileFigureMA,format=formatFigure);

#MA plot after normalization
	print >> stderr,"plotting normalized MA plot"
	figure(figsize=(8, 8));
	
	UpA,UpM=selectVectorByMatchingFlag([TMM_normA,TMM_normM],DEFlags,["Up","Up2"]) ##
	MiddleA,MiddleM=selectVectorByMatchingFlag([TMM_normA,TMM_normM],DEFlags,["Middle"]) ##
	DownA,DownM=selectVectorByMatchingFlag([TMM_normA,TMM_normM],DEFlags,["Down","Down2"]) ##
	#NEA,NEM=selectVectorByMatchingFlag([TMM_normA,TMM_normM],DEFlags,["NotExpressed"]) ##

	HKA,HKM=selectVectorByMatchingFlag([TMM_normA,TMM_normM],GeneIDs,houseKeepingGenes) ##

	if len(UpA)>0:
		scatter(UpA,UpM,marker='o',color='red',s=sizeDot);

	if len(MiddleA)>0:
		scatter(MiddleA,MiddleM,marker='o',color='blue',s=sizeDot);
	
	if len(DownA)>0:
		scatter(DownA,DownM,marker='o',color='green',s=sizeDot);

	#if len(NEA)>0:
	#	scatter(NEA,NEM,marker='o',color='grey',s=sizeDot);
	
	if len(HKA)>0:
		scatter(HKA,HKM,marker='o',color='grey',s=sizeDot)
		medianM=median(HKM)
		#meanM=mean(HKM)
		axhline(y=medianM, linewidth=1, color='cyan', linestyle="-")
		#axhline(y=meanM, linewidth=1, color='gray', linestyle=":")

	axhline(linewidth=1, color='gray', linestyle="-")

	for geneName in toPlotGenes:
		geneNameS=geneName.split("=");
		geneID=geneNameS[0];
		
		lspliton=len(geneNameS);
		geneDisplay="";		
		color='k';		
	
		if(lspliton>1 and len(geneNameS[1])>0):
			geneDisplay=geneNameS[1];
		
		if(lspliton>2 and len(geneNameS[2])>0):
			color=geneNameS[2];
		
		
		drawGeneByName(GeneIDs,TMM_normA,TMM_normM,geneID,geneDisplay,color); ##



	title(plotName);
	xlabel("TMM-normalized A" );
	ylabel("TMM-normalized M");

	savefig(outFileFigurenormMA,format=formatFigure);


	middlecolor=[0.0,0.0,1.0]
	up2color=[1.0,0.0,0.0]
	down2color=[0.0,1.0,0.0]
	necolor="grey"

#RG plot with Audic Claverie after TMM-Normalized
	Up2PvalueAC,Up2FDRAC,Up2RFlagAC,Up2GFlagAC,Up2GeneIDAC,Up2logGAC,Up2logRAC,Up2RCounts,Up2GCounts=selectVectorByMatchingFlag([ACpvalue,ACFDR,RFlag,GFlag,GeneIDs,logGValues,logNormRValues,RtoUse,GReads],DEFlags,["Up2"]) ##	
	UpPvalueAC,UpFDRAC,UpRFlagAC,UpGFlagAC,UpGeneIDAC,UplogGAC,UplogRAC,UpRCounts,UpGCounts=selectVectorByMatchingFlag([ACpvalue,ACFDR,RFlag,GFlag,GeneIDs,logGValues,logNormRValues,RtoUse,GReads],DEFlags,["Up"]) ##
	MiddlePvalueAC,MiddleFDRAC,MiddleRFlagAC,MiddleGFlagAC,MiddleGeneIDAC,MiddlelogGAC,MiddlelogRAC,MiddleRCounts,MiddleGCounts=selectVectorByMatchingFlag([ACpvalue,ACFDR,RFlag,GFlag,GeneIDs,logGValues,logNormRValues,RtoUse,GReads],DEFlags,["Middle"]) ##
	DownPvalueAC,DownFDRAC,DownRFlagAC,DownGFlagAC,DownGeneIDAC,DownlogGAC,DownlogRAC,DownRCounts,DownGCounts=selectVectorByMatchingFlag([ACpvalue,ACFDR,RFlag,GFlag,GeneIDs,logGValues,logNormRValues,RtoUse,GReads],DEFlags,["Down"]) ##
	Down2PvalueAC,Down2FDRAC,Down2RFlagAC,Down2GFlagAC,Down2GeneIDAC,Down2logGAC,Down2logRAC,Down2RCounts,Down2GCounts=selectVectorByMatchingFlag([ACpvalue,ACFDR,RFlag,GFlag,GeneIDs,logGValues,logNormRValues,RtoUse,GReads],DEFlags,["Down2"]) ##
	NEPvalueAC,NEFDRAC,NERFlagAC,NEGFlagAC,NEGeneIDAC,NElogGAC,NElogRAC,NERCounts,NEGCounts=selectVectorByMatchingFlag([ACpvalue,ACFDR,RFlag,GFlag,GeneIDs,logGValues,logNormRValues,RtoUse,GReads],DEFlags,["NotExpressed"]) ##

	for x in [alpha]:
		for y in [beta]:
#RG plot with Audic Claverie
			print >> stderr,"working on AC x=",x,"y=",y
			X=1-x
			Y=1-y
			mr,mg,mb=middlecolor
			ur,ug,ub=up2color
			dr,dg,db=down2color
			upcolor=[mr*x+ur*X,mg*x+ug*X,mb*x+ub*X]
			downcolor=[mr*y+dr*Y,mg*y+dg*Y,mb*y+db*Y]
			
			figure(figsize=(8, 8));
			#scatter(SNR,MValues,marker='o',color='black',s=sizeDot);
			
			
			

			
			if len(MiddlelogGAC)>0:
				scatter(MiddlelogGAC,MiddlelogRAC,marker='o',color=middlecolor,s=sizeDot);

			if len(NElogGAC)>0:
				scatter(NElogGAC,NElogRAC,marker='o',color=necolor,s=sizeDot);

			if len(UplogGAC)>0:
				scatter(UplogGAC,UplogRAC,marker='o',color=upcolor,s=sizeDot);

	

			if len(DownlogGAC)>0:
				scatter(DownlogGAC,DownlogRAC,marker='o',color=downcolor,s=sizeDot);


			if len(Up2logGAC)>0:
				scatter(Up2logGAC,Up2logRAC,marker='o',color=up2color,s=sizeDot);


			if len(Down2logGAC)>0:
				scatter(Down2logGAC,Down2logRAC,marker='o',color=down2color,s=sizeDot);

			title(plotName);
			xlabel("log"+GLabel);
			ylabel("log"+RLabel);

			for geneName in toPlotGenes:
				geneNameS=geneName.split("=");
				geneID=geneNameS[0];
		
				lspliton=len(geneNameS);
				geneDisplay="";		
				color='k';		
	
				if(lspliton>1 and len(geneNameS[1])>0):
					geneDisplay=geneNameS[1];
		
				if(lspliton>2 and len(geneNameS[2])>0):
					color=geneNameS[2];
		
		
				drawGeneByName(GeneIDs,logGValues,logNormRValues,geneID,geneDisplay,color);
	
	
			#if showNumber:
			#	lUplogGAC=len(UplogGAC)
			#	lMiddlelogGAC=len(MiddlelogGAC)
			#	lDownlogGAC=len(DownlogGAC)

			#	text(min(logGValues),max(logNormRValues),str(lUplogGAC),color='red',horizontalalignment='left',verticalalignment='top')
			#	text(max(logGValues),max(logNormRValues),str(lMiddlelogGAC),color='black',horizontalalignment='right',verticalalignment='top')
			#	text(max(logGValues),min(logNormRValues),str(lDownlogGAC),color='green',horizontalalignment='right',verticalalignment='bottom')

			xlim(xmin=min(logGValues)-minValueAdd,xmax=max(logGValues)*1.1)
			ylim(ymin=min(logNormRValues)-minValueAdd,ymax=max(logNormRValues)*1.1)

			savefig(outFileFigureRGAC,format=formatFigure);
			#savefig(outFileFigureRGAC+"_"+str(x)+"_"+str(y)+".png",format=formatFigure);

#OutputFiles
	foutData=open(outFileDataAC,"w")
	
	foutData.write("GeneID\t"+GLabel+".Flag\tlog"+GLabel+"\t"+RLabel+".Flag\tlog"+RLabel+"\tlogRatioRPKM\tRatio\tFlag\tAudicClaveriePvalue\tAudicClaverieFDR\tGCount\tRCount\n")
	#normalized
	writeOutDataAC(foutData,Up2GeneIDAC,Up2RFlagAC,Up2logRAC,Up2GFlagAC,Up2logGAC,"Up2",logData,Up2PvalueAC,Up2FDRAC,Up2GCounts,Up2RCounts)	
	writeOutDataAC(foutData,UpGeneIDAC,UpRFlagAC,UplogRAC,UpGFlagAC,UplogGAC,"Up",logData,UpPvalueAC,UpFDRAC,UpGCounts,UpRCounts)	
	writeOutDataAC(foutData,Down2GeneIDAC,Down2RFlagAC,Down2logRAC,Down2GFlagAC,Down2logGAC,"Down2",logData,Down2PvalueAC,Down2FDRAC,Down2GCounts,Down2RCounts)	
	writeOutDataAC(foutData,DownGeneIDAC,DownRFlagAC,DownlogRAC,DownGFlagAC,DownlogGAC,"Down",logData,DownPvalueAC,DownFDRAC,DownGCounts,DownRCounts)	
	writeOutDataAC(foutData,MiddleGeneIDAC,MiddleRFlagAC,MiddlelogRAC,MiddleGFlagAC,MiddlelogGAC,"Middle",logData,MiddlePvalueAC,MiddleFDRAC,MiddleGCounts,MiddleRCounts)	
	writeOutDataAC(foutData,NEGeneIDAC,NERFlagAC,NElogRAC,NEGFlagAC,NElogGAC,"NotExpressed",logData,NEPvalueAC,NEFDRAC,NEGCounts,NERCounts)
	
	foutData.close()


def printUsageAndExit(programName):
	print >> stderr,programName,"[options] filename GeneIdCol RCol[RPKM] RFlagCol RReadCol totalRReads GCol[RPKM] GFlagCol GReadCol totalGReads outFilePrefix"
	print >> stderr,"[options]"
	print >> stderr,"--Acutoff x [-1e10]  minimal average log2 value for TMM trimming"
	print >> stderr,"--not-doWeighting    turn off weighting by inverse of variance in TMM normalization"
	print >> stderr,"--sumTrim p [0.05]  double trim the lowest and highest p percentile of sum (A) in TMM normalization"
	print >> stderr,"--logratioTrim p [0.3] double trim the lowest and highest p percentile of log ratio (M) in TMM normalization"
	print >> stderr,"--sep c  [TAB]  set separator symbol"
	print >> stderr,"--minValueAdd x [0.02] add a small value for xlim and ylim low bound in order to show a whole dot for the lowest points"
	print >> stderr,"--alpha x [0.5] specify the alpha for the RG AC plot"
	print >> stderr,"--beta y [0.5] specify the beta for the RG AC plot"
	print >> stderr,"--showNumber   turn on the show number feature"
	print >> stderr,"--sizeDot x [5] set the size of dots in graphs"
	print >> stderr,"--foldCutOff f [2.0] set the fold cutoff"
	print >> stderr,"--FDRCutOff F [0.05] set the FDR cutoff for AC statistics"
	print >> stderr,"--geneNamesToPlotFile filename [""] plot gene listsed in the file"
	print >> stderr,"--plotName name [""] explicitly specify the plot name"
	print >> stderr,"--startRow r [2] data starts from row r"
	print >> stderr,"--formatFigure ext [png] the format for figure output"
	print >> stderr,"--not-useSmartName turn off the use smart name feature"
	print >> stderr,"--logData logbase [10] log the data by logbase"
	print >> stderr,"--houseKeppingGenes filename [""] load the housekeeping gene list for diagnostic purposes"
	explainColumns(stderr)
	exit()

if __name__=="__main__":
	alpha=0.5
	beta=0.5
	minValueAdd=0.02
	showNumber=False
	logratioTrim=0.3
	sumTrim=0.05
	doWeighting=True
	Acutoff=-1e10
	startRow=2
	sep="\t"
	foldCutOff=2.0

	sizeDot=5
	FDRCutOff=0.05

	geneNamesToPlotFile=""
	plotName=""
	formatFigure="png"
	useSmartName=True

	logData=10
	#geneNamesToPlotFileName=""
	houseKeepingGeneList=""

	programName=argv[0]
	opts,args=getopt(argv[1:],'',['Acutoff=','not-doWeighting','sumTrim=','logratioTrim=','sep=','minValueAdd=','alpha=','beta=','showNumber','sizeDot=','foldCutOff=','FDRCutOff=','geneNamesToPlotFile=','plotName=','startRow=','formatFigure=','not-not-useSmartName','logData=','houseKeepingGenes='])
	
	try:
		filename,GeneIdCol,RCol,RFlagCol,RReadCol,totalRReads,GCol,GFlagCol,GReadCol,totalGReads,outFilePrefix=args
	except:
		printUsageAndExit(programName)
	
	for o,v in opts:
		if o=="--Acutoff":
			Acutoff=float(v)
		elif o=="--not-doWeighting":
			doWeighting=False
		elif o=="--sumTrim":
			sumTrim=float(v)
		elif o=="--logratioTrim":
			logratioTrim=float(v)
		elif o=="--sep":
			sep=replaceSpecialChar(v)
		elif o=="--minValudAdd":
			minValudAdd=float(v)
		elif o=="--alpha":	
			alpha=float(v)
		elif o=="--beta":
			beta=float(v)
		elif o=="--showNumber":
			showNumber=True
		elif o=="--sizeDot":
			sizeDot=float(v)
		elif o=="--foldCutOff":
			foldCutOff=float(v)
		elif o=="--FDRCutOff":
			FDRCutOff=float(v)
		elif o=="--geneNamesToPlotFile":
			geneNamesToPlotFile=v
		elif o=="--plotName":
			plotName=v
		elif o=="--startRow":
			startRow=int(v)
		elif o=="--formatFigure":
			formatFigure=v
		elif o=="--not-useSmartName":
			useSmartName=False
		elif o=="--logData":
			logData=float(v)
		elif o=="--houseKeepingGenes":
			houseKeepingGeneList=v



	headerRow=startRow-1
	
	header,prestarts=getHeader(filename,headerRow,startRow,sep)
		
	GeneIdCol=getCol0ListFromCol1ListStringAdv(header,GeneIdCol)[0]
	RCol=getCol0ListFromCol1ListStringAdv(header,RCol)[0]
	RFlagCol=getCol0ListFromCol1ListStringAdv(header,RFlagCol)[0]	
	RReadCol=getCol0ListFromCol1ListStringAdv(header,RReadCol)[0]	
	GCol=getCol0ListFromCol1ListStringAdv(header,GCol)[0]
	GFlagCol=getCol0ListFromCol1ListStringAdv(header,GFlagCol)[0]	
	GReadCol=getCol0ListFromCol1ListStringAdv(header,GReadCol)[0]
	
	totalRReads=int(totalRReads)
	totalGReads=int(totalGReads)

	plotTMMGeneExpression(filename,houseKeepingGeneList,geneNamesToPlotFile,plotName,startRow,GeneIdCol,RCol,RFlagCol,GCol,GFlagCol,outFilePrefix,formatFigure,useSmartName,logData,RReadCol,GReadCol,totalRReads,totalGReads,FDRCutOff,foldCutOff,sizeDot,showNumber,alpha,beta,minValueAdd,sep,logratioTrim,sumTrim,doWeighting,Acutoff)



