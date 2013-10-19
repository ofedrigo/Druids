#! /usr/bin/python

from Bio.Nexus import Nexus
from Bio import SeqIO
import math,random,sets,getopt,sys,os

# not used yet
def getFileData(fileData):
	symbol="\t"
	thresh={}
	if fileData[1].find(symbol)==-1:
		symbol=" "
	legends=fileData[1].split(symbol)[1:]
	threshValues=fileData[0].split(symbol)[1:]
	for i in range(len(legends)):
		if len(threshValues)>0:
			thresh[legends[i]]=threshValues[i]
		else:
			thresh[legends[i]]=""
	vector={}
	print legends
	for line in fileData[2:]:
		items=line.split(symbol)
		for i in range(len(legends)):
			vector.setdefault(legends[i],[]).append(float(items[i+1]))
	return vector,legends,thresh

# codon bootstrap
def bootstrap(matrix,nreps,windowSize,percent,function,nSites):
	listOfScores=[]
	for rep in range(nreps):
		listOfSites=[]
		for i in range(int(windowSize/3)): #watch out! this is only for /3 multiple
			site=random.randint(0,nSites-1)
			if getPosition(site)==1:
				listOfSites.append(site)
				listOfSites.append(site+1)
				listOfSites.append(site+2)
			if getPosition(site)==2:
				listOfSites.append(site-1)
				listOfSites.append(site)
				listOfSites.append(site+1)
			if getPosition(site)==3:
				listOfSites.append(site-2)
				listOfSites.append(site-1)
				listOfSites.append(site)
		listOfScores.append(function(matrix,listOfSites))
	listOfScores.sort()
	thresh=listOfScores[int(len(listOfScores)*(1.0-percent))]
	return thresh

# calculate mean
def getMean(vector):
	if len(vector)>0: return float(sum(vector))/float(len(vector))
	return 0.0

# output results in a list/site/attribute: a site is assigned the maximum F-test value (DFS) of every windows it belongs to
def saveList(outputFileName,vector,thresh):
	fileHandle = open(outputFileName, 'w')
	nsites=len(vector[vector.keys()[0]])
	fileHandle.write("thresh")
	for name in vector.keys():
		fileHandle.write("\t"+str(thresh[name]))
	fileHandle.write("\nsites\t"+"\t".join(vector.keys()))
	for site in range(nsites):
		fileHandle.write("\n"+str(site))
		for name in vector.keys():
			fileHandle.write("\t"+str(vector[name][site]))
	fileHandle.close()

# output significant windows
def saveWindows(outputFileName,vector):
	fileHandle = open(outputFileName, 'w')
	for name in vector.keys():
		fileHandle.write(name+"\t"+vector[name]+"\n")
	fileHandle.close()

# convert regions coordinates to a string	
def convertToString(windows):
	newLargeList=[]
	for window in windows:
		newList=[str(x) for x in window]
		newLargeList.append(newList)
	return ",".join(["-".join(v) for v in newLargeList])

# create regions with all the significant sites
def makeWindow(vector):
	results=[]
	vector.sort()
	start=vector[0]
	end=vector[0]
	for site in vector[1:len(vector)]:
		if site==end+1:
			end=site
		else:
			if start==end:
				results.append([start])
			else:
				results.append([start,end])
				start=site
				end=site
	return results

# get the list of sites that are above a threshhold
def getSig(vector,windowSize,thresh):
	results=sets.Set([])
	for site in range(len(vector)):
		if vector[site]>=thresh:
				for i in range(windowSize):
					results.add(site+i)
	return list(results)

#get codon position, assume your data is in frame and /3 multiple length
def getPosition(site):
	if (site+1+1)%3==0: return 2
	if (site+1+2)%3==0: return 1
	if (site+1)%3==0: return 3

# not used yet
def baseComp(sequences,nsites,speciesList,windowSize):
	allChi=[]
	for site in range(nsites): allChi.append(0.0)
	for site in range(nsites-windowSize):
		freqs={"A":0,"C":0,"G":0,"T":0}
		alphabet=freqs.keys()
		allFreqs={}
		Tsize={}
		T=0.0
		for species in speciesList:
			Sfreqs={"A":0,"C":0,"G":0,"T":0}
			cc=0
			for i in range(windowSize):
				base=sequences.matrix[species].tostring()[site+i]
				if base in alphabet:
					Sfreqs[base]=Sfreqs[base]+1
					freqs[base]=freqs[base]+1
					cc=cc+1
			Tsize[species]=cc
			T=T+float(cc)
			allFreqs[species]=Sfreqs
			
		chi_Square=0.0
		for base in alphabet:
			freqs[base]=(float(freqs[base])/float(len(speciesList)))/float(T)
			for species in speciesList:
				expected=Tsize[species]*freqs[base]
				observed=allFreqs[species][base]
				if expected>0:
					chi_Square=chi_Square+math.pow(observed-expected,2)/expected
		allChi[site]=chi_Square
	return allChi

# assgin attribute value for each site
def getProfiles(alphabet,seq,properties,genetic):
	profile={"GC":{},"GT":{},"GA":{},"AC":{},"AT":{},"TC":{},"H":{},"V":{},"C":{}}
	newSeqs=[]
	nchar=len(seq[seq.keys()[0]])
	for site in range(nchar):
		gc=0
		for species in seq.keys():
			base=seq[species][site]
			if base in alphabet:
				if (base in ["G","C"]):
					gc=gc+1
					profile["GC"].setdefault(species,[]).append(1)
				else:
					profile["GC"].setdefault(species,[]).append(0)
				if (base in ["G","T"]):
					profile["GT"].setdefault(species,[]).append(1)
				else:
					profile["GT"].setdefault(species,[]).append(0)
				if (base in ["G","A"]):
					profile["GA"].setdefault(species,[]).append(1)
				else:
					profile["GA"].setdefault(species,[]).append(0)
				if (base in ["A","C"]):
					profile["AC"].setdefault(species,[]).append(1)
				else:
					profile["AC"].setdefault(species,[]).append(0)
				if (base in ["A","T"]):
					profile["AT"].setdefault(species,[]).append(1)
				else:
					profile["AT"].setdefault(species,[]).append(0)
				if (base in ["T","C"]):
					profile["TC"].setdefault(species,[]).append(1)
				else:
					profile["TC"].setdefault(species,[]).append(0)
			else:
				profile["GC"].setdefault(species,[]).append(-1000)
				profile["GT"].setdefault(species,[]).append(-1000)
				profile["GA"].setdefault(species,[]).append(-1000)
				profile["AC"].setdefault(species,[]).append(-1000)
				profile["AT"].setdefault(species,[]).append(-1000)
				profile["TC"].setdefault(species,[]).append(-1000)
		newSeqs.append(gc)
		if getPosition(site)==3:
			for species in seq.keys():
				codon="".join(seq[species][site-2:site+1])
				if codon in genetic.keys():
					aa=genetic[codon]
					for i in range(3):
						for name in properties.keys(): profile[name].setdefault(species,[]).append(properties[name][aa])
				else:
					for i in range(3):
						for name in  properties.keys(): profile[name].setdefault(species,[]).append(-1000)
	return profile

# ignore unknown character by assigning the mean value of the windows to it
def cleanData(matrix):
	matrix2={}
	nspecies=len(matrix.keys())
	nSites=len(matrix[matrix.keys()[0]])
	mean=[]
	for i in range(nSites):
		thisMean=[]
		for species in matrix.keys():
			siteX=matrix[species][i]
			if siteX!=-1000: thisMean.append(siteX)
		mean.append(getMean(thisMean))
	for species in matrix.keys():
		matrix2[species]=[]
		for i in range(nSites):
			siteX=matrix[species][i]
			if siteX!=-1000:
				matrix2[species].append(siteX)
			else:
				matrix2[species].append(mean[i])
	return matrix2

# F-test (DFS)
def dfsScore(matrix):
	meanSp={}
	nspecies=len(matrix.keys())
	nSites=len(matrix[matrix.keys()[0]])
	matrix=cleanData(matrix)
	ccSp={} # if the data is perfect then ccSp[species]=windowSize
	ccp_within=0.0 # if the data is perfect then ccp_within=(nspecies-1)*windowSize
	ccp=0 #if the data is perfect, then ccp=nspecies
	for species in matrix.keys():
		meanSp[species]=0
		ccSp[species]=0
		for i in range(nSites):
			siteX=matrix[species][i]
			if siteX!=-1000:
				meanSp[species]=meanSp[species]+siteX
				ccSp[species]=ccSp[species]+1
		if ccSp[species]>0:
			meanSp[species]=float(meanSp[species])/float(ccSp[species])
			ccp=ccp+1
	meanT=float(sum(meanSp.values()))/float(ccp)
	between=0.0
	within=0.0
	for species in matrix.keys(): 
		for i in range(nSites):
			siteX=matrix[species][i]
			if siteX!=-1000:
				within=within+math.pow(siteX-meanSp[species],2)
				ccp_within=ccp_within+1
		if ccSp[species]>0: between=between+math.pow(meanSp[species]-meanT,2)
	MSb=0
	MSw=0
	if ccp-1!=0: MSb=float(between)/float(ccp-1)
	if ccp_within-1!=0: MSw=float(within)/float(ccp_within-1)
	thisDFS=0.0
	if float(MSw)!=0: thisDFS=float(MSb)/float(MSw)
	return thisDFS

# extract windows
def DFS(matrix,nSites,speciesList,windowSize):
	dfs=[]
	for site in range(nSites-windowSize):
		thisDFS=dfsCompare(matrix,[x for x in range(site,site+windowSize)])
		dfs.append(thisDFS)
	return dfs

# setup matrix for F-test (DFS)
def dfsCompare(matrix,listOfSites):
	matrix2={}
	for species in matrix.keys():
		for sampledSite in listOfSites:
			matrix2.setdefault(species,[]).append(matrix[species][sampledSite])
	return dfsScore(matrix2)


# deault values and attributes
percent=0.05
resample=100
windowSize=21
toDo=["H"]
modify=False
genetC="universal"
outputFileName="DruidsOutput"

properties={
"H":{"A":1.8 ,"C":2.5  ,"D":-3.5 ,"E":-3.5 ,"F":2.8  ,"G":-0.4,"H":-3.2 ,"I":4.5  ,"K":-3.9 ,"L":3.8  ,"M":1.9  ,"N":-3.5 ,"P":-1.6 ,"Q":-3.5 ,"R":-4.5 ,"S":-0.8,"T":-0.7 ,"V":4.2  ,"W":-0.9 ,"Y":-1.3 ,"*":0},
"V":{"A":88.6,"C":108.5,"D":111.1,"E":138.4,"F":189.9,"G":60.1,"H":153.2,"I":166.7,"K":168.6,"L":166.7,"M":162.9,"N":114.1,"P":112.7,"Q":143.8,"R":173.4,"S":89.0,"T":116.1,"V":140.0,"W":227.8,"Y":193.6,"*":0},
"C":{"A":0   ,"C":0    ,"D":-1   ,"E":-1   ,"F":0    ,"G":0   ,"H":1    ,"I":0    ,"K":1    ,"L":0    ,"M":0    ,"N":0    ,"P":0    ,"Q":0    ,"R":1    ,"S":0   ,"T":0    ,"V":0    ,"W":0    ,"Y":0    ,"*":0}
}
geneticCode={"mammalMt":{"GCT":"A","GCA":"A","GCG":"A","GCC":"A","TGT":"C","TGC":"C","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TTT":"F","TTC":"F","GGT":"G","GGA":"G","GGG":"G","GGC":"G","CAC":"H","CAT":"H","ATC":"I","ATT":"I","AAG":"K","AAA":"K","CTT":"L","CTA":"L","CTG":"L","CTC":"L","TTA":"L","TTG":"L","ATG":"M","ATA":"M","AAC":"N","AAT":"N","CCT":"P","CCA":"P","CCG":"P","CCC":"P","CAG":"Q","CAA":"Q","CGG":"R","CGA":"R","CGC":"R","CGT":"R","AGT":"S","AGC":"S","AGA":"*","AGG":"*","TCC":"S","TCG":"S","TCA":"S","TCT":"S","ACG":"T","ACA":"T","ACC":"T","ACT":"T","GTG":"V","GTA":"V","GTC":"V","GTT":"V","TGG":"W","TGA":"W","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*"},"universal":{"GCT":"A","GCA":"A","GCG":"A","GCC":"A","TGT":"C","TGC":"C","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TTT":"F","TTC":"F","GGT":"G","GGA":"G","GGG":"G","GGC":"G","CAC":"H","CAT":"H","ATC":"I","ATT":"I","ATA":"I","AAG":"K","AAA":"K","CTT":"L","CTA":"L","CTG":"L","CTC":"L","TTA":"L","TTG":"L","ATG":"M","AAC":"N","AAT":"N","CCT":"P","CCA":"P","CCG":"P","CCC":"P","CAG":"Q","CAA":"Q","CGG":"R","CGA":"R","CGC":"R","CGT":"R","AGT":"S","AGC":"S","AGA":"R","AGG":"R","TCC":"S","TCG":"S","TCA":"S","TCT":"S","ACG":"T","ACA":"T","ACC":"T","ACT":"T","GTG":"V","GTA":"V","GTC":"V","GTT":"V","TGG":"W","TGA":"*","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*"}}
labels={"GT":"G-T","GA":"G-A","AT":"A-T","AC":"A-C","TC":"T-C","GC":"G-C","H":"Hydrophobicity","V":"Volume","C":"Charge"}
alphabet=["A","C","G","T"]

# MAIN PROGRAM
# arguments and options parser
try:
	opts,args=getopt.getopt(sys.argv[1:],"f:w:a:o:b:p:g:m:",['fileName','windowSize',"attribute","outputFile","bootstrap","percent","geneC","modify"])
except getopt.GetoptError:
	print "Usage Error"
	sys.exit(2)
for opt,arg in opts:
	if opt in ("-f", "--fileName"):
		fileName=arg
	if opt in ("-w", "--windowSize"):
		windowSize=max(3,int(arg))
	if opt in ("-a", "--attribute"):
		toDo=[]
		names=arg.split(",")
		for name in names:
			if name in labels.keys():
				toDo.append(name)
			else:
				print "Property "+name+" unknown"
		if len(toDo)==0:
			print "no property were recognized: H used by default"
			toDo=["H"]
	if opt in ("-o", "--outputFile"):
		outputFileName=arg
	if opt in ("-b", "--bootstrap"):
		resample=int(arg)
	if opt in ("-p", "--percent"):
		percent=float(arg)
	if opt in ("-g", "--geneC"):
		if arg in geneticCode.keys():
			genetC=arg
		else:
			print arg+" not recognized, universal genetic code used by default"
	if opt in ("-m", "--modify"):
		if arg.upper() in ["TRUE","FALSE"]:
			if arg.upper()=="TRUE": modify==True
			if arg.upper()=="FALSE": modify==False
		else:
			print arg+" not recognized, False used by default for modification flag"

# introduction
if fileName=="":
	print "Please specify a file name"
else:
	print "\npyDRUIDS June 2007 Olivier Fedrigo and Gavin Naylor"
	print "File name: "+fileName
	print "Output file: "+outputFileName
	print "Properties: "+",".join([labels[x] for x in toDo])
	print "Window size: " +str(windowSize)
	print "Bootstrap replicates: " +str(resample)
	print "Percentage significance: " +str(percent)
	print "Genetic code: "+genetC
	print "Modification: "+str(modify)
	print "\n"

# read data file, requires biopython
seq={}
if fileName.endswith('nex') or fileName.endswith('nexus'):
	#sequences= Nexus.Nexus()
	#sequences.read(fileName)
	handle = open(fileName, "r")
	seq = SeqIO.to_dict(SeqIO.parse(handle, "nexus"))
	handle.close()
if fileName.endswith('phy') or fileName.endswith('phylip'):
	handle = open(fileName, "r")
	seq = SeqIO.to_dict(SeqIO.parse(handle, "phylip"))
	handle.close()
if fileName.endswith('fasta'):
	handle = open(fileName, "r")
	seq = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()

# populate sequences hash from the data read with different formats
sequences={}
for species in seq.keys():
	sequences[species]=seq[species].seq
nchar=len(sequences[sequences.keys()[0]])
ntaxa=len(sequences.keys())
print str(nchar)+" sites analyzed for "+str(ntaxa)+" taxa"

windowSize=min(windowSize,nchar/2) # limit windosize to thalf the data set length

if modify==True: # modify the data by the most common base if the base is not A, C , G or T
	mostCommon="A"
	maximum=0
	alpha={}
	for site in range(nchar):
		for species in sequences.keys():
			thisSequence=sequences[species]
			base=thisSequence[site]
			alpha.setdefault(base,[]).append(1)
		for base in alpha.keys():
			if len(alpha[base])>maximum:
				maximum=len(alpha[base])
				mostCommon=base
		for species in sequences.keys():
			thisSequence=sequences[species]
			base=thisSequence[site]
			if (base in alphabet)==False:
				base=mostCommon
			sequences[species][base]=base


profile=getProfiles(alphabet,sequences,properties,geneticCode[genetC]) # get attributes' profile per site (all 3 bases get the same value if it is an amino acid property)
dfs={}
thresh={}
resultsWindows={}
for name in toDo: #loop through attributes
	print labels[name]
	dfs[labels[name]]=DFS(profile[name],nchar,ntaxa,windowSize) # calculate F-test in a sliding window
	thresh[labels[name]]=bootstrap(profile[name],resample,windowSize,percent,dfsCompare,nchar) # codon bootstrap the dataset and get the F-test threshhold for a 0.05 cutoff
	resultsWindows[labels[name]]=convertToString(makeWindow(getSig(dfs[labels[name]],windowSize,thresh[labels[name]]))) # get the list of regions that deviates from stationarity for a given attribute

# save output
saveWindows(outputFileName+"_windows.out",resultsWindows)
saveList(outputFileName+"_DFS.out",dfs,thresh)
