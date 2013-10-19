Druids
======

This software allows to identify sequence regions of multiple protein-coding gene alignments that exhibit local deviation from stationarity for different physico-chemical amino acid properties and that can be responsible for misleading phylogenetic inferences.

Please cite:
Fedrigo, O., D.C. Adams, and G.J.P. Naylor. (2005) DRUIDS – Detection of Regions with Unexpected Internal Deviation from Stationarity.  J. Exp. Zoolog. Part B Mol. Dev. Evol. 304(2): 119-128.

Three versions are available:
<ul>
<li> Windows vesion (Druids_Setup.exe), not tested after windows XP</li>
<li> Mac OSX version (Druids.dmg), not tested after OSX 10.5</li>
<li> Python version (pyDruids.zip), command line version. See below</li>
</ul>

pyDruids.py -f <inputFileName> -w <windowSize> -a <attribute> -o <outputFileName> -b <nreplicates> -p <percent> -g <geneticCode> -m <fillGaps> 

requires Bio Python to be installed.

inputFileName
	any of the following formats: phylip, fasta, nexus
	(should ends with .phy or .phylip or .fasta or .nex or.nexus)
	bases should be all upper case: A ,C ,G or T
	assume your data is in frame and its length is /3 multiple

windowSize
	minimum window size=3
	maximum wiwdow size= half the the sequence length
	default= 21

attribute
	comma separated attributes
	any of  GA,GT,AT,AC,TC,GC for %base content
	H: hydrophobicity
	V: colume
	C: charge
	default= H

outputFileName
	will output:
		<outputFileName>_windows.out: only significant windows
		<outputFileName>_DFS.out: all the values
	default= DruidsOutput

nreplicates
	number of codon bootstrap replicates
	default= 100

percent
	significant cutoff in the resampled distribution to call significant windows
	default= 0.05

geneticCode
	for now, are only available: mammalMt and universal
	default= universal

fillGaps
	FALSE: it will assign the mean attribute value of the window to this site for unknown characters (e.g. gaps, lower case, N)
	TRUE: it will replace gaps by the most common base if the character is unknown
	default= FALSE
