#! /usr/bin/env python
__author__="Andres Perez-Figueroa"
__date__ ="$20-mar-2012 13:18:20$"
__version_="1.1"
# Converter from expanded FSTAT format (like in quantiNemo) to METAPOP2
# For 3 digit without sex

from sys import exit
import sys


if __name__ == "__main__":
	#print "F2M 1.0 - File conversor from expanded FSTAT to Metapop"
#    print " "
	if len(sys.argv)<2:
		print "Argument missing! USAGE: python %s <inputFile>" % sys.argv[0]
 		sys.exit(-1)

	inputfile=sys.argv[1]
	
	indcode = []
	indsex = []
	indPop = []
	popSize = []
	popCode = []
	indgeno1 = []
	indgeno2 = []
	hashPop = []
	reg=[]
	indiv=[]
	indivT=[]
	lociName = []

	homolog=False
	try:
		
		handle = open(inputfile)
		#print inputfile
		pop=-1
		index=0
		
		temp_i=0
		temp_p=0
		popcount=0
		indexLine=1
    	    
	
		for line in handle.xreadlines():
			
			if indexLine==1:
				#First line: npops  nloci  maxall  digits
				
				fline=line.split()
				
				pops=int(fline[0])
				loci=int(fline[1])	
				start=loci+1
			#	indexLine=indexLine+1
				
				#continue
			if indexLine>start:
				data=line.split()
				
				Pop=int(data[0])
				
				geno=[]
				#sex=int(data[loci+2])
				sex=1
				for i in xrange(1,loci+1):
					geno.append((data[i][:3]) ) #
					geno.append((data[i][3:]) )
				#print geno[0]
				if sex==0:
					sex=1
				else:
					sex=0
				indiv.append((Pop, sex, geno) )
	
	
				index = index+1
			indexLine=indexLine+1
	
	


	except IOError:
		print "Can not open file " + inputfile + " to read. Exiting"
		exit(2)
	except:
#        print "Error"
#        print homolog, indexLine
		exit(3)

	#print index
	#print indiv[0]
	#print indiv[0][2][2]
	popSize=[0]*pops
	for i in xrange(index):
		
		popSize[indiv[i][0]-1]+=1
	#print popSize
	
	rind=[]
	cont=0
	for x in xrange(index):
	    #print x, indiv[x][0]
		rind.append(str(x+1)+" ["+str(index)+"] "+ str(indiv[x][1]))
		for l in xrange(loci):

			rind[x]=rind[x]+" "+str(indiv[x][2][(l*2)])+" "+str(indiv[x][2][(l*2)+1])
	#print rind





	#loci=len(reg)
	#pops=7
	#
	#    print""
	metapopfile=inputfile+".mtp"	
	try:
		f = open(metapopfile, "w")
	except IOError:
		print "Can not open file " + metapopfile + " to write. Exiting"
		exit(2)
	#
	f.write("n\t%d\n" % (pops))
	f.write("nloci\t%d\n" % (loci))
	for p in xrange(pops):
		f.write("sub%d\n" % (p+1))
		f.write("N\t%d\n" % (popSize[p]))
		f.write("Ne\t%d\n" % (popSize[p]))
		f.write("newN\t%d %d\n" % (popSize[p]/2, popSize[p]/2))
		for i in xrange(popSize[p]):
			f.write("%s\n" % (rind[cont]))
			cont+=1
	f.close()

