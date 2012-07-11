#!/usr/bin/env pythons

import os, sys, argparse, subprocess
import datetime, random, string
import numpy as ny
import pylab as pl

#################################################################
	
def eventsStats(inputFileName, outputDir, percentList=[1], reuse=False):
	print "\n>> start:eventsStats - full file size=", os.path.getsize(inputFileName), "bytes"
	
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	
	statsDir = os.path.join(outputDir, "stats/")
	if not os.path.exists(statsDir):
		os.makedirs(statsDir)
	#~ print "statsDir = "+statsDir
	
	sampleFileNames, statsFileNames = [], []
	for pct in percentList:
		
		## determine and store data file name
		sampleFileName = inputFileName + "-%" + string.replace(str(pct), '.', '-') if pct != 1 else inputFileName
		sampleFileNames.append(sampleFileName)
		
		## go through file + create sample if required
		nbEntries, nbFeatures = browseEvents(sampleFileName, pct, reuse, inputFileName)
		
		## write statistics
		statsFileName = statsDir + os.path.basename(sampleFileName) + ".stats.events"
		statsFileNames.append(statsFileName)
		
		if verbose:
			print "   Writing statistics in " + statsFileName + "\n"
		open(statsFileName, 'w').write(\
			str(os.path.getsize(sampleFileName)) + " Size_Feature_Entries\n" +\
			str(nbEntries) + " Total_Entries\n" +\
			str(nbFeatures) + " Total_Features\n" +\
			str(1.0 * nbFeatures/nbEntries) + " Average_Features\n")
		
	print ">> end:eventsStats"
	return sampleFileNames, statsFileNames


def browseEvents(sampleFileName, pct, reuse=False, inputFileName=None):
	if verbose:
		print "\n   >> start:browseEvents"
	## initialise variables for stats
	nbEntries, nbFeatures = 0, 0
	
	## go through file
	if pct == 1 or reuse and os.path.isfile(sampleFileName):
		if verbose:
			print "   Reusing file for", 100.0*pct, "%"
		sampleFile = open(sampleFileName, 'r')
		
		for entry in sampleFile:
			## statistics
			nbEntries += 1
			nbFeatures += len(entry.split()) -1
			
	else:
		if verbose:
			print "   Creating file for", 100.0*pct, "%"
		sampleFile = open(sampleFileName, 'w')
		for entry in open(inputFileName, 'r'):
			if random.random() < pct: 	
				## new sample file
				sampleFile.write(entry)
				## statistics
				nbEntries += 1
				nbFeatures += len(entry.split()) -1
		sampleFile.close()
		
	if verbose:
		print "  ", nbEntries, "entries and", nbFeatures, "features"
		print "   >> end:browseEvents"
	return nbEntries, nbFeatures


def bybloStats(sampleFileNames, outputDir, bybloDir, bybloParams="", reuse=False, verbose=False):
	print "\n>> start:bybloStats"
	
	## prepare for Byblo output
	thesauriDir = os.path.join(outputDir, "thesauri/")
	if not os.path.exists(thesauriDir):
		os.makedirs(thesauriDir)
	statsDir = os.path.join(outputDir, "stats/")
	statsFileNames = []
	
	#~ print "statsDir = "+statsDir
	#~ print "thesauriDir = "+thesauriDir
	
	for fileName in sampleFileNames:
		
		## run Byblo for this sample file
		runTime = runByblo(os.path.abspath(fileName), os.path.abspath(thesauriDir), \
			os.path.abspath(bybloDir), bybloParams, verbose)
		
		## write statistics
		thesaurusName = thesauriDir + os.path.basename(fileName)
		statsFileName = statsDir + os.path.basename(fileName) + ".stats.byblo"
		statsFileNames.append(statsFileName)
		statsFile = open(statsFileName, 'w')
		
		for suffix in [".sims.neighbours.strings",".events.filtered",".entries.filtered"]:
			size = os.path.getsize(thesaurusName + suffix)
			nbLines = sum(1 for line in open(thesaurusName + suffix, 'r'))
			statsFile.write(str(size) + " File_Size_Bytes" + suffix + '\n' +\
				str(nbLines) + " Lines_in_file" + suffix + '\n')
		statsFile.write(str(runTime) + " Byblo_Run_Time\n")
		
	print ">> end:bybloStats"
	return statsFileNames


def runByblo(inputFileName, outputDir,  bybloDir, bybloParams, verbose=False):
	if verbose:
		print "\n   >> start:singleBybloRun "
		print  "Input =", inputFileName
		print  "Output=", outputDir
	
	## move to BYblo directory, start timer and run it in a subprocess
	startDir=os.path.abspath(os.getcwd())
	logFile = open(os.devnull, 'w') if not verbose else None
	os.chdir(bybloDir)
	if verbose:
		print "  Moved to " + os.getcwd()
	
	stime = datetime.datetime.now()
	out = subprocess.call(os.path.abspath("./byblo.sh ") + " -i " + inputFileName + " -o " + outputDir +\
		" "+ bybloParams, shell = True, stdout = logFile, stderr = logFile)
	etime = datetime.datetime.now()
	
	if  logFile != None:
		logFile.close()
	os.chdir(startDir)
	if verbose:
		print "  Moved back to " + os.getcwd()
	
	## fail?
	if(not out == 0):
		print "   Byblo failed on input file: " + inputFileName + "\n   Fail Code: " + str(out)
		sys.exit()
		
	runTime = 1.0*(etime - stime).total_seconds() #/ 3600
	if verbose:
		print "   Run time =", runTime
		print "   >> end:singleBybloRun \n"
	return runTime
	
def generateHistograms(sampleFileNames, outputDir, verbose):
	print "\n>> start:generateHistograms "
	
	thesauriDir = os.path.join(outputDir, "thesauri/")
	graphsDir = os.path.join(outputDir, "graphs/")
	if not os.path.exists(graphsDir):
		os.makedirs(graphsDir)
	
	#~ print "thesauriDir = "+thesauriDir
	#~ print "graphsDir = "+graphsDir
	
	for fileName in sampleFileNames:
		for suffix in ['.entries.filtered','.events.filtered']:
			createOccurenceHistogram(thesauriDir, os.path.basename(fileName), suffix, graphsDir)
		
		suffix = '.sims.neighbours'
		createSimilarityHistogram(thesauriDir, os.path.basename(fileName), suffix, graphsDir)
	
	print ">> end:generateHistograms"


def createOccurenceHistogram(directory, fileName, suffix, outputDir, verbose=False):
	if verbose:
		print "\n     >> start:createOccurenceHistogram "+ suffix
	
	bins = ny.arange(1, 1000, 10)
	values = [int( line.split()[1] ) for line in open(directory + fileName + suffix, 'r')]
	#~ print "out = ", out
	createHistogram("Occurence", values, bins, outputDir, fileName+suffix)
	
	if verbose:
		print "     >> end:createOccurenceHistogram\n"

def createSimilarityHistogram(directory, fileName, suffix, outputDir, verbose=False):
	if verbose:
		print "\n     >> start:createSImilarityHistogram "+ suffix

	bins = ny.arange(0.0,1.0,0.1)
	values=[]
	for line in open(directory + fileName + suffix,'r'):
		values += [value for i, value in enumerate(line.split()[2:]) if i%2 == 0 ]
	#~ print "out = ", out
	createHistogram("Similarity", values, bins, outputDir, fileName+suffix)
	
	if verbose:
		print "     >> end:createSimilarityHistogram\n"


def createHistogram(label, values, bins, outputDir, fileName):
	xlabel, ylabel = label+" values", "Frequency"
	pl.xlabel(xlabel)
	pl.ylabel(ylabel)
	pl.title(label+' histogram for: ' + fileName)
	pl.hist(values, bins=bins)
	pl.savefig(outputDir + xlabel + '-' + ylabel + '-' + fileName + '.pdf')


def statsToDictionary(fileNames):
	dictList = []
	## create dictionary for stats in each file
	for f in sorted(fileNames):
		dict = {}
		for line in open(f, 'r'):
			fields = line.split()
			dict[fields[1]] = float(fields[0])
		dictList.append(dict)	
	
	## for each key, list in tuple all values (in all dictionaries)  for this key
	finalDict = {}
	for dict in dictList:
		for k in dict.iterkeys():
			finalDict[k] = (tuple(d[k] for d in dictList if k in d))
		
	return finalDict
	
def generatePlots(dictOfPlots, outputDir, verbose=False, yAxis="", xAxis=""):
	print "\n >> start:generatePlots "
	
	graphsDir = os.path.join(outputDir, "graphs/")
	## for a constant value on the yAxis
	if(not yAxis == ""):
		for key, value in dictOfPlots.iteritems():
			if(not yAxis == key):
				createPlot(value, dictOfPlots[yAxis][1], key, yAxis, graphsDir)
	## for a constant value on the xAxis
	elif(not xAxis == ""):
		for key, value in dictOfPlots.iteritems():
			if(not xAxis == key):
				createPlot(value, dictOfPlots[xAxis][1], key, xAxis, graphsDir)
	## no constant value (plot everything)
	else:
		for key, value in dictOfPlots.iteritems():
			for key2, value2 in dictOfPlots.iteritems():
				if(not key == key2):
					if verbose:
						print "Plot " + key + " vs. " + key2
					createPlot(value, value2, key, key2, graphsDir)
	
	print " >> end:generatePlots"

def createPlot(x, y, xlabel, ylabel, outputDirectory):
        fig = pl.figure(xlabel + '-' + ylabel)
        pl.plot(x,y)
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        pl.savefig(outputDirectory + xlabel+'-'+ylabel+'.pdf')



"""
def deleteOnExit():
        ## If deleteOnExit is "True" then all files except the output are deleted
        if(self.deleteOnExit):
            directory = self.listdir_nohidden(self.mainDirectory)
            for d in directory:
                if not (d.startswith("output")):
                    os.system("rm -r " + self.mainDirectory + d)

        print "Output can be found here: " + outputDirectory
        print "Process complete!"
"""

def print_lines(list, min=0, max=None, line_max=None, title="List"):
	if max is None:
		max = len(list)
	print "\n" + title
	
	if type(list) is dict:
		for i, key in enumerate(list.iterkeys()) :
			if i<max:
				print key, list[key]
			else:
				break
	else:
		for index in range(min, max):
			print str(list[index])[:line_max], "..."
	print "\n"


if __name__=='__main__':
	
	##############################################################################
	parser = argparse.ArgumentParser(description='Generate statistics for Byblo.')

	# data file name
	parser.add_argument('data', metavar='file', nargs=1, \
		action='store', 
		help='data file containing  a feature / entry set')
	# percentages of entries to take into account
	parser.add_argument('-p', '--percentage', metavar='x', type=float, nargs='*', dest='p', \
		action='store', default=None,
		help='approximate percentage (0-1) of entries to keep  for a set of statistics')
	# output directory
	parser.add_argument('-O', '--output-dir', metavar='path', nargs=1, dest='O', \
		action='store', default=None,
		help='output directory where "stats", "thesauri" and "graphs" sub-directories might be created '+\
		'(default: ".")')
	# Byblo directory
	parser.add_argument('-B', '--Byblo-dir', metavar='path', nargs=1, dest='B', \
		action='store', default=None,
		help='Byblo directory, usually named "Byblo-x.x.x" '+\
		'(default: "../Byblo-2.0.1")')
	# BYblo parameters
	parser.add_argument('-P', '--Byblo-params', metavar='string', nargs=1, dest='P', \
		action='store', default=None,
		help='Byblo parameters, read Byblo doc for more info '+\
		'(default: "-fef 10 -fff 10 -t 6 -Smn 0.1")')
	# reuse option
	parser.add_argument('-r', '--reuse', dest='reuse', action='store_const', \
		const=True, default=False, 
		help='reuse sample files previously created')
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_const', \
		const=True, default=False, 
		help='display Byblo output and statistics information for each step')
	
	args = parser.parse_args()
	##############################################################################
	
	pctList = [p for p in args.p if p != 0 and p<=1] if args.p != None else [1]
	outputDir = args.O[0] if args.O != None else "."
	bybloDir = args.B[0] if args.B !=None else "../Byblo-2.0.1"
	bybloParams = args.P[0] if args.P != None else "-fef 10 -fff 10 -t 6 -Smn 0.1"
	reuse = args.reuse
	verbose = args.verbose
	
	stime = datetime.datetime.now()
	print "***************************************************************************"
	print "BYBLO STATISTICS TOOL"
	print "***************************************************************************\n"
	
	## EVENTS STATS + sample file creation when necessary
	sampleFileNames, statsFileNames = eventsStats(args.data[0], outputDir, pctList, reuse)
	if verbose:
		print_lines(statsFileNames, title="Statistics files after feature extraction")
	
	## BYBLO STATS
	statsFileNames += bybloStats(sampleFileNames, outputDir, bybloDir, bybloParams, reuse, verbose)
	if verbose:
		print_lines(statsFileNames, title="Statistics files after byblo run")
	
	## HISTOGRAMS
	generateHistograms(sampleFileNames, outputDir, verbose)
	
	## PLOTS
	statsDict = statsToDictionary(statsFileNames)
	if verbose:
		print_lines(statsDict, title="Table containing the statistics")
	generatePlots(statsDict, outputDir, verbose)
	
	
	etime = datetime.datetime.now()
	print "\n>Execution took", etime-stime, "hours"
	
	# pause at the end
	#raw_input()

	    