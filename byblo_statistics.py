#!/usr/bin/env pythons

import os, sys, argparse, subprocess
import datetime, random, string

#################################################################
	
def eventsStats(inputFileName, percentList=[1], reuse=False):
	print "\n>> start:eventsStats - full file size=", os.path.getsize(inputFileName), "bytes"
	
	sampleFileNames, statsFileNames = [], []
	if not os.path.exists("./stats"):
		os.makedirs("./stats")
		
	for pct in percentList:
		
		## determine and store data file name
		sampleFileName = inputFileName + "-%" + string.replace(str(pct), '.', '-') if pct != 1 else inputFileName
		sampleFileNames.append(sampleFileName)
		
		## go through file + create sample if required
		nbEntries, nbFeatures = browseEvents(sampleFileName, pct, reuse, inputFileName)
		
		## write statistics
		statsFileName = "./stats/" + os.path.basename(sampleFileName) + ".stats.events"
		statsFileNames.append(statsFileName)
		
		print "   Writing statistics in " + statsFileName + "\n"
		open(statsFileName, 'w').write(\
			str(os.path.getsize(sampleFileName)) + " Size_Feature_Entries\n" +\
			str(nbEntries) + " Total_Entries\n" +\
			str(nbFeatures) + " Total_Features\n" +\
			str(1.0 * nbFeatures/nbEntries) + " Average_Features\n")
		
	print ">> end:eventsStats"
	return sampleFileNames, statsFileNames


def browseEvents(sampleFileName, pct, reuse=False, inputFileName=None):
	print "\n   >> start:browseEvents"
	## initialise variables for stats
	nbEntries, nbFeatures = 0, 0
	
	## go through file
	if pct == 1 or reuse and os.path.isfile(sampleFileName):
		print "   Reusing file for", 100.0*pct, "%"
		sampleFile = open(sampleFileName, 'r')
		
		for entry in sampleFile:
			## statistics
			nbEntries += 1
			nbFeatures += len(entry.split()) -1
			
	else:
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
		
	print "  ", nbEntries, "entries and", nbFeatures, "features"
	print "   >> end:browseEvents"
	return nbEntries, nbFeatures


def bybloStats(sampleFileNames, bybloParams="", reuse=False):
	print "\n>> start:bybloStats"
	
	## prepare for Byblo output
	thesauriDir = "./thesauri/"
	statsFileNames = []
	if not os.path.exists(thesauriDir):
		os.makedirs(thesauriDir)
	
	for fileName in sampleFileNames:
		
		## run Byblo for this sample file
		runTime = runByblo(os.path.abspath(fileName), os.path.abspath(thesauriDir), bybloParams)
		
		## write statistics
		thesaurusName = thesauriDir + os.path.basename(fileName)
		statsFileName = "./stats/" + os.path.basename(fileName) + ".stats.byblo"
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


def runByblo(inputFileName, outputDir,  bybloParams):
	print "\n   >> start:singleBybloRun "

	# temporary
	bybloLoc = "../Byblo-2.0.1"
	print os.path.abspath(bybloLoc+"/byblo.sh ")
	print os.path.abspath(inputFileName)
	print os.path.abspath(outputDir)
	
	stime = datetime.datetime.now()
	#fnull = open(os.devnull, 'w')
	os.chdir(bybloLoc)
	print os.getcwd()
	out = subprocess.call(os.path.abspath("./byblo.sh ") + \
		" -i " + inputFileName + " -o " + outputDir +\
		" "+ bybloParams, shell = True\
		#, stdout = fnull, stderr = fnull\
		)
	#fnull.close()
	os.chdir("../bybloStats")
	if(not out == 0):
		print "   Byblo failed on input file: " + inputFileName + "\n   Fail Code: " + str(out)
		sys.exit()
		
	
	etime = datetime.datetime.now()
	totalTime = 
	runTime = 1.0*(etime - stime).total_seconds() / 3600
	print "runTime=", runTime
	print "   >> end:singleBybloRun "
	return runTime
	
def generateHistograms(sampleFileNames):
	print "\n   >> start:generateHistograms "
	
	thesauriDir = "./thesauri/"
	# don4t forget to creae this one!
	outputDirectory = "./graphs/"
	
	for suffix in ['.entries.filtered','.events.filtered']:
		for fileName in sampleFileNames:
			self.createOccurenceHistogram(thesauriDir, os.path.basename(fileName), suffix, outputDirectory)

	for suff in ['.sims.neighbours']:
		for fileName in sampleFileNames:
			self.createSimilaritiesHistogram(thesauriDir, os.path.basename(fileName), suffix, outputDirectory)
	
	
	print "   >> end:generateHistograms"
"""
	...
	
	## Compiles the stats and builds the output plots...
	results = {}
	for loc in statsLocations:
		self.resultsTable(loc, results)

        
        self.createPlots(results, outputDirectory)
        
        for suff in suffixesOccHist:
            for pair in nameTimePairs:
                self.createOccurenceHistogram(outputLocationByblo, suff, pair[0], outputDirectory)


        for suff in suffixesSimHist:
            for pair in nameTimePairs:
                self.createSimilaritiesHistogram(outputLocationByblo, suff, pair[0], outputDirectory)
        
        ## If deleteOnExit is "True" then all files except the output are deleted
        if(self.deleteOnExit):
            directory = self.listdir_nohidden(self.mainDirectory)
            for d in directory:
                if not (d.startswith("output")):
                    os.system("rm -r " + self.mainDirectory + d)

        print "Output can be found here: " + outputDirectory
        print "Process complete!"

    ## Takes the table/dictionary of compiled stats and plots them against y and/or x.
    ## If both x and y are empty then it plots all items against themselves
    def createPlots(self, dictOfPlots, outputDirectory, yAxis="", xAxis=""):
        print "Building Plots..."
        
        if(not yAxis == ""):
            for key, value in dictOfPlots.iteritems():
                if(not yAxis == key):
                    self.createPlot(value, dictOfPlots[yAxis][1], key, yAxis)
        elif(not xAxis == ""):
            for key, value in dictOfPlots.iteritems():
                if(not xAxis == key):
                    self.createPlot(value, dictOfPlots[xAxis][1], key, xAxis)
        else:
            for key, value in dictOfPlots.iteritems():
                for key2, value2 in dictOfPlots.iteritems():
                    if(not key == key2):
                        self.createPlot(value, value2, key, key2,outputDirectory)
                        
        print "Plots built!"

    def createPlot(self, x, y, xlabel, ylabel, outputDirectory):
        fig = pl.figure(xlabel + '-' + ylabel)
        pl.plot(x,y)
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        pl.savefig(outputDirectory + xlabel+'-'+ylabel+'.pdf')

    ## Creates the histograms for the frequency counts of a given file
    def createOccurenceHistogram(self, directory, suffix, fileName, outputLocation):
        print "Building occurence histogram for: " + suffix
        inputLoc = directory + fileName + suffix
        text = open(inputLoc, 'r')
        out = []
        ran = ny.arange(1, 1000, 10)
        for line in text:
            freq = line.split()[1]
            out.append(int(freq))
        xlabel = "Occurences"
        ylabel = "Frequency"
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        pl.title('Occurence histogram for: ' + fileName + suffix)
        histo = pl.hist(out, bins=ran)
        pl.savefig(outputLocation + xlabel + '-' + ylabel + fileName + suffix + '.pdf')

        print "Occurence histogram complete!"

    def createSimilaritiesHistogram(self, directory, suffix, fileName, outputLocation):
        print "Building similarities histogram for: " + suffix
        
        inputLoc = directory + fileName + suffix
        text = open(inputLoc,'r')
        out = []
        ran = ny.arange(0.0,1.0,0.1)
        for line in text:
            spline = line.split()
            i = 2
            while i < len(spline):
                if(i % 2 == 0):
                    out.append(spline[i])
                i = i + 1
        xlabel = 'Similarity value'
        ylabel = 'Frequency'
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        pl.title('Similarity histogram for: ' + fileName + suffix)
        histo = pl.hist(out,bins=ran)
        pl.savefig(outputLocation + xlabel + '-' + ylabel + fileName + suffix + '.pdf')
        
        print "Similarities histogram complete!"

"""

def statsToTable(fileNames):
	dictList = []
	## create dictionary for stats in each file
	for f in sorted(fileNames):
		dict = {}
		for line in open(f, 'r'):
			fields = line.split()
			dict[fields[1]] = float(fields[0])
		dictList.append(dict)
	
	
	## for each key, list (as tuples) all entries (in all dictionnaries) that contain this key
	table = []
	for key in [dict.iterkeys() for dict in dictList]:
		## list tuples in all dictionaries for a given key
		table[key] = (tuple(d[k] for d in dictList if k in d))
		
	return table



def print_lines(list, min=0, max=None, line_max=None, title="List"):
	if max is None:
		max = len(list)
	print "\n" + title
	for index in range(min, max):
		print str(list[index])[:line_max], "..."
	print "\n"


if __name__=='__main__':
	
	##############################################################################
	parser = argparse.ArgumentParser(description='Generate statistics for Byblo.')

	# percentages of entries to take into account
	parser.add_argument('-p', '--percentage', type=float, nargs='*', dest='pct', \
		action='store', default=None)
	# reuse option
	parser.add_argument('-r', '--reuse', dest='reuse', action='store_const', \
		const=True, default=False)
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_const', \
		const=True, default=False)
	# file name
	parser.add_argument('file', metavar='file', type=file, nargs=1, \
		action='store', help='input file (feature / entry set)')
	
	args = parser.parse_args()
	##############################################################################
	## Validate list of size percentages
	pctList = [p for p in args.pct if p != 0 and p<=1] if args.pct != None else [1]
	## Parameters for Byblo
	bybloParams = "-fef 10 -fff 10 -t 6 -Smn 0.1"
	
	stime = datetime.datetime.now()
	print "***************************************************************************"
	print "BYBLO STATISTICS TOOL"
	print "***************************************************************************\n"
	
	## EVENTS STATS + sample file creation when necessary
	sampleFileNames, statsFileNames = eventsStats(args.file[0].name, pctList, args.reuse)
	print_lines(statsFileNames, title="Statistics files after feature extraction")
	
	## BYBLO STATS
	statsFileNames += bybloStats(sampleFileNames, bybloParams, args.reuse)
	print_lines(statsFileNames, title="Statistics files after byblo run")
	
	## STATS PLOTTING
	
	
	
	
	etime = datetime.datetime.now()
	print ">Execution took", etime-stime, "hours"
	
	# pause at the end
	#raw_input()

	    