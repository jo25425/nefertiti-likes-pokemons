
import os, sys, argparse, subprocess
import datetime, random, string

#################################################################
	
def eventsStats(inputFileName, percentList=[1], reuse=False):
	print "\n>> start:eventsStats - full file size=", os.path.getsize(inputFileName), "bytes"
	
	sampleFileNames, statsFileNames = [], []
	if not os.path.exists(".\\stats"):
		os.makedirs(".\\stats")
		
	for pct in percentList:
		
		## determine and store data file name
		sampleFileName = inputFileName + "-sample" + string.replace(str(pct), '.', '-') if pct != 1 else inputFileName
		sampleFileNames.append(sampleFileName)
		
		## go through file + create sample if required
		nbEntries, nbFeatures = browseEvents(sampleFileName, pct, reuse, inputFileName)
		
		## write statistics
		statsFileName = ".\\stats\\" + sampleFileName[string.rfind(sampleFileName, '\\'):] + ".stats.events"
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
	if reuse and os.path.isfile(sampleFileName):
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
	thesauriDir = ".\\thesauri\\"
	if not os.path.exists(thesauriDir):
		os.makedirs(thesauriDir)
		
	for fileName in sampleFileNames:
		
		## initialise variables for stats
		
		## run Byblo for this sample file
		runTime = runByblo(fileName, thesauriDir, bybloParams)
		
		## write statistics
		thesaurusName = thesauriDir + sampleFileName[string.rfind(sampleFileName, '\\'):]
		statsFileName = ".\\stats\\" + sampleFileName[string.rfind(sampleFileName, '\\'):] + ".stats.byblo"
		statsFileNames.append(statsFileName)
		statsFile = open(statsFileName, 'a')
		
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

	stime = datetime.datetime.now()
	#fnull = open(os.devnull, 'w')
	bybloLoc = "..\\Byblo-2.0.1"
	print bybloLoc+"\\byblo.sh "
	print os.path.abspath(inputFileName)
	print os.path.abspath(outputDir)
	out = subprocess.call(bybloLoc+"\\byblo.sh " + "-i " + os.path.abspath(inputFileName) + " -o " + os.path.abspath(outputDir) +\
		" "+ bybloParams, shell = True\
		#, stdout = fnull, stderr = fnull\
		)
	#fnull.close()
	if(not out == 0):
		print "   Byblo failed on input file: " + inputFileName + "\n   Fail Code: " + str(out)
		sys.exit()
		
	
	etime = datetime.datetime.now()
	totalTime = timeFinish - timeStart
	# use implicit conversion!!!
	runTime = float(totalTime.total_seconds())/float(3600)
	print "   >> end:singleBybloRun "
	return runTime
	

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
	
	## BYBLO STATS + sample file creation when necessary
	statsFileNames += bybloStats(sampleFileNames, bybloParams, args.reuse)
	print_lines(statsFileNames, title="Statistics files after byblo run")
	

	etime = datetime.datetime.now()
	print ">Execution took", etime-stime, "hours"
	
	# pause at the end
	#raw_input()

	    