#!/usr/bin/env python -c

import os, sys, argparse, subprocess
import datetime, random, string
import numpy as np
import pylab as pl

from os.path import\
	basename as basename, \
	exists as exists, \
	getsize as getsize, \
	join as join, \
	abspath as abspath, \
	isfile as isfile
from string import find as contains

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import *
from scipy import optimize, stats

#################################################################

## Determines the set of sample files to generate statistics for using an original input file and a list of percentages,
## creates these files when required, then goes through all of them and writes some of their carateristics (size, 
## entries, events) in statistics files put in the subfolder "stats"
## @return list of sample files, list of statistics files
def eventsStats(inputFileName, outputDir, percentList=[100], reuse=[], verbose=False):
	print "\n>> start:eventsStats - full file size=", getsize(inputFileName), "bytes"
	
	if not exists(outputDir):
		os.makedirs(outputDir)
	
	statsDir = join(outputDir, "stats")
	if not exists(statsDir):
		os.makedirs(statsDir)
	
	sampleFileNames, statsFileNames = [], []
	for pct in percentList:
		
		## determine and store data file name for input
		sampleFileName = inputFileName + "-%" + string.replace(str(pct), '.', '-') if pct != 100 else inputFileName
		sampleFileNames.append(sampleFileName)
		
		## determine and store statistics file name for output
		statsFileName = join(statsDir, basename(sampleFileName) + ".stats.events")
		statsFileNames.append(statsFileName)
		
		if "events_stats" in reuse and isfile(statsFileName):
			if verbose:
				print "   Reusing statistics file for", pct, "%"
		else:
			if verbose:
				print "   Creating statistics file for", pct, "%\n"

			## go through file + create sample if required
			nbEntries, nbEvents = browseEvents(sampleFileName, pct, inputFileName, reuse,  verbose)
			
			## write statistics
			if verbose:
				print "   Writing statistics in " + statsFileName + "\n"
			open(statsFileName, 'w').write(\
				str(getsize(sampleFileName)) + " Size_In_Bytes_Of_Input_File\n" +\
				str(nbEntries) + " Total_Number_Of_Distinct_Entries\n" +\
				str(nbEvents) + " Total_Number_Of_Distinct_Events\n" +\
				str(1.0 * nbEvents/nbEntries) + " Average_Number_Of_Events_By_Entry\n")
		
	print ">> end:eventsStats"
	return sampleFileNames, statsFileNames

## Counts entries (lines) and events (fields) for the specified sample file and, if required, at the same time creates
## this sample with a size approximately equal to the specified percentage of the original file,
## @return number of entries, number of events
def browseEvents(sampleFileName, pct, inputFileName=None, reuse=[], verbose=False):
	if verbose:
		print "   >> start:browseEvents"
	## initialise variables for stats
	nbEntries, nbEvents = 0, 0
	
	## go through the file
	if pct == 100 or "samples" in reuse and isfile(sampleFileName):
		if verbose:
			print "      Reusing file for", pct, "%"
		sampleFile = open(sampleFileName, 'r')
		
		for entry in sampleFile:
			## statistics
			nbEntries += 1
			nbEvents += len(entry.split()) -1
			
	else:
		if verbose:
			print "      Creating file for", pct, "%"
		sampleFile = open(sampleFileName, 'w')
		for entry in open(inputFileName, 'r'):
			if random.random() < pct/100.: 	
				## new sample file
				sampleFile.write(entry)
				## statistics
				nbEntries += 1
				nbEvents += len(entry.split()) -1
	sampleFile.close()
		
	if verbose:
		print "     ", nbEntries, "entries, ", nbEvents, "events."
		print "   >> end:browseEvents\n"
	return nbEntries, nbEvents

##
##
def bybloStats(sampleFileNames, outputDir, bybloDir, bybloParams="", reuse=[], verbose=False):
	print "\n>> start:bybloStats"
	
	## prepare for Byblo output
	thesauriDir = join(outputDir, "thesauri")
	if not exists(thesauriDir):
		os.makedirs(thesauriDir)
	statsDir = join(outputDir, "stats")
	statsFileNames = []
	
	#~ print "statsDir = "+statsDir
	#~ print "thesauriDir = "+thesauriDir
	
	for fileName in sampleFileNames:
		
		## determine and store statistics file name for output
		statsFileName = join(statsDir, basename(fileName) + ".stats.byblo")
		statsFileNames.append(statsFileName)
		
		if "byblo_stats" in reuse and isfile(statsFileName):
			if verbose:
				print "   Reusing statistics file " + basename(statsFileName)
		else:
			if verbose:
				print "   Creating statistics file for " + basename(statsFileName)
		
			## run Byblo for this sample file
			runTime = runByblo(abspath(fileName), abspath(thesauriDir), \
				abspath(bybloDir), bybloParams, verbose)
			
			## write statistics
			resultFileName = join(thesauriDir, basename(fileName))
			statsFile = open(statsFileName, 'w')
			
			for suffix in [".sims.neighbours",".events.filtered",".entries.filtered"]:
				size = getsize(resultFileName + suffix)
				nbLines = sum(1 for line in open(resultFileName + suffix, 'r'))
				statsFile.write(str(size) + " Size_In_Bytes_Of_File_" + suffix + '\n' +\
					str(nbLines) + " Number_Of_Lines_In_File_" + suffix + '\n')
			statsFile.write(str(runTime) + " Byblo_Run_Time\n")
		
	print ">> end:bybloStats"
	return statsFileNames

##
##
def runByblo(inputFileName, outputDir,  bybloDir, bybloParams, verbose=False):
	if verbose:
		print "\n   >> start:singleBybloRun "
		#~ print  "Input =", inputFileName
		#~ print  "Output=", outputDir
	
	## move to BYblo directory, start timer and run it in a subprocess
	startDir=abspath(os.getcwd())
	logFile = open(os.devnull, 'w') if not verbose else None
	os.chdir(bybloDir)
	if verbose:
		print "  Moved to " + os.getcwd()
	
	stime = datetime.datetime.now()
	out = subprocess.call(abspath("./byblo.sh ") + " -i " + inputFileName + " -o " + outputDir +\
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
	
	runTime = 1.0*(etime - stime).seconds / 3600
	if verbose:
		print "   Run time =", runTime
		print "   >> end:singleBybloRun \n"
	return runTime


##
##
def generateStringsFiles(sampleFileNames, outputDir, bybloDir, reuse=[], verbose=False):
	print "\n>> start:generateStringsFiles"
	
	thesauriDir = abspath(join(outputDir, "thesauri"))
	#~ print  "Thesauri dir =", thesauriDir			
	#~ print  "Byblo dir =", bybloDir
		
	for fileName in sampleFileNames:
		for typeSuffix in [".entries", ".features", ".events", ".sims"]:
			inputFileName = abspath(join(thesauriDir, basename(fileName)))
			
			## move to Byblo directory
			startDir=abspath(os.getcwd())
			os.chdir(bybloDir)
			if verbose:
				print "\n   Moved to " + os.getcwd()
			
			## convert both filtered and unfiltered versions
			for filterSuffix in (["", ".filtered"] if typeSuffix != ".sims" else [""]):
				sourceFileName = inputFileName+typeSuffix+filterSuffix
				
				## ...but only if needed
				if "byblo_stats" in reuse and isfile(sourceFileName+".strings"):
					if verbose:
						print "   Reusing strings output file " + basename(sourceFileName)
				else:
					if verbose:
						print "   Creating strings output file " + basename(sourceFileName)
				
					logFile = open(os.devnull, 'w') if not verbose else None
					
					out = subprocess.call(abspath("./tools.sh") + " unindex-" + typeSuffix[1:] \
						+ " -i " + sourceFileName\
						+ " -o " + sourceFileName + ".strings" \
						+ (" -Xe "+inputFileName+".entry-index" if typeSuffix != ".features" else "") \
						+ (" -Xf "+inputFileName+".feature-index" if typeSuffix != ".entries" else "") \
						+ " -et JDBC",\
						shell = True, stdout = logFile, stderr = logFile)
					
					if  logFile != None:
						logFile.close()
			
			## move back to execution directory
			os.chdir(startDir)
			if verbose:
				print "   Moved back to " + os.getcwd()
	print ">> end:generateStringsFiles"

##
##
def generateHistograms(sampleFileNames, outputDir, reuse=[], verbose=False):
	print "\n>> start:generateHistograms "
	
	thesauriDir = join(outputDir, "thesauri")
	graphsDir = join(outputDir, "graphs")
	if not exists(graphsDir):
		os.makedirs(graphsDir)
	#~ print "thesauriDir = "+thesauriDir
	#~ print "graphsDir = "+graphsDir
	
	for fileName in sampleFileNames:
		for suffix in ['.entries','.events', '.sims']:
			fileBaseName = basename(fileName) + suffix
			if "graphs" in reuse and isfile(join(graphsDir, "Histogram-" + fileBaseName + ".pdf")):
				if verbose:
					print "   Reusing histogram for " + fileBaseName
			else:
				if verbose:
					print "   Creating histogram for " + fileBaseName + "\n"
				
				## determine label and type
				label = "Similarity" if suffix == '.sims' else "Occurence"
				
				## create histogram with a custom step value when required
				if suffix == '.entries' or suffix == '.events':
					createHistogram(label, fileBaseName, thesauriDir, graphsDir, 1, verbose)
				else:
					createHistogram(label, fileBaseName, thesauriDir, graphsDir, 0.01, verbose)
	
	print ">> end:generateHistograms"
	
##
##
def createHistogram(label, fileName, thesauriDir, graphsDir, step=1, verbose=False):
	
	## generate the bins
	reducedFileSuffix = ".filtered" if label == "Occurence" else ".neighbours"
	bins = extractRowsValues(join(thesauriDir, fileName), step, verbose)
	reducedBins = extractRowsValues(join(thesauriDir, fileName + reducedFileSuffix), step, verbose)
	
	## prepare for graphs
	pl.suptitle(label+' histogram for file ' + fileName[string.rfind(fileName, '.'):] + ' files', fontsize=14)
	if label == "Occurence":
		pl.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.5)
		pl.subplot(211)
		pl.title('Normal scale', fontsize = 12)
	
	## BAR CHART
	## representing the bins
	left = [i*step for i in xrange(len(bins))]
	noThreshold = pl.bar(left, bins, width=step, color = "orange")
	withThreshold = pl.bar(left, reducedBins, width=step, color = "red")
	
	pl.legend((noThreshold[0], withThreshold[0]), ("No threshold", "With threshold"))
	pl.xlabel(label+" values", fontsize=10)
	pl.ylabel("Frequency", fontsize=10)
	
	## LOG-LOG LINE
	## when dealing with numbers of occurences, add a log scale representation
	if label == "Occurence":
		pl.subplot(212)
		pl.title('Log scale', fontsize = 12)
		noThreshold = pl.loglog(left, bins, color = "orange")
		withThreshold = pl.loglog(left, reducedBins, color = "red")
		
		pl.legend((noThreshold[0], withThreshold[0]), ("No threshold", "With threshold"))
		pl.xlabel(label+" values", fontsize=10)
		pl.ylabel("Frequency", fontsize=10)
		
	pl.savefig(join(graphsDir, 'Histogram-' + fileName + '.pdf'))
	pl.close()
	
	if label == "Occurence":
		## define functions 
		powerlaw = lambda x, amp, index: amp * (x**index)
		zipf = lambda x, a: x**(-a)/special.zetac(a)
		
		## adjust data
		removeZeroValues = lambda L, L2: [L[i] for i in xrange(len(L)) if L[i]>0 and L2[i]>0]
		x = removeZeroValues(left, bins)
		y = removeZeroValues(bins,left)
		y = [1.*val/max(y) for val in y]
		
		## figure set up
		f, (normScale, logScale) = pl.subplots(2, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.35) ## add margins
		f.suptitle('Occurence histogram fitting', fontsize=14, fontweight='bold')
		labelXPos, labelYPos = -0.1, 0.5 
		
		#!>>MODIFY HERE THE METHODS TO USE<<!#
		methods = [0, 4]
		
		## NORMAL SCALE
		normScale.set_title('Normal scale', fontsize = 12, fontstyle = 'italic')
		normScale.set_ylabel("frequency")
		normScale.yaxis.set_label_coords(labelXPos, labelYPos)
		normScale.set_xlabel("number of occurences")
		
		## LOG SCALE
		logScale.set_title('Log scale', fontsize = 12, fontstyle = 'italic')
		logScale.set_ylabel("frequency")
		logScale.yaxis.set_label_coords(labelXPos, labelYPos)
		logScale.set_xlabel("number of occurences")
		
		## data
		normScale.fill_between(x, y, 0, color="lightgrey")
		logScale.fill_between(x, 1e-50, y, color="lightgrey", closed=False)
		
		## fit functions
		for m in methods:
			mFit, mLabel, mColor= zipfFittingMethod(x, y, m)
			normScale.plot(x, mFit, label=mLabel, color=mColor)
			logScale.loglog(x, mFit, label=mLabel, color=mColor)
			
		normScale.legend(loc="upper right", prop={'size':8})
		logScale.legend(loc="upper right", prop={'size':8})
		pl.xlim(xmin=min(x), xmax=max(x))
		pl.ylim(ymin=min(y), ymax=max(y))
		
		pl.savefig(join(graphsDir, 'Histogram-fit-' + fileName + '.pdf'))
		pl.close()
		print ""
	
##
##
def extractRowsValues(fileName, step=1, verbose=False):
	if verbose:
		print "   >> start:extractRowsValues from " + basename(fileName)
	
	## initialisation
	bins = [0]
	lastBin, nbLine = 0, 0
	name = ""
	
	lineOffset = 0 if fileName.endswith('.entries') or fileName.endswith('.entries.filtered') else 1
	
	## go through the file to increment the counts
	for line in open(fileName+'.strings','r'):
		nbLine += 1
		for i, field in enumerate(line.split()[lineOffset:]):
			## name: store to verify later
			if i%2 == 0:
				name = field
			## value: unless for ___FILTERED___, determine, create (if needed) and increment bin
			elif name != "___FILTERED___":
				value = float(field)
				binVal = int(value / step)
				if lastBin < binVal:
					## add empty bins up to needed index
					bins, lastBin = bins + [0 for b in xrange(lastBin, binVal)], binVal
				bins[binVal]  += 1
	
	if verbose:
		print "      Size bins = " +  str(len(bins)) + ", total values = " + str(sum(bins)) + "."
		print "   >> end:extractRowsValues\n"
	return bins


##
##
def zipfFittingMethod(xdata, ydata, method=0):
	"""
	Code from Cookbook / FittingData
	@ scipy.org
	"""
	## default return values
	fit, label, color = [], "", ""
	colors = ["red", "green", "blue", "purple", "orange", "magenta", "cyan"]
	
	## define functions
	powerlaw = lambda x, amp, index: amp * (x**index)
	zipf = lambda x, a: x**(-a)/special.zetac(a)
	diffErr = lambda p, x, y, f: (y - f(p, x))
	
	## fitting often best done by first converting to a linear equation and then fitting to a straight line:
	##  y = a * x^b   <=>   log(y) = log(a) + b*log(x)
	affine = lambda p, x: p[0] + p[1] * x
	
	## curve_fit on powerlaw
	if method == 0:
		label = "   curve_fit on powerlaw"
		pfinal, covar = optimize.curve_fit(powerlaw, xdata, ydata)
		amp, index = pfinal[0], pfinal[1]
		print label + " >> amp, index =", amp, index
		fit = powerlaw(xdata, amp, index)
		
	## curve_fit on zipf
	elif method == 1:
		label = "   curve_fit on zipf"
		pfinal, covar = optimize.curve_fit(zipf, xdata, ydata)
		s = pfinal[0]
		print label + " >> s =", s
		fit = zipf(xdata, s)
	
	## polyfit on log with powerlaw
	elif method == 2:
		label = "   polyfit on log with powerlaw"
		(b, a) = polyfit(log10(xdata), log10(ydata), 1)
		a = 10. ** a
		print label + " >> a, b =", a, b
		fit = powerlaw(xdata, a, b)
	
	## polyfit on log with zipf
	elif method == 3:
		label = "   polyfit on log with zipf"
		(b, a) = polyfit(log10(xdata), log10(ydata), 1)
		print label + " >> s =", b
		fit = zipd(xdata, b)
	
	## leastsq on affine log with powerlaw
	elif method == 4:
		label = "   leastsq on affine log with powerlaw"
		logx, logy = log10(xdata), log10(ydata)
		pinit = [1.0, -1.0]
		out = optimize.leastsq(diffErr, pinit, args=(logx, logy, affine), full_output=1)
		amp, index = 10.0**out[0][0], out[0][1]
		print label + " >> amp, index =", amp, index
		fit = powerlaw(xdata, amp, index)
	## polyfit + regression on affine log with powerlaw
	elif method == 5:
		label = "   polyfit + regression on affine log with powerlaw"
		logx, logy = log10(xdata), log10(ydata)
		(b, a) = polyfit(logx, logy, 1)
		a = 10. ** a
		print label + " >> a, b =", a, b
		
		(b, a, r, tt, stderr)=stats.linregress(logx, logy)
		print label + " after regression >> a, b =", a, b
		fit = powerlaw(xdata, amp, index)

	return fit, label, colors[method]
	
##
##
def statsToDictionary(fileNames):
	dictList = []
	## create dictionary for stats in each file
	for f in sorted(fileNames):
		d = {}
		for line in open(f, 'r'):
			fields = line.split()
			d[fields[1]] = float(fields[0])
		dictList.append(d)
		
	## for each key, list in tuple all values (in all dictionaries)  for this key
	finalDict = {}
	for dict in dictList:
		for k in dict.iterkeys():
			finalDict[k] = (tuple(d[k] for d in dictList if k in d))
	return finalDict
	
##
##
def createPlotInputVariationForFiles(statsDictionary, outputDirectory, originalInputFile=""):
	if verbose:
		print "   >> start:createPlotInputVariationForFiles"
	
	## figure set up
	f, (entries, lines, sizes) = pl.subplots(3, 1)
	f.set_size_inches(8.3, 11.7) ## set figure size to A4
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.35) ## add margins
	f.suptitle('Variation of the input', fontsize=14, fontweight='bold')
	sizesInputFile, sizeUnit = convertFileSize(statsDictionary["Size_In_Bytes_Of_Input_File"])
	labelXPos, labelYPos = -0.1, 0.5 
	
	## ENTRIES and EVENTS
	numberEntries = statsDictionary["Total_Number_Of_Distinct_Entries"]
	plot1 = entries.plot(sizesInputFile, numberEntries, color='royalblue', linestyle="solid", marker='o', label="entries", alpha=0.9)
	entries.set_title('Entries and events', fontsize = 12, fontstyle = 'italic')
	entries.set_ylabel("number of entries")
	entries.yaxis.set_label_coords(labelXPos, labelYPos)
	entries.set_xlabel("size ("+sizeUnit+")")
	
	events = entries.twinx()
	numberEvents = statsDictionary["Total_Number_Of_Distinct_Events"]
	plot2 = events.plot(sizesInputFile, numberEvents, color='red', linestyle="solid", marker='o', label="events", alpha=0.9)
	events.legend((plot1[0], plot2[0]), ('entries', 'events'), loc="upper left", prop={'size':8})
	events.set_ylabel("number of events")
	
	nbEventsByEntry = statsDictionary["Average_Number_Of_Events_By_Entry"]
	avgVal = sum(nbEventsByEntry) / len(nbEventsByEntry)
	minDiff, maxDiff = 1.*min(nbEventsByEntry)/avgVal -1, 1.*max(nbEventsByEntry)/avgVal -1
	infoString = "Average # of events by entry =\n" + "%.5f" % (avgVal) \
		+ " [" + "%.2f" % (minDiff*100) + "%, +" + "%.2f" % (maxDiff*100) + "%]"
	entries.text(0.75, 0.25, infoString, horizontalalignment='center', verticalalignment='center', \
		transform = entries.transAxes, fontsize=10, color="dimGrey")
	
	## draw lines and size plots for each result file
	suffixes = ['.entries.filtered', '.events.filtered', '.sims.neighbours']
	colors1 = determineColors(suffixes, statsDictionary, "Size_In_Bytes_Of_File_")
	colors2 = determineColors(suffixes, statsDictionary, "Number_Of_Lines_In_File_")
	for suffix, color1, color2 in zip(suffixes, colors1, colors2):
		linesResultFile = statsDictionary["Number_Of_Lines_In_File_"+suffix]
		lines.plot(numberEvents, linesResultFile, color=color2, linestyle="dashed", marker='o', label=suffix, alpha=0.9)
		sizesResultFile = convertFileSize(statsDictionary["Size_In_Bytes_Of_File_"+suffix], sizeUnit)
		sizes.plot(numberEvents, sizesResultFile, color=color1, linestyle="solid", marker='o', label=suffix, alpha=0.9)

	## LINES in the result files
	lines.set_title('Lines in the result files', fontsize = 12, fontstyle = 'italic')
	lines.legend(loc="upper left", prop={'size':8})
	lines.set_ylabel("number of lines")
	lines.yaxis.set_label_coords(labelXPos, labelYPos)
	lines.set_xlabel("number of distinct observed events")
	
	## SIZE of the result files
	sizes.set_title('Size of the result files', fontsize = 12, fontstyle = 'italic')
	sizes.legend(loc="upper left", prop={'size':8})
	sizes.set_ylabel("size ("+sizeUnit+")")
	sizes.yaxis.set_label_coords(labelXPos, labelYPos)
	sizes.set_xlabel("number of distinct observed events")
		
	## save figure
	f.savefig(join(outputDir, "graphs", 'Input-variation-files-'+basename(originalInputFile)+'.pdf'))
	pl.close()
	if verbose:
		print "   >> end:createPlotInputVariationForFiles\n"

##
##
def createPlotInputVariationForTime(statsDictionary, outputDirectory, originalInputFile=""):
	if verbose:
		print "   >> start:createPlotInputVariationForTime"
	
	## figure set up
	f, (runTime) = pl.subplots(1, 1, sharex=True)
	f.set_size_inches(8.3, 11.7) ## set figure size to A4
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=None, bottom=0.5) ## add margins
	f.suptitle('Variation of the input file size\nImpact on run time', fontsize=14, fontweight='bold')
	
	numberEvents = statsDictionary["Total_Number_Of_Distinct_Events"]
	labelXPos, labelYPos = -0.1, 0.5 
	
	## BYBLO RUN TIME
	times = statsDictionary["Byblo_Run_Time"]
	print times
	runTime.fill_between(numberEvents, 0, times, color='aquamarine', label="runTime", alpha=0.9)
	runTime.set_xlabel("number of distinct observed events")
	runTime.set_ylabel("time (hours)")
	runTime.yaxis.set_label_coords(labelXPos, labelYPos)
	pl.xlim(xmin=min(numberEvents), xmax=max(numberEvents))
	pl.ylim(ymin=0)
		
	## save figure
	f.savefig(join(outputDir, "graphs", 'Input-variation-time-'+basename(originalInputFile)+'.pdf'))
	pl.close()
	if verbose:
		print "   >> end:createPlotInputVariationForTime\n"
		
## Converts a list of file sizes to a target unit (bytes system), either specified or determined based on the average size
## @return list of converted sixes, with the determined unit if it wasn't specified
def convertFileSize(sizesList, fixedUnit=None):
	units = ['bytes', 'kB', 'MB', 'GB', 'TB']
	
	## find and convert to most appropriate unit
	if fixedUnit == None:
		avg = lambda L: sum(L) / len(L)
		avgSize = avg(sizesList)
		
		for i, unit in enumerate(units):
			if i == len(units)-1 or avgSize < 1000**(i+1):
				return sizesList if i==0 else [s / 1000**i for s in sizesList], unit
	## convert to specified unit
	else:
		i = units.index(fixedUnit)
		return sizesList if i==0 else [s / 1000**i for s in sizesList]    

## Merges colors used in a graph for sets that are identical (so that they can be seen as such)
##@return list of colors to use
## //TODO: variable list size, map color additions
def determineColors(suffixes, dict, stat):
	toUse = ['red', 'RoyalBlue', 'gold']
	avail = toUse + ['MediumVioletRed', 'DarkOrange', 'Chartreuse']
	for i1, s1 in enumerate(suffixes):
		for i2, s2 in enumerate(suffixes):
			if i2 > i1 and dict[stat+s1] == dict[stat+s2]:
				## base colors => merge
				if avail.index(toUse[i1]) < 3 and avail.index(toUse[i2]) < 3:
					toUse[i1] = toUse[i2] = avail[i1+i2+2]
				## at least one "complex" color => everything to grey
				else:
					toUse = ['grey']*3
					return toUse
	return toUse
	
	
## Deletes intermediary files needed for the statistics generation when they are in the list of item types (samples, 
## thesauri, stats) that the user specified for deletion
def deleteOnExit(deleteList, outputDirectory, sampleFiles=[], originalInputFile=""):
	print "\n>> start:deleteOnExit "
	if(deleteList != []):
		for directory in deleteList:
			## sample files
			if directory == "samples":
				for file in sampleFiles:
					if file != originalInputFile:
						os.system("rm -r " + file)
						if verbose:
							print "   Deleted file " + file
			## directories
			else:
				os.system("rm -r " + join(outputDirectory, directory))
				if verbose:
					print "   Deleted directory " + directory
	print ">> end:deleteOnExit"


##
##
def print_lines(list, min=0, max=0, line_max=0, title="List"):
	if max == 0:
		max = len(list)
	print "\n" + title
	if type(list) is dict:
		for i, key in enumerate(sorted(list.iterkeys())) :
			if i<max:
				print "   " + (key + str(list[key]))[:line_max] +  "..." \
					if line_max != 0 and len(key + str(list[key])) > line_max else "   " + key + str(list[key])
			else:
				break
	else:
		for index in range(min, max):
			print "   " + str(list[index])[:line_max] + "..." \
				if  line_max != 0 and len(str(list[index])) > line_max else "   " + str(list[index])
	print ""


if __name__=='__main__':
	
	##############################################################################
	parser = argparse.ArgumentParser(description='Generate statistics for Byblo.')

	# data file name
	parser.add_argument('data', metavar='file', nargs=1, 
		action='store', 
		help='data file containing  a feature / entry set')
	# percentages of entries to take into account
	parser.add_argument('-p', '--percentage', metavar='x', type=float, nargs='*', dest='p', 
		action='store', default=None,
		help='approximate percentage (]0.-100.]) of entries to keep  for a set of statistics')
	# output directory
	parser.add_argument('-O', '--output-dir', metavar='path', nargs=1, dest='O', 
		action='store', default=None,
		help='output directory where sub-directories "stats", "thesauri" and "graphs" can be found'+\
		'(default: ".")')
	# Byblo directory
	parser.add_argument('-B', '--Byblo-dir', metavar='path', nargs=1, dest='B', 
		action='store', default=None,
		help='Byblo directory, usually named "Byblo-x.x.x" '+\
		'(default: "../Byblo-2.0.1")')
	# BYblo parameters
	parser.add_argument('-P', '--Byblo-params', metavar='string', nargs=1, dest='P', 
		action='store', default=None,
		help='Byblo parameters, read Byblo doc for more info '+\
		'(default: "-fef 10 -fff 10 -t 6 -Smn 0.1")')
	# delete option (for files and directories created in the proces of the graph generation)
	parser.add_argument('-d', '--delete', metavar='string', nargs='*', dest='delete',
		action='store', choices=('samples', 'thesauri', 'stats'),  
		help='delete specified items from "samples", "thesauri" and/or "stats" before exiting, '+\
		' everything when none specifed (default: delete nothing)')
	# reuse option (for files and directories created in the proces of the graph generation)
	parser.add_argument('-r', '--reuse', metavar='string', nargs='*', dest='reuse',
		action='store', choices=('samples', 'events_stats', 'byblo_stats', 'graphs'),  
		help='reuse specified items from "samples", "events_stats", "byblo_stats" and/or "graphs" during execution, '+\
		' everything when none specifed (default: reuse nothing)')
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_const', 
		const=True, default=False, 
		help='display Byblo output and statistics information for each step (default: False)')
	
	args = parser.parse_args()
	##############################################################################
	
	pctList = [p for p in args.p if p != 0 and p<=100] if args.p != None else [100]
	outputDir = args.O[0] if args.O != None else "."
	bybloDir = args.B[0] if args.B !=None else "../Byblo-2.0.1"
	bybloParams = args.P[0] if args.P != None else "-fef 10 -fff 10 -t 6 -Smn 0.1"
	delList = (set(args.delete) if args.delete != [] else ['samples', 'thesauri', 'stats']) \
		if args.delete != None else []
	reuseList = (set(args.reuse) if args.reuse != [] else ['samples', 'events_stats', 'byblo_stats', 'graphs']) \
		if args.reuse != None else []
	verbose = args.verbose
	
	stime = datetime.datetime.now()
	print "***************************************************************************"
	print "BYBLO STATISTICS TOOL"
	print "***************************************************************************\n"
	
	## EVENTS STATS + sample file creation when necessary
	sampleFileNames, statsFileNames = eventsStats(args.data[0], outputDir, pctList, reuseList, verbose)
	if verbose:
		print_lines(statsFileNames, title="Statistics files after feature extraction")
	
	## BYBLO STATS
	statsFileNames += bybloStats(sampleFileNames, outputDir, bybloDir, bybloParams, reuseList, verbose)
	if verbose:
		print_lines(statsFileNames, title="Statistics files after byblo run")
	
	## FILE CONVERSION [RESTORE STRINGS]
	generateStringsFiles(sampleFileNames, outputDir, bybloDir, reuseList, verbose)
	
	
	## HISTOGRAMS
	generateHistograms(sampleFileNames, outputDir, reuseList, verbose)
	"""
	## PLOTS
	statsDict = statsToDictionary(statsFileNames)
	if verbose:
		print_lines(statsDict, title="Table containing the statistics")
	
	createPlotInputVariationForFiles(statsDict, outputDir, args.data[0])
	createPlotInputVariationForTime(statsDict, outputDir, args.data[0])
	"""
	## CLEAN UP
	deleteOnExit(delList, outputDir, sampleFileNames, args.data[0])
	
	etime = datetime.datetime.now()
	print "\n>Execution took", etime-stime, "hours"

	    
