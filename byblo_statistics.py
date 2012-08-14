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
from scipy import optimize, stats, special

#################################################################

## Determines the set of sample files to generate statistics for using an original input file and a list of percentages,
## creates these files in the subfolder "samples" when required, then goes through all of them and writes some of 
## their carateristics (size, entries, events) in statistics files put in the subfolder "stats"
## @return list of sample files, list of statistics files
def eventsStats(inputFileName, outputDir, percentList=[100], reuse=[], verbose=False):
	print "\n>> start:eventsStats - full file size=", getsize(inputFileName), "bytes"
	
	if not exists(outputDir):
		os.makedirs(outputDir)
	
	statsDir = join(outputDir, "stats")
	if not exists(statsDir):
		os.makedirs(statsDir)
		
	samplesDir = join(outputDir, "samples")
	if not exists(samplesDir):
		os.makedirs(samplesDir)
	
	sampleFileNames, statsFileNames = [], []
	for pct in percentList:
		
		## determine and store data file name for input
		sampleFileName = join(samplesDir, basename(inputFileName) + "-%" + string.replace(str(pct), '.', '-')) \
			if pct != 100 else inputFileName
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


## Runs Byblo for each sample file / parameter string association and writes information about its behaviour (size
## and number of lines of result files, total run time) in new statistics files put in the subfolder "stats"
## @return list of statistics files
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


## Runs Byblo for the specified input file with all of the parameter strings, putting the result files (indexes, counts,
## thesauri) in the subfolder "thesauri" and recording the running times
## @return list of running times
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
	
	## fail?
	if(not out == 0):
		print "   Byblo failed on input file: " + inputFileName + "\n   Fail Code: " + str(out)
		sys.exit()
		
	if  logFile != None:
		logFile.close()
	os.chdir(startDir)
	if verbose:
		print "  Moved back to " + os.getcwd()
	
	runTime = 1.0*(etime - stime).seconds
	if verbose:
		print "   Run time =", runTime
		print "   >> end:singleBybloRun \n"
	return runTime


## Converts Byblo result files that use skip indexing so that the strings represented by the indexes are restored 
## (slower and heavier but better for result analysis, readability and adaptability)
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


## Makes a graph nicer and clearer by adding a title, axes labels, a legend, limits on the axes (for an exact fit or a 
## "zoom"), a scientific notation (for large numbers)
def decorateGraph(subplot, title="", xLabel="", yLabel="", yLabelPos=None, legendPos=None, \
	data=None, largeX=False, largeY=False):
	## text: title, legend and axes labels
	if title != "":
		subplot.set_title(title, fontsize = 12, fontstyle = 'italic')
	if legendPos != None:
		subplot.legend(loc=legendPos, prop={'size':8})
	if xLabel != "":
		subplot.set_xlabel(xLabel, fontsize=10)
	if yLabel != "":
		if yLabelPos != None:
			subplot.yaxis.set_label_coords(yLabelPos[0], yLabelPos[1])
		subplot.set_ylabel(yLabel, fontsize=10)
	## axes limits (e.g. log scale or histograms)
	if data != None:
		subplot.set_xlim(xmin=min(data[0]), xmax=max(data[0]))
		subplot.set_ylim(ymin=min(data[1]))
	## scientific notation for large values
	if largeX:
		subplot.ticklabel_format(style='sci', scilimits=(0,0), axis='x')
	if largeY:
		subplot.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
	
	
## Offers a selection of fitting methods and model functions that allow to model some of the relations that exist within
## the data handled here
## @return model data produced with the model function, method name, associated colour (for graphs)
def fittingMethod(xdata, ydata, method=0):
	## default return values
	fit, label, color = [], "", ""
	colors = ["red", "green", "blue", "purple", "orange", "magenta", "cyan", "limegreen", "gold"]
	
	## fitting functions
	powerlaw = lambda x, amp, index: amp * (x**index)
	zipf = lambda x, a: x**(-a)/special.zetac(a)
	diffErr = lambda p, x, y, f: (y - f(p, x))
	quadratic = lambda p, x: p[0] * (x**np.array(2))
	#~ nlogn = lambda  p, x: p[1] * p[0] * x * np.log(p[0] * x)
	nlogn = lambda  a, b, x: b * a * x * np.log(a * x)
	
	## fitting often best done by first converting to a linear equation and then fitting to a straight line:
	##  y = a * x^b   <=>   log(y) = log(a) + b*log(x)
	affine = lambda p, x: p[0] + p[1] * x
	
	## conversion functions (needed when handling probabilities)
	toProbabilities = lambda data: [1.*val/sum(data) for val in data]
	toFrequencies = lambda data, originalData: [1.*val*sum(originalData) for val in data] 
	
	## curve_fit on powerlaw
	if method == 0:
		label = "curve_fit on powerlaw"
		pfinal, covar = optimize.curve_fit(powerlaw, xdata, ydata)
		amp, index = pfinal[0], pfinal[1]
		print "   ["+str(method)+"] " + label + " >> amp, index =", amp, index
		fit = powerlaw(xdata, amp, index)
		
	## curve_fit on zipf
	elif method == 1:
		y = toProbabilities(ydata) ## convert to probabilities
		label = "curve_fit on zipf"
		pfinal, covar = optimize.curve_fit(zipf, xdata, y)
		s = pfinal[0]
		print "   ["+str(method)+"] " + label + " >> s =", s
		fit = zipf(xdata, s)
		fit = toFrequencies(fit, ydata) ## restore frequencies
	
	## polyfit on log with powerlaw
	elif method == 2:
		label = "polyfit on log with powerlaw"
		(a, b) = polyfit(log10(xdata), log10(ydata), 1)
		b = 10. ** b
		print "   ["+str(method)+"] " + label + " >> a, b =", a, b
		fit = powerlaw(xdata, b, a)
	
	## polyfit on log with zipf
	elif method == 3:
		y = toProbabilities(ydata) ## convert to probabilities
		label = "polyfit on log with zipf"
		(a, b) = polyfit(log10(xdata), log10(y), 1)
		print "   ["+str(method)+"] " + label + " >> s =", -a
		fit = zipf(xdata, -a)
		fit = toFrequencies(fit, ydata) ## restore frequencies
	
	## leastsq on log with affine + powerlaw
	elif method == 4:
		label = "leastsq on log with affine + powerlaw"
		logx, logy = log10(xdata), log10(ydata)
		pinit = [1.0, -1.0]
		out = optimize.leastsq(diffErr, pinit, args=(logx, logy, affine), full_output=1)
		amp, index = 10.0**out[0][0], out[0][1]
		print "   ["+str(method)+"] " + label + " >> amp, index =", amp, index
		fit = powerlaw(xdata, amp, index)
		
	## leastsq on log with affine + zipf
	elif method == 5:
		y = toProbabilities(ydata) ## convert to probabilities
		label = "leastsq on log with affine + zipf"
		pinit = [1.0, -1.0]
		out = optimize.leastsq(diffErr, pinit, args=(log10(xdata), log10(y), affine), full_output=1)
		a =  out[0][1]
		print "   ["+str(method)+"] " + label + " >> s =", -a
		fit = zipf(xdata, -a)
		fit = toFrequencies(fit, ydata) ## restore frequencies
	
	## polyfit / polyval - 2nd degree
	elif method == 6:
		label = "polyfit / polyval - 2nd degree"
		params = polyfit(xdata, ydata, 2)
		print "   ["+str(method)+"] " + label + " >> parameters =", params
		fit = polyval(params, xdata)
		return fit, label, params
		
	## leastsq on quadratic
	elif method == 7:
		label = "leastsq on quadratic"
		pinit = [1.0]
		out = optimize.leastsq(diffErr, pinit, args=(xdata, ydata, quadratic), full_output=1)
		print "   ["+str(method)+"] " + label + " >> parameters =", out[0]
		fit = quadratic(out[0], xdata)
		return fit, label, out[0]
	
	## curve_fit on nlogn
	elif method == 8:
		label = "curve_fit on nlogn"
		pfinal, covar = optimize.curve_fit(nlogn, xdata, ydata)
		print "   ["+str(method)+"] " + label + " >> parameters =", pfinal
		fit = nlogn(pfinal[0], pfinal[1], xdata)
		return fit, label, pfinal

	return fit, label, colors[method]


## Generates histograms for Byblo result files (all sizes and parameter strings) containing counts of entries, features 
## and events, as well as thesauri containing similarity values, always combining filtered and unfiltered versions of the
## same file
def generateHistograms(sampleFileNames, outputDir, reuse=[], verbose=False):
	print "\n>> start:generateHistograms "
	
	thesauriDir = join(outputDir, "thesauri")
	graphsDir = join(outputDir, "graphs")
	if not exists(graphsDir):
		os.makedirs(graphsDir)
	#~ print "thesauriDir = "+thesauriDir
	#~ print "graphsDir = "+graphsDir
	
	for fileName in sampleFileNames:
		for suffix in ['.entries', '.features', '.events', '.sims']:
		#~ for suffix in ['.entries', '.events', '.sims']:
			fileBaseName = basename(fileName) + suffix
			if "graphs" in reuse and isfile(join(graphsDir, "Histogram-" + fileBaseName + ".pdf")):
				if verbose:
					print "   Reusing histogram for " + fileBaseName
			else:
				if verbose:
					print "   Creating histogram for " + fileBaseName + "\n"
				
				## create histogram with a custom step value when required
				if suffix == '.sims':
					createSimilarityHistogram("Similarity", fileBaseName, thesauriDir, graphsDir, verbose)
				else:
					createOccurenceHistogram("Occurence", fileBaseName, thesauriDir, graphsDir, verbose)
	
	print ">> end:generateHistograms"


## Creates a histogram showing the distribution of frequencies for the chosen element (entry, feature, event)
## Both normal and logarithmic scale are created, together with a (for now (very) false) model attempt
def createOccurenceHistogram(label, fileName, thesauriDir, graphsDir, verbose=False):
	## generate the histograms
	XBASE, YBASE = 2, 10
	LIMITS = [1, 2.0 ** 64]
	reducedFileSuffix = ".filtered"
	bins, hist = extractRowsValues(join(thesauriDir, fileName), LIMITS, [XBASE, YBASE], verbose=verbose)
	reducedBins, reducedHist = extractRowsValues(join(thesauriDir, fileName + reducedFileSuffix), LIMITS, [XBASE, YBASE], verbose=verbose)
	
	## figure set up
	f, (normScale, logScale, fitNormScale, fitLogScale) = pl.subplots(4, 1)
	f.set_size_inches(8.3, 11.7) ## set figure size to A4
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.4) ## add margins
	f.suptitle('Occurence histogram for  ' + fileName[string.rfind(fileName, '.'):] + ' file', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5]
	
	## REPRESENT THE DATA
	## linear bar chart of 99% of the mass
	noThreshold = normScale.bar(bins[:-1], hist, width=bins[1:] - bins[:-1], color="orange", label="No threshold")
	withThreshold = normScale.bar(reducedBins[:-1], reducedHist, width=reducedBins[1:] - bins[:-1], color="red", label="With threshold")
	## log-log line
	noThreshold = logScale.loglog(bins[:-1], hist, 'x-', color="orange", label="No threshold")
	withThreshold = logScale.loglog(reducedBins[:-1], reducedHist, 'x-', color="red", label="With threshold")
	
	pp = max([i for i,x in enumerate((np.cumsum(hist) / sum(hist)) < 0.99) if x])
	decorateGraph(normScale, 'Norm scale', "number of occurences", "frequency", yLabelPos, "upper right", largeY=True)
	decorateGraph(logScale, 'Log scale', "number of occurences", "frequency", yLabelPos, "upper right")
	normScale.set_xlim(xmax=XBASE ** max(pp, XBASE ** 2))
	
	## FIT THE DATA
	## data
	removeZeroValues = lambda L, L2: [L[i] for i in xrange(len(L)) if L[i]>0 and L2[i]>0]
	x, y = removeZeroValues(bins[:-1], hist), removeZeroValues(hist, bins[:-1])
	fitNormScale.fill_between(x, y, 0, color="lightgrey")
	fitLogScale.fill_between(x, 1e-50, y, color="lightgrey")
	## fit functions
	#!>>MODIFY HERE THE METHODS TO USE<<!#
	methods = [1, 3, 4]
	for m in methods:
		mFit, mLabel, mColor= fittingMethod(x, y, m)
		fitNormScale.plot(x, mFit, label=mLabel, color=mColor)
		fitLogScale.loglog(x, mFit, label=mLabel, color=mColor)
		
	decorateGraph(fitNormScale, 'Norm scale', "number of occurences", "frequency",  yLabelPos, "upper right", data=(x, y), largeY=True)
	decorateGraph(fitLogScale, 'Log scale', "number of occurences", "frequency", yLabelPos, "upper right", data=(x, y))
	fitNormScale.set_xlim(xmax=XBASE** max(pp, XBASE ** 2))
	
	f.savefig(join(graphsDir, 'Histogram-fit-' + fileName + '.pdf'))
	pl.close()
	print ""

## Creates a histogram showing the distribution of similarity scores that appear in a thesaurus
## Only normal scale (values between 0 and 1), no model
def createSimilarityHistogram(label, fileName, thesauriDir, graphsDir, verbose=False):
	## generate the histograms
	XBASE, YBASE = 10, 10
	LIMITS = [0, 1]
	step = 0.01
	reducedFileSuffix = ".neighbours"
	bins, hist = extractRowsValues(join(thesauriDir, fileName), LIMITS, [XBASE, YBASE], step, verbose)
	reducedBins, reducedHist = extractRowsValues(join(thesauriDir, fileName + reducedFileSuffix), LIMITS, [XBASE, YBASE], step, verbose)
	
	## figure set up
	f, (sims) = pl.subplots()
	f.set_size_inches(8.3, 5.8) ## set figure size to A5
	f.subplots_adjust(left=0.15, right=0.85) ## add margins
	f.suptitle('Similarity histogram for  ' + fileName[string.rfind(fileName, '.'):] + ' file', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5]
	
	## REPRESENT THE DATA
	## bar chart
	noThreshold = sims.bar(bins[:-1], hist, width=bins[1:] - bins[:-1], color="orange", label="No threshold")
	withThreshold = sims.bar(reducedBins[:-1], reducedHist, width=reducedBins[1:] - bins[:-1], color="red", label="With threshold")
	decorateGraph(sims, 'Norm scale', "similarity score", "frequency", yLabelPos, "upper right", largeY=True)

	f.savefig(join(graphsDir, 'Histogram-' + fileName + '.pdf'))
	pl.close()
	print ""
	

## Creates bins of width "step" using all of the values from the specified file that are relevant for the corresponding histogram
## @return array of bins
def extractRowsValues(fileName, limits, bases=None, step=None, verbose=False):
	if verbose:
		print "   >> start:extractRowsValues2 from " + basename(fileName)
	
	## initialisation
	name = ""
	chunk = []
	MIN_VALUE, MAX_VALUE = limits[0], limits[1]
	
	if step:
		bins = np.arange(MIN_VALUE, MAX_VALUE, step)
	else:
		# calculate the min and max powers
		start_power = np.floor(np.log(MIN_VALUE) / np.log(bases[0]))
		end_power = np.ceil(np.log(MAX_VALUE) / np.log(bases[0]))
		#  generate a range of delimiters in log space
		num_bins = (end_power - start_power) + 1
		bins = np.logspace(start_power, end_power, num_bins, base=bases[0])
		
	hist = np.zeros(len(bins) - 1)
	lineOffset = 0 if '.entries' in fileName or '.features' in fileName else 1
	
	## go through the file to increment the counts
	for line in open(fileName+'.strings','r') :
		for i, field in enumerate(line.split()[lineOffset:]):
			## name: store to verify later
			if i%2 == 0:
				name = field
			## value: unless for ___FILTERED___, determine, create (if needed) and increment bin
			elif name != "___FILTERED___":
				chunk.append(float(field))
				##  Iteratively generate the histogram in 10k chunks
				if len(chunk) == 1e5:
					(tmp,_) = np.histogram(chunk, bins=bins)
					hist += tmp
					chunk[:] = []
	
	##  don't forget the last uncomplete chunk
	if len(chunk) != 0:
		(tmp,_) = np.histogram(chunk, bins=bins)
		hist += tmp				

	## slice all the empty bins of the end
	last_idx = max([i for i,h in enumerate(hist) if h])
	hist = hist[:last_idx+1]
	bins = bins[:last_idx+2]
	
	if verbose:
		print "      Size hist = " +  str(len(hist)) + ", total values = " + str(sum(hist)) + "."
		print "   >> end:extractRowsValues2\n"
	return bins, hist
	

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
	yLabelPos = [-0.1, 0.5] 
	
	## ENTRIES and EVENTS
	numberEntries = statsDictionary["Total_Number_Of_Distinct_Entries"]
	plot1 = entries.plot(sizesInputFile, numberEntries, color='royalblue', linestyle="solid", marker='o', label="entries", alpha=0.9)
	decorateGraph(entries, 'Entries and events', "size ("+sizeUnit+")", "number of entries", yLabelPos, largeY=True)
	
	events = entries.twinx()
	numberEvents = statsDictionary["Total_Number_Of_Distinct_Events"]
	plot2 = events.plot(sizesInputFile, numberEvents, color='red', linestyle="solid", marker='o', label="events", alpha=0.9)
	decorateGraph(events, yLabel="number of events", largeY=True)
	events.legend((plot1[0], plot2[0]), ('entries', 'events'), loc="upper left", prop={'size':8})
	
	## display average number of observed events by entry (with error)
	nbEventsByEntry = statsDictionary["Average_Number_Of_Events_By_Entry"]
	avgVal = sum(nbEventsByEntry) / len(nbEventsByEntry)
	minDiff, maxDiff = 1.*min(nbEventsByEntry)/avgVal -1, 1.*max(nbEventsByEntry)/avgVal -1
	infoString = "Average # of events by entry =\n" + "%.5f" % (avgVal) \
		+ " [" + "%.2f" % (minDiff*100) + "%, +" + "%.2f" % (maxDiff*100) + "%]"
	entries.text(0.75, 0.25, infoString, horizontalalignment='center', verticalalignment='center', \
		transform = entries.transAxes, fontsize=10, color="dimGrey")
	
	## draw LINES and SIZE plots for each result file
	suffixes = ['.entries.filtered', '.events.filtered', '.sims.neighbours']
	colors1 = determineColors(suffixes, statsDictionary, "Size_In_Bytes_Of_File_")
	colors2 = determineColors(suffixes, statsDictionary, "Number_Of_Lines_In_File_")
	for suffix, color1, color2 in zip(suffixes, colors1, colors2):
		linesResultFile = statsDictionary["Number_Of_Lines_In_File_"+suffix]
		lines.plot(numberEvents, linesResultFile, color=color2, linestyle="dashed", marker='o', label=suffix, alpha=0.9)
		sizesResultFile = convertFileSize(statsDictionary["Size_In_Bytes_Of_File_"+suffix], sizeUnit)
		sizes.plot(numberEvents, sizesResultFile, color=color1, linestyle="solid", marker='o', label=suffix, alpha=0.9)
	decorateGraph(lines, 'Lines in the result files', "number of distinct observed events", "number of lines", \
		yLabelPos, "upper left", largeX=True, largeY=True)
	decorateGraph(sizes, 'Size of the result files', "number of distinct observed events", "size ("+sizeUnit+")", \
		yLabelPos, "upper left", largeX=True)
		
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
	f, (entries, events) = pl.subplots(2, 1)
	f.set_size_inches(8.3, 11.7) ## set figure size to A5
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.35) ## add margins
	f.suptitle('Variation of the input file size\nImpact on run time', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5] 
	
	## REPRESENT THE DATA
	times,  timeUnit = convertTimeRange(statsDictionary["Byblo_Run_Time"])
	numberEntries = statsDictionary["Total_Number_Of_Distinct_Entries"]
	entries.plot(numberEntries, times, color='aquamarine', label="run time function of entries", marker='o', linestyle='None')
	numberEvents = statsDictionary["Total_Number_Of_Distinct_Events"]
	events.plot(numberEvents, times, color='aquamarine', label="run time function of events", marker='o', linestyle='None')
	
	## FIT THE DATA
	quadratic = lambda p, x: p[1] * ((p[0] *x)**np.array(2))
	nlogn = lambda  p, x: p[1] * p[0] * x * np.log(p[0] * x)
	## entries [method 6: polyfit / polyval - 2nd degree]
	mFit, mLabel, mParams1= fittingMethod(numberEntries, times, 6)
	entries.plot(numberEntries, mFit, label=mLabel, color='firebrick', linestyle='dashed')
	nMin, nMax = min(numberEntries), max(numberEntries)
	steps = np.arange(nMin, nMax, (nMax-nMin)/100.) ## for a smoother model
	entries.plot(steps, polyval(mParams1, steps), label="smoothed model", color='red')
	
	## events [method 8: ...]
	#~ mFit, mLabel, mParams2= fittingMethod(numberEvents, times, 8)
	#~ events.plot(numberEvents, mFit, label=mLabel, color='firebrick', linestyle='dashed')
	#~ nMin, nMax = min(numberEvents), max(numberEvents)
	#~ steps = np.arange(nMin, nMax, (nMax-nMin)/100.) ## for a smoother model
	#~ events.plot(steps, nlogn(mParams2, steps), label="smoothed model", color='red')
	
	decorateGraph(entries, "", "number of distinct entries", "time ("+timeUnit+")", yLabelPos, "lower right", largeX=True, largeY=True)
	decorateGraph(events, "", "number of distinct observed events", "time ("+timeUnit+")", yLabelPos, "lower right", largeX=True, largeY=True)
	
	## DISPLAY APPROXIMATIONS
	
	## entries
	a, b, c = mParams1
	ratio = sum(times) / sum(numberEntries)
	infoString = "Run time approximation:\nt(n) = " \
		+ "%.3f" % (a/(ratio**2))  + " x N**2 + " +  "%.3f" % (b/ratio)  + " x N + " + "%.3f" % (c) \
		+ "\nwith N = n x " + "%.5e" % (ratio)
	entries.text(0.4, 0.85, infoString, horizontalalignment='center', verticalalignment='center', \
		transform = entries.transAxes, fontsize=10, color="dimGrey")
		
	#~ ## events
	#~ a, b = mParams2
	#~ infoString = "Run time approximation from events:" \
		#~ + "\nt(n) = " + "%.5e" % (a) + " x N log(N) + " + "%.5e" % (b) \
		#~ + "\nwith N = n x " + "%.5e" % (ratio)
	#~ events.text(0.4, 0.85, infoString, horizontalalignment='center', verticalalignment='center', \
		#~ transform = events.transAxes, fontsize=10, color="dimGrey")
	
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


## Converts a list of time ranges to a target unit either specified or determined based on the average size
## @return list of converted time ranges, with the determined unit if it wasn't specified
def convertTimeRange(timesList, fixedUnit=None):
	units = ['seconds', 'minutes', 'hours']
	
	## find and convert to most appropriate unit
	if fixedUnit == None:
		avg = lambda L: sum(L) / len(L)
		avgTime = avg(timesList)
		
		for i, unit in enumerate(units):
			if i == len(units)-1 or avgTime < 60**(i+1):
				return timesList if i==0 else [s / 60**i for s in timesList], unit
	## convert to specified unit
	else:
		i = units.index(fixedUnit)
		return timesList if i==0 else [s / 60**i for s in timesList]    


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
## thesauri, event_stats, byblo_stats) that the user specified for deletion
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
	# Byblo parameters
	parser.add_argument('-P', '--Byblo-params', metavar='string', nargs=1, dest='P', 
		action='store', default=None,
		help='Byblo parameters, read Byblo doc for more info '+\
		'(default: "-fef 10 -fff 10 -t 6 -Smn 0.1")')
	# delete option (for files and directories created in the process of the graph generation)
	parser.add_argument('-d', '--delete', metavar='string', nargs='*', dest='delete',
		action='store', choices=('samples', 'thesauri', 'stats'),  
		help='delete specified items from "samples", "thesauri" and/or "stats" before exiting, '+\
		' everything when none specifed (default: delete nothing)')
	# reuse option (for files and directories created in the process of the graph generation)
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
	
	pctList = sort([p for p in args.p if p != 0 and p<=100] if args.p != None else [100])
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
	
	## PLOTS
	statsDict = statsToDictionary(statsFileNames)
	if verbose:
		print_lines(statsDict, title="Table containing the statistics")
	
	createPlotInputVariationForFiles(statsDict, outputDir, args.data[0])
	createPlotInputVariationForTime(statsDict, outputDir, args.data[0])
	
	## CLEAN UP
	deleteOnExit(delList, outputDir, sampleFileNames, args.data[0])
	
	etime = datetime.datetime.now()
	print "\n>Execution took", etime-stime, "hours"

	    
