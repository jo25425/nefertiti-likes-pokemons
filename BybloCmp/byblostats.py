#~ !/usr/bin/env python -c

#~ Copyright (c) 2012, University of Sussex
#~ All rights reserved.

#~ Redistribution and use in source and binary forms, with or without 
#~ modification, are permitted provided that the following conditions are met:

 #~ * Redistributions of source code must retain the above copyright notice, this
   #~ list of conditions and the following disclaimer.

 #~ * Redistributions in binary form must reproduce the above copyright notice, 
   #~ this list of conditions and the following disclaimer in the documentation 
   #~ and/or other materials provided with the distribution.

 #~ * Neither the name of the University of Sussex nor the names of its 
   #~ contributors may be used to endorse or promote products  derived from this
   #~ software without specific prior written permission.

#~ THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#~ AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#~ IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#~ DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY DIRECT, 
#~ INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#~ BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#~ DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY 
#~ OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
#~ NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
#~ EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#~ -----------------------------------------------------------------------------
__author__ = "Joanne Robert"
__copyright__ = "Copyright (c) 2012, University of Sussex"
__credits__ = ["Joanne Robert", "Hamish Morgan"]
__license__ = "3-clause BSD"
__version__ = "1.0.0"
__maintainer__ = "Joanne Robert"
__email__ = "jr317@sussex.ac.uk"
__status__ = "Development"
#~ -----------------------------------------------------------------------------

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

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import *
from scipy import optimize, stats, special


## custom type for non zero percentage
def nonZeroPct(text):
	value = float(text)
	if not 0 < value <= 100:
		raise argparse.ArgumentTypeError("non zero percentage expected")
	return value

## custom type for positive integer (frequency filters)
def positiveInt(text):
	value = int(text)
	if not 0 <= value:
		raise argparse.ArgumentTypeError("positive integer expected")
	return value

## custom type for float between 0 and 1 (similarity score)
def zeroToOneFloat(text):
	value = float(text)
	if not 0 <= value <= 1:
		raise argparse.ArgumentTypeError("float between zero and one expected")
	return value
	
	
## Produces a string that gives information about size (percentage) used and can be in a file name
## @return formatted pct string
def sizeSubstring(pct):
	return '%' + string.replace(str(pct), '.', '_')

## Produces a string that gives information about parameters used and can be in a file name
## @return formatted parameters string
def paramSubstring(str):
	return '#' + string.replace(str, ' ', '_')
	
## Extracts from a file name the substrings produced by the 2 methods above
## @return string reflecting size and parameters chosen
def getPctParamSubstring(name):
	end = basename(name)
	index = string.find(end, '%')
	return "..." + end[index:]

## 
## @return complete name of the parameter
def fullParamName(initials):
	return { 'fef': 'entry frequency filter',
			'fff': 'feature frequency filter',
			'fvf': 'event frequency filter',
			'smn': 'minimum similarity '
			}[string.lower(initials)]


##
##
def generateParameterStrings(fef, fff, fvf, Smn, verbose=False):
	if verbose:
		print "\n>> start:generateParameterStrings"
	
	bybloParams = []
	if fef + fff + fvf + Smn:
		for fefStr in ['-fef '+str(n)+' ' for n in fef] if fef else ['']:
			for fffStr in ['-fff '+str(n)+' ' for n in fff] if fff else ['']:
					for fvfStr in ['-fvf '+str(n)+' ' for n in fvf] if fvf else ['']:
							for SmnStr in ['-Smn '+str(x)+' ' for x in Smn] if Smn else ['']:
								bybloParams.append(fefStr + fffStr + fvfStr + SmnStr)
	if verbose:
		print ">> end:generateParameterStrings\n"
	return bybloParams


## Displays Byblo help to give information about the possible parameter strings, then allows the user to choose
## new ones
## @return new list of parameter strings for Byblo
def displayBybloHelp(bybloDir, bybloParams, verbose=False):
	if verbose:
		print "\n>> start:displayBybloHelp"

	## move to Byblo directory, and run it in a subprocess
	startDir=abspath(os.getcwd())
	os.chdir(bybloDir)
	if verbose:
		print "   Moved to " + os.getcwd() + "\n"
	out = subprocess.call(abspath("./byblo.sh ") + " --help", shell = True)
	## fail?
	if(out != 0 and out != 255):
		print "   Byblo help failed.\n   Fail Code: " + str(out)
	os.chdir(startDir)
	if verbose:
		print "   Moved back to " + os.getcwd()
		
	## offer user to change the list of parameters
	print_lines(bybloParams, title="   The list of parameter strings to use is for now:")
	print "   Change it? (y/n)",
	x = None
	while True:
		x = raw_input()
		if x == 'y':
			## scan a new list of parameter strings
			print "   Type new list (syntax: string[, string ...]).\n  ",
			while True:
				input = raw_input()
				bybloParams = input.split(', ')
				## accepted input: non-empty list of strings
				if type(bybloParams) is list and bybloParams is not [] \
					and len(bybloParams) == len([1 for s in bybloParams if type(s) is str]):
					break
				else:
					print "   Invalid input. Syntax: string[, string ...]).\n  ",
			break
		elif x == 'n':
			break
		else:
			print "   Invalid choice. Change it? (y/n)",
	
	print_lines(bybloParams, title="   The list of parameter strings to use is now:")
	
	if verbose:
		print ">> end:displayBybloHelp\n"
	return bybloParams


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
		sampleFileName = join(samplesDir, basename(inputFileName) + sizeSubstring(pct)) \
			if pct != 100 else inputFileName
		sampleFileNames.append(sampleFileName)
		
		## determine and store statistics file name for output
		statsFileName = join(statsDir, basename(sampleFileName) + (sizeSubstring(pct) if pct == 100 else '') +  ".stats.events")
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
				print "   Writing statistics in " + getPctParamSubstring(statsFileName) + "\n"
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
def bybloStats(sampleFileNames, outputDir, bybloDir, paramList=[], reuse=[], verbose=False):
	print "\n>> start:bybloStats"
	
	## prepare for Byblo output
	thesauriDir = join(outputDir, "thesauri")
	if not exists(thesauriDir):
		os.makedirs(thesauriDir)
	statsDir = join(outputDir, "stats")
	statsFileNames = []
	
	for fileName in sampleFileNames:
		for paramStr in paramList:
			## determine and store statistics file name for output
			fullInputFileStr = sizeSubstring(100.0) if string.find(fileName, '%') == -1 else '' 
			statsFileName = join(statsDir, basename(fileName) + fullInputFileStr + paramSubstring(paramStr) + ".stats.byblo")
			statsFileNames.append(statsFileName)
			
			if "byblo_stats" in reuse and isfile(statsFileName):
				if verbose:
					print "   Reusing statistics file " + getPctParamSubstring(statsFileName)
			else:
				if verbose:
					print "   Creating statistics file for " + getPctParamSubstring(statsFileName)
			
				## run Byblo for this sample file
				runTime = runByblo(abspath(fileName), abspath(thesauriDir), \
					abspath(bybloDir), paramStr, verbose)
					
				## write statistics
				resultFileName = join(thesauriDir, basename(fileName) + fullInputFileStr + paramSubstring(paramStr))
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
		print "\n   >> start:runByblo"
	
	## move to Byblo directory and temporary rename input file (so that output files names reflect the parameters used)
	startDir=abspath(os.getcwd())
	fullInputFileStr = sizeSubstring(100.0) if string.find(inputFileName, '%') == -1 else '' 
	tmpInputFileName = inputFileName + fullInputFileStr + paramSubstring(bybloParams)
	os.chdir(bybloDir)
	if verbose:
		print "  Moved to " + os.getcwd()
	os.system("mv " + inputFileName + " " + tmpInputFileName)
	
	## run Byblo in a subprocess and time its execution
	logFile = open(os.devnull, 'w') if not verbose else None
	stime = datetime.datetime.now()
	out = subprocess.call(abspath("./byblo.sh ") + " -i " + tmpInputFileName + " -o " + outputDir +\
		" "+ bybloParams, shell = True, stdout = logFile, stderr = logFile)
	etime = datetime.datetime.now()
	
	## whatever the result, give the input file its name back file and close the log
	os.system("mv " + tmpInputFileName + " " + inputFileName)
	if logFile != None:
		logFile.close()
	if(not out == 0):
		print "   Byblo failed on input file: " + tmpInputFileName + "\n   Fail Code: " + str(out)
		sys.exit()
	
	## move back to initial directory
	os.chdir(startDir)
	if verbose:
		print "  Moved back to " + os.getcwd()
	
	runTime = 1.0*(etime - stime).seconds
	if verbose:
		print "   Run time =", runTime
		print "   >> end:runByblo\n"
	return runTime


## Converts Byblo result files that use skip indexing so that the strings represented by the indexes are restored 
## (slower and heavier but better for result analysis, readability and adaptability)
def generateStringsFiles(sampleFileNames, paramList, thesauriDir, bybloDir, reuse=[], verbose=False):
	if verbose:
		print "\n   >> start:generateStringsFiles"
		
	for fileName in sampleFileNames:
		for paramStr in [paramSubstring(s) for s in paramList]:
			for typeSuffix in [".entries", ".features", ".events", ".sims"]:
				fullInputFileStr = sizeSubstring(100.0) if string.find(fileName, '%') == -1 else '' 
				inputFileName = abspath(join(thesauriDir, basename(fileName)) + fullInputFileStr)
				
				## move to Byblo directory
				startDir=abspath(os.getcwd())
				os.chdir(bybloDir)
				if verbose:
					print "\n   Moved to " + os.getcwd()
				
				## convert both filtered (when it exists and hasn't already been converted) and unfiltered versions
				hasFiltered = isfile(inputFileName+paramStr+typeSuffix)
				for filterSuffix in ([''] if typeSuffix == ".sims" else (['', '.filtered']) if hasFiltered else ['.filtered'] ):
					sourceFileName = inputFileName + paramStr + typeSuffix + filterSuffix
					
					## ...but only if needed
					if "byblo_stats" in reuse and isfile(sourceFileName+".strings"):
						if verbose:
							print "   Reusing strings output file " + getPctParamSubstring(sourceFileName)
					else:
						if verbose:
							print "   Creating strings output file " + getPctParamSubstring(sourceFileName)
					
						logFile = open(os.devnull, 'w') if not verbose else None
						
						out = subprocess.call(abspath("./tools.sh") + " unindex-" + typeSuffix[1:] \
							+ " -i " + sourceFileName\
							+ " -o " + sourceFileName + ".strings" \
							+ (" -Xe "+inputFileName+paramStr+".entry-index" if typeSuffix != ".features" else "") \
							+ (" -Xf "+inputFileName+paramStr+".feature-index" if typeSuffix != ".entries" else "") \
							+ " -et JDBC",\
							shell = True, stdout = logFile, stderr = logFile)
						
						if  logFile != None:
							logFile.close()
				
				## move back to execution directory
				os.chdir(startDir)
				if verbose:
					print "   Moved back to " + os.getcwd()
	if verbose:
		print "   >> end:generateStringsFiles"


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
def fittingMethod(xdata, ydata, method, initialParameters=None, verbose=False):
	colors = ["red", "green", "blue", "purple", "orange", "magenta", "cyan", "limegreen", "gold"]
	## default return values
	fit, params, label, color = [], [], "", ""
	
	## prepare data
	xdata = np.array(xdata)
	ydata = np.array(ydata)
	
	## fitting functions
	powerlaw = lambda x, amp, index: amp * (x**index)
	zipf = lambda x, a: x**(-a)/special.zetac(a)
	diffErr = lambda p, x, y, f: (y - f(p, x))
	
	quadratic = lambda p, x: p[0] * (p[3]*x)**2 + p[1] * p[3]*x + p[2]
	quadraticErr = lambda p, x, y: abs(quadratic(p,x) -y) / abs( (quadratic(p,x) + y) / 2.)
	
	nlogn = lambda p, x: p[0] * (p[3]*x * log(p[3]*x)) + p[1] * p[3]*x + p[2]
	nlognErr = lambda p, x, y: abs(nlogn(p,x) -y) / abs( (nlogn(p,x) + y) / 2.)
	
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
		fit = powerlaw(xdata, amp, index)
		
	## curve_fit on zipf
	elif method == 1:
		y = toProbabilities(ydata) ## convert to probabilities
		label = "curve_fit on zipf"
		pfinal, covar = optimize.curve_fit(zipf, xdata, y)
		s = pfinal[0]
		fit = zipf(xdata, s)
		fit = toFrequencies(fit, ydata) ## restore frequencies
	
	## polyfit on log with powerlaw
	elif method == 2:
		label = "polyfit on log with powerlaw"
		(a, b) = polyfit(log10(xdata), log10(ydata), 1)
		b = 10. ** b
		fit = powerlaw(xdata, b, a)
	
	## polyfit on log with zipf
	elif method == 3:
		y = toProbabilities(ydata) ## convert to probabilities
		label = "polyfit on log with zipf"
		(a, b) = polyfit(log10(xdata), log10(y), 1)
		fit = zipf(xdata, -a)
		fit = toFrequencies(fit, ydata) ## restore frequencies
	
	## leastsq on log with affine + powerlaw
	elif method == 4:
		label = "leastsq on log with affine + powerlaw"
		logx, logy = log10(xdata), log10(ydata)
		pinit = [1.0, -1.0] if initialParameters == None else initialParameters
		out = optimize.leastsq(diffErr, pinit, args=(logx, logy, affine), full_output=1)
		amp, index = 10.0**out[0][0], out[0][1]
		fit = powerlaw(xdata, amp, index)
		
	## leastsq on log with affine + zipf
	elif method == 5:
		y = toProbabilities(ydata) ## convert to probabilities
		label = "leastsq on log with affine + zipf"
		pinit = [1.0, -1.0] if initialParameters == None else initialParameters
		out = optimize.leastsq(diffErr, pinit, args=(log10(xdata), log10(y), affine), full_output=1)
		a =  out[0][1]
		fit = zipf(xdata, -a)
		fit = toFrequencies(fit, ydata) ## restore frequencies
	
	## polyfit / polyval - 2nd degree
	elif method == 6:
		label = "polyfit / polyval - 2nd degree"
		params= polyfit(xdata, ydata, 2)
		fit = polyval(params, xdata)
		
	## leastsq on quadratic
	elif method == 7:
		label = "leastsq on quadratic"
		pinit = [1., 1., 1., 1.]  if initialParameters == None else initialParameters
		out = optimize.leastsq(quadraticErr, pinit, args=(xdata, ydata), full_output=1)
		params = out[0]
		fit = quadratic(params, xdata)
	
	## curve_fit on nlogn
	elif method == 8:
		label = "leastsq on nlogn"
		pinit = [1., 1., 1., 1.]  if initialParameters == None else initialParameters
		out = optimize.leastsq(nlognErr, pinit[:], args=(xdata, ydata), full_output=1) 
		params= out[0] 
		covar = out[1] 
		fit = nlogn(params, xdata)
	
	if verbose:
		print "      Fitting with ["+str(method)+"] " + label
	if method <6:
		return fit, label, colors[method]
	else:
		return fit, label, params


## Generates histograms for Byblo result files (all sizes and parameter strings) containing counts of entries, features 
## and events, as well as thesauri containing similarity values, always combining filtered and unfiltered versions of the
## same file
def generateHistograms(sampleFileNames, paramList, outputDir, bybloDir, reuse=[], verbose=False, cut=False):
	print "\n>> start:generateHistograms "
	
	thesauriDir = join(outputDir, "thesauri")
	graphsDir = join(outputDir, "graphs")
	if not exists(graphsDir):
		os.makedirs(graphsDir)
	
	## file conversion (restore strings to ease comparison)
	generateStringsFiles(sampleFileNames, paramList, thesauriDir, bybloDir, reuse, verbose)
	
	for fileName in sampleFileNames:
		for paramStr in [paramSubstring(s) for s in paramList]:
			for suffix in ['.entries', '.features', '.events', '.sims']:
				fullInputFileStr = sizeSubstring(100.0) if string.find(fileName, '%') == -1 else '' 
				fileBaseName = basename(fileName) + fullInputFileStr + paramStr + suffix
				if "graphs" in reuse and isfile(join(graphsDir, "Histogram-" + fileBaseName + ".pdf")):
					if verbose:
						print "   Reusing histogram for " + getPctParamSubstring(fileBaseName)
				else:
					if verbose:
						print "   Creating histogram for " + getPctParamSubstring(fileBaseName)
					
					## create histogram with a custom step value when required
					if suffix == '.sims':
						createSimilarityHistogram("Similarity", fileBaseName, thesauriDir, graphsDir, verbose)
					else:
						createOccurenceHistogram("Occurence", fileBaseName, thesauriDir, graphsDir, verbose, cut)
	
	print ">> end:generateHistograms"


## Creates a histogram showing the distribution of frequencies for the chosen element (entry, feature, event)
## Both normal and logarithmic scale are created, together with a (for now (very) false) model attempt
def createOccurenceHistogram(label, fileName, thesauriDir, graphsDir, verbose=False, cut=False):
	## generate the histograms
	XBASE, YBASE = 2, 10
	LIMITS = [1, 2.0 ** 64]
	reducedFileSuffix = ".filtered"
	hasFiltered = isfile(join(thesauriDir, fileName))
	bins, hist = extractRowsValues(join(thesauriDir, fileName + (reducedFileSuffix if not hasFiltered else '')), \
		LIMITS, [XBASE, YBASE], verbose=verbose)
	if hasFiltered:
		reducedBins, reducedHist = extractRowsValues(join(thesauriDir, fileName + reducedFileSuffix), 	
			LIMITS, [XBASE, YBASE], verbose=verbose)
	
	## figure set up
	if not cut:
		f, (linearScale, logScale, fitLinearScale, fitLogScale) = pl.subplots(4, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.4) ## add margins
		f.suptitle('Occurence histogram for  ' + fileName[string.rfind(fileName, '.'):] + ' file', fontsize=14, fontweight='bold')
	else:
		individualPlots = [pl.subplots(1, 1) for x in xrange(4)]
		figures = [p[0] for p in individualPlots]
		linearScale, logScale, fitLinearScale, fitLogScale = [p[1] for p in individualPlots]
		for f in figures:
			f.set_size_inches(8.3, 5.8) ## set figure size to A5
			f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2) ## add margins
			f.suptitle('Occurence histogram for  ' + fileName[string.rfind(fileName, '.'):] + ' file', fontsize=14, fontweight='bold')
		
	yLabelPos = [-0.1, 0.5]
	
	## REPRESENT THE DATA
	## linear bar chart of 99% of the mass
	noThreshold = linearScale.bar(bins[:-1], hist, width=bins[1:] - bins[:-1], color="orange", label="No threshold")
	if hasFiltered:
		withThreshold = linearScale.bar(reducedBins[:-1], reducedHist, width=reducedBins[1:] - bins[:-1], \
			color="red", label="With threshold")
	## log-log line
	noThreshold = logScale.loglog(bins[:-1], hist, 'x-', color="orange", label="No threshold")
	if hasFiltered:
		withThreshold = logScale.loglog(reducedBins[:-1], reducedHist, 'x-', color="red", label="With threshold")
	
	pp = max([i for i,x in enumerate((np.cumsum(hist) / sum(hist)) < 0.99) if x])
	decorateGraph(linearScale, 'Linear scale', "number of occurences", "frequency", \
		yLabelPos, "upper right" if hasFiltered else None, largeY=True)
	decorateGraph(logScale, 'Log scale', "number of occurences", "frequency", yLabelPos, "upper right" if hasFiltered else None)
	linearScale.set_xlim(xmax=XBASE ** max(pp, XBASE ** 2))
	
	## FIT THE DATA
	## data
	removeZeroValues = lambda L, L2: [L[i] for i in xrange(len(L)) if L[i]>0 and L2[i]>0]
	x, y = removeZeroValues(bins[:-1], hist), removeZeroValues(hist, bins[:-1])
	fitLinearScale.fill_between(x, y, 0, color="lightgrey")
	fitLogScale.fill_between(x, 1e-50, y, color="lightgrey")
	
	## fit functions (if enough values)
	if len(x) > 4:
		#!>>MODIFY HERE THE METHODS TO USE<<!#
		methods = [1, 3, 4]
		for m in methods:
			mFit, mLabel, mColor= fittingMethod(x, y, m, verbose=verbose)
			fitLinearScale.plot(x, mFit, label=mLabel, color=mColor)
			fitLogScale.loglog(x, mFit, label=mLabel, color=mColor)
	else:
		infoString = "(At least 5 values required to fit a model.)" 
		fitLinearScale.text(0.7, 0.1, infoString, horizontalalignment='center', verticalalignment='center', \
			transform = fitLinearScale.transAxes, fontsize=10, color="firebrick")
		fitLogScale.text(0.7, 0.1, infoString, horizontalalignment='center', verticalalignment='center', \
			transform = fitLogScale.transAxes, fontsize=10, color="firebrick")
		
		
	decorateGraph(fitLinearScale, 'LInear scale - zipfian model', "number of occurences", "frequency", \
		yLabelPos, "upper right", data=(x, y), largeY=True)
	decorateGraph(fitLogScale, 'Log scale - zipfian model', "number of occurences", "frequency", \
		yLabelPos, "upper right", data=(x, y))
	fitLinearScale.set_xlim(xmax=XBASE** max(pp, XBASE ** 2))
	
	if not cut:
		f.savefig(join(graphsDir, 'Histogram-' + fileName + '.pdf'))
	else:
		for i, f in enumerate([p[0] for p in individualPlots]):
			f.savefig(join(graphsDir, 'Histogram-' + fileName + '-' + str(i+1) + '.pdf'))
	pl.close()
	if verbose:
		print ""


## Creates a histogram showing the distribution of similarity scores that appear in a thesaurus
## Only normal scale (values between 0 and 1), no model
def createSimilarityHistogram(label, fileName, thesauriDir, graphsDir, verbose=False):
	## generate the histograms
	XBASE, YBASE = 10, 10
	LIMITS = [0, 1]
	step = 0.01
	reducedFileSuffix = ".neighbours"
	hasFiltered = isfile(join(thesauriDir, fileName + reducedFileSuffix))
	bins, hist = extractRowsValues(join(thesauriDir, fileName + (reducedFileSuffix if not hasFiltered else '')), \
		LIMITS, [XBASE, YBASE], step, verbose)
	if hasFiltered:
		reducedBins, reducedHist = extractRowsValues(join(thesauriDir, fileName + reducedFileSuffix), \
			LIMITS, [XBASE, YBASE], step, verbose)
	
	## figure set up
	f, (sims) = pl.subplots()
	f.set_size_inches(8.3, 5.8) ## set figure size to A5
	f.subplots_adjust(left=0.15, right=0.85) ## add margins
	f.suptitle('Similarity histogram for  ' + fileName[string.rfind(fileName, '.'):] + ' file', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5]
	
	## REPRESENT THE DATA
	## bar chart
	noThreshold = sims.bar(bins[:-1], hist, width=bins[1:] - bins[:-1], color="orange", label="No threshold")
	if hasFiltered:
		withThreshold = sims.bar(reducedBins[:-1], reducedHist, width=reducedBins[1:] - bins[:-1], color="red", label="With threshold")
	decorateGraph(sims, 'LInear scale', "similarity score", "frequency", yLabelPos, "upper right" if hasFiltered else None, largeY=True)

	f.savefig(join(graphsDir, 'Histogram-' + fileName + '.pdf'))
	pl.close()
	if verbose:
		print ""
	

## Creates bins of width "step" using all of the values from the specified file that are relevant for the corresponding histogram
## @return array of bins
def extractRowsValues(fileName, limits, bases=None, step=None, verbose=False):
	## initialisation
	name = ""
	chunk = []
	MIN_VALUE, MAX_VALUE = limits[0], limits[1]
	
	if step:
		bins = np.arange(MIN_VALUE, MAX_VALUE, step)
	else:
		## calculate the min and max powers
		start_power = np.floor(np.log(MIN_VALUE) / np.log(bases[0]))
		end_power = np.ceil(np.log(MAX_VALUE) / np.log(bases[0]))
		##  generate a range of delimiters in log space
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
		print "      Size hist = " +  str(len(hist)) + ", total values = " + str(int(sum(hist))) + "."
	return bins, hist
	

## Generates plots for the impact of changes in the way Byblo is run (input file sizes and parameter strings)
## Representation with either linear plots or bar charts depending on the nature of the element varying
def generatePlots(statsFileNames, paramList, paramListUser, paramListAuto, paramValuesLists,\
	pctList, fileName, outputDir, reuse=[], verbose=False, cut=False):
	
	print "\n>> start:generatePlots"
	graphsDir = join(outputDir, "graphs")
	if not exists(graphsDir):
		os.makedirs(graphsDir)
		
	## lists of percentages and parameter strings used
	pctSubstringsList = [sizeSubstring(p) for p in pctList]
	paramSubstringsList = [paramSubstring(s) for s in paramList]
	
	## create directories to hold all of the statistics information
	sizesVarDict, paramsVarDict= generateStatsDictionaries(statsFileNames, pctSubstringsList, paramSubstringsList)
	if verbose:
		print_lines(sizesVarDict, title="Statistics for varying sizes")
		print_lines(paramsVarDict, title="Statistics for varying parameters")
		
	## study input size variation when several sizes
	if len(pctSubstringsList) > 1:
		for paramStr in paramSubstringsList:
			for plotType in ['files', 'time']:
				fileBaseName = basename(fileName) + paramStr
				if "graphs" in reuse and isfile(join(graphsDir, 'Input-vs-'+plotType+'-'+ fileBaseName+'.pdf')):
					if verbose:
						print "   Reusing " + plotType + " plot for " + paramStr
				else:
					if verbose:
						print "   Creating " + plotType + " plot for " + paramStr
					plotFunction = createPlotInputVsFiles if plotType == "files" else createPlotInputVsTime
					plotFunction(sizesVarDict[paramStr], outputDir, fileBaseName, verbose, cut)
	
	## study Byblo parameter strings variation when several parameter strings
	if len(paramSubstringsList) > 1:
		for pctStr in pctSubstringsList:
			fileBaseName = basename(fileName) + pctStr
			
			## create appropriate dictionary for each graph
			userStringsDict, singleParamDicts = decomposeDictionary(paramsVarDict[pctStr], paramList, paramListUser, paramListAuto, paramValuesLists)
			
			## PARAMETER STRINGS SPECIFIED BY THE USER
			if userStringsDict:
				for plotType in ['files', 'time']:
					if "graphs" in reuse and isfile(join(graphsDir, 'Param-strings-vs-'+plotType+'-'+ fileBaseName+'.pdf')):
						if verbose:
							print "   Reusing param-strings-vs-" + plotType + " plot for " + pctStr
					else:
						if verbose:
							print "   Creating param-strings-vs-" + plotType + " plot for " + pctStr
						plotFunction = createPlotParametersVsFiles if plotType == "files" else createPlotParametersVsTime
						plotFunction(paramListUser, userStringsDict, outputDir, fileBaseName, verbose, cut)
				
			## PARAMETERS VARYING INDIVIDUALLY
			if singleParamDicts:
				# each parameter varying...
				for paramKey, paramValues in zip(singleParamDicts.iterkeys(), [values for values in paramValuesLists if values]):
					# ... in each different configuration
					for confKey in singleParamDicts[paramKey].iterkeys():
						
						for plotType in ['files', 'time']:
							if "graphs" in reuse and isfile(join(graphsDir, paramKey+'-vs-'+plotType+'['+confKey+']'+'-'+ fileBaseName+'.pdf')):
								if verbose:
									print "   Reusing " + paramKey + "-vs-" + plotType + " plot with " + confKey+ " for  " + pctStr 
							else:
								if verbose:
									print "   Creating " + paramKey + "-vs-" + plotType + " plot with " + confKey + " for " + pctStr
								plotFunction = createPlotSingleParamVsFiles if plotType == "files" else createPlotSingleParamVsTime
								plotFunction(paramValues, singleParamDicts[paramKey][confKey], outputDir, fileBaseName, paramKey, confKey, verbose, cut)
		
	print ">> end:generatePlots"
	
	
## Produces a dictionary for all the statistics that appear in the files listed, associating their name with a tuple containing
## all their values, always in the same order
## @return dictionary from the files listed
def statsFileToDictionary(fileNames):
	dictList = []
	## create dictionary for stats in each file
	for f in fileNames:
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


## Generates statistics dictionaries for the complete list of all statistics files by making sublists of files obtained with some common
## element and creating dictionaries for each of these sublists
## @return 2 dictionaries, 1 for each element varying (input file size and Byblo parameters)
def generateStatsDictionaries(fileNames, pctStrings, paramStrings):
	## input sizes varies, parameters are fixed
	sizesVarDict = {}
	for paramStr in paramStrings:
		# don't forget the statistics on input files, that don't contain any parameter string (created BEFORE Byblo run)
		fileList = [fileName for fileName in fileNames if paramStr in fileName or '#' not in fileName]
		thisDict = statsFileToDictionary(fileList)
		sizesVarDict[paramStr] = thisDict
	
	## parameters vary, input file size is fixed
	paramsVarDict= {}
	for pctStr in pctStrings:
		fileList = [fileName for fileName in fileNames if pctStr in fileName]
		thisDict = statsFileToDictionary(fileList)
		paramsVarDict[pctStr] = thisDict
	
	return sizesVarDict, paramsVarDict


##
##
def decomposeDictionary(mainDict, allParamStrings, paramStringsUser, paramStringsAuto, paramValuesLists):
	## PARAMETER STRINGS SPECIFIED BY THE USER
	userDict = {}
	if paramStringsUser:
		for key in mainDict.iterkeys():
			if len(mainDict[key]) == 1:
				userDict[key] = mainDict[key]
			else :
				userDict[key] = (tuple(value for value, paramStr in zip(mainDict[key], allParamStrings) if paramStr in paramStringsUser))
	
	## PARAMETERS VARYING INDIVIDUALLY
	singleParamVarDicts = {}
	if paramStringsAuto:
		paramNames = ["fef", "fff", "fvf", "Smn"]
		for targetValuesList, targetName in zip(paramValuesLists, paramNames):
			if len(targetValuesList) > 1:
				
				## determine substrings that represent fixed settings in each configuration
				list1, list2, list3 = [l for l in paramValuesLists if l != targetValuesList]
				name1, name2, name3 = [n for n in paramNames if n != targetName]
				
				confSettings = []
				if list1 or list2 or list3:
					for substring1 in ['-'+name1+' '+str(n) for n in list1] if list1 else ['']:
						for substring2 in ['-'+name2+' '+str(n) for n in list2] if list2 else ['']:
							for substring3 in ['-'+name3+' '+str(n) for n in list3] if list3 else ['']:
								confSettings.append([s for s in [substring1, substring2, substring3] if s != ''])
				else:
					confSettings.append('')
				
				## build a dictionary for each of these configurations
				singleParamDict = {}
				for settings in confSettings:
					confDict = {}
					for key in mainDict.iterkeys():
						if len(mainDict[key]) == 1:
							confDict[key] = mainDict[key]
						else :
							# select fields corresponding to automatically determined strings containing all of the fixed settings for this conf
							confDict[key] = [value for value, paramStr in zip(mainDict[key], allParamStrings)\
								if paramStr in paramStringsAuto and sum([substring in paramStr for substring in settings]) == len(settings)]
					singleParamDict['"'+' '.join(settings)+'"'] = confDict
				singleParamVarDicts[string.upper(targetName)] = singleParamDict

	return userDict, singleParamVarDicts

## Creates plots showing the relation between the input file size and the result files obtained by running Byblo (size and # of lines)
## Also includes a plot for the relation between entries and events in the input
def createPlotInputVsFiles(statsDictionary, outputDirectory, graphName="", verbose=False, cut=False):
	## figure set up
	if not cut:
		f, (entries, lines, sizes) = pl.subplots(3, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.35) ## add margins
		f.suptitle(graphName+'\nVariation of the input - impact on result files', fontsize=14, fontweight='bold')
	else:
		individualPlots = [pl.subplots(1, 1) for x in xrange(3)]
		figures = [p[0] for p in individualPlots]
		entries, lines, sizes = [p[1] for p in individualPlots]
		for f in figures:
			f.set_size_inches(8.3, 5.8) ## set figure size to A5
			f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2) ## add margins
			f.suptitle(graphName+'\nVariation of the input - impact on result files', fontsize=14, fontweight='bold')
	sizesInputFile, sizeUnit = convertFileSize(statsDictionary["Size_In_Bytes_Of_Input_File"])
	yLabelPos = [-0.1, 0.5] 
	
	## ENTRIES and EVENTS
	numberEntries = statsDictionary["Total_Number_Of_Distinct_Entries"]
	plot1 = entries.plot(sizesInputFile, numberEntries, color='royalblue', linestyle="solid", marker='o', label="entries", alpha=0.6)
	decorateGraph(entries, 'Entries and events', "size ("+sizeUnit+")", "number of entries", yLabelPos, largeY=True)
	
	events = entries.twinx()
	numberEvents = statsDictionary["Total_Number_Of_Distinct_Events"]
	plot2 = events.plot(sizesInputFile, numberEvents, color='red', linestyle="solid", marker='o', label="events", alpha=0.6)
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
	colors = ['RoyalBlue', 'red', 'gold', 'MediumVioletRed', 'DarkOrange', 'Chartreuse']
	for suffix, color in zip(suffixes, colors):
		linesResultFile = statsDictionary["Number_Of_Lines_In_File_"+suffix]
		lines.plot(numberEvents, linesResultFile, color=color, linestyle="dashed", marker='o', label=suffix, alpha=0.6)
		sizesResultFile = convertFileSize(statsDictionary["Size_In_Bytes_Of_File_"+suffix], sizeUnit)
		sizes.plot(numberEvents, sizesResultFile, color=color, linestyle="solid", marker='o', label=suffix, alpha=0.6)
	decorateGraph(lines, 'Lines in the result files', "number of distinct observed events", "number of lines", \
		yLabelPos, "upper left", largeX=True, largeY=True)
	decorateGraph(sizes, 'Size of the result files', "number of distinct observed events", "size ("+sizeUnit+")", \
		yLabelPos, "upper left", largeX=True)
		
	## save figure
	
	if not cut:
		f.savefig(join(outputDirectory, "graphs", 'Input-vs-files-' + graphName + '.pdf'))
	else:
		for i, f in enumerate([p[0] for p in individualPlots]):
			f.savefig(join(outputDirectory, "graphs", 'Input-vs-files-' + graphName + '-' + str(i+1) + '.pdf'))
	pl.close()
	if verbose:
		print ""


## Creates plots showing the relation between the input file size and the run time of Byblo 
## Also includes a modelling attempt of the run time depending on the # of entries and then events
def createPlotInputVsTime(statsDictionary, outputDirectory, graphName="", verbose=False, cut=False):
	## figure set up
	if not cut:
		f, (entries, events) = pl.subplots(2, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.35) ## add margins
		f.suptitle('Variation of the input file size\nImpact on run time', fontsize=14, fontweight='bold')
	else:
		individualPlots = [pl.subplots(1, 1) for x in xrange(2)]
		figures = [p[0] for p in individualPlots]
		entries, events = [p[1] for p in individualPlots]
		for f in figures:
			f.set_size_inches(8.3, 5.8) ## set figure size to A5
			f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2) ## add margins
			f.suptitle('Variation of the input file size\nImpact on run time', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5] 
	
	## REPRESENT THE DATA
	times,  timeUnit = convertTimeRange(statsDictionary["Byblo_Run_Time"])
	numberEntries = statsDictionary["Total_Number_Of_Distinct_Entries"]
	entries.plot(numberEntries, times, color='aquamarine', label="run time function of entries", marker='o', linestyle='dashed')
	numberEvents = statsDictionary["Total_Number_Of_Distinct_Events"]
	events.plot(numberEvents, times, color='aquamarine', label="run time function of events", marker='o', linestyle='dashed')
	
	## FIT THE DATA (if enough values)
	if len(times) > 3:
		quadratic = lambda p, x: p[0] * (p[3]*x)**2 + p[1] * p[3]*x + p[2]
		nlogn = lambda p, x: p[0] * (p[3]*x * log(p[3]*x)) + p[1] * p[3]*x + p[2]
		
		## entries >> method 7 to fit a quadratic shape
		pInit = [1., 1., 1., 1./max(numberEntries)] # inital guess
		mFit, mLabel, mParams1= fittingMethod(numberEntries, times, 7, pInit, verbose)
		entries.plot(numberEntries, mFit, label=mLabel, color='firebrick', linestyle='dashed')
		nMin, nMax = min(numberEntries), max(numberEntries)
		steps = np.arange(nMin, nMax, (nMax-nMin)/100.) ## for a smoother model
		entries.plot(steps, quadratic(mParams1, steps), label="smoothed model", color='red')
		
		## events >> method 8 to fit an "n log n" shape
		pInit = [1., 1., 1., 1./max(numberEvents)] # inital guess
		mFit, mLabel, mParams2= fittingMethod(numberEvents, times, 8, pInit, verbose) # expected curve
		events.plot(numberEvents, mFit, label=mLabel, color='firebrick', linestyle='dashed')
		nMin, nMax = min(numberEvents), max(numberEvents)
		steps = np.arange(nMin, nMax, (nMax-nMin)/100.) # and also a smoother curve
		events.plot(steps, nlogn(mParams2, steps), label="smoothed model", color='red')
	
		## DISPLAY APPROXIMATIONS
		## entries
		a, b, c, ratio = mParams1
		infoString = "Run time approximation:\nt(n) = " \
			+ "%.3f" % (a)  + " x N**2 + " +  "%.3f" % (b)  + " x N + " + "%.3f" % (c) \
			+ "\nwith N = n x " + "%.5e" % (ratio)
		entries.text(0.4, 0.85, infoString, horizontalalignment='center', verticalalignment='center', \
			transform = entries.transAxes, fontsize=10, color="dimGrey")
			
		## events
		a, b, c, ratio = mParams2
		infoString = "Run time approximation from events:" \
			+ "\nt(n) = " + "%.3f" % (a) + " x N log(N) + " + "%.3f" % (b) + " x N + " + "%.3f" % (c) \
			+ "\nwith N = n x " + "%.5e" % (ratio)
		events.text(0.4, 0.85, infoString, horizontalalignment='center', verticalalignment='center', \
			transform = events.transAxes, fontsize=10, color="dimGrey")
	else:
		infoString = "(At least 4 values required to fit a model.)" 
		entries.text(0.7, 0.1, infoString, horizontalalignment='center', verticalalignment='center', \
			transform = entries.transAxes, fontsize=10, color="firebrick")
		events.text(0.7, 0.1, infoString, horizontalalignment='center', verticalalignment='center', \
			transform = events.transAxes, fontsize=10, color="firebrick")
		
		
	decorateGraph(entries, "Run time based on entries", "number of distinct entries", "time ("+timeUnit+")", yLabelPos, \
		"lower right", largeX=True, largeY=True)
	decorateGraph(events, "Run time based on events", "number of distinct observed events", "time ("+timeUnit+")", yLabelPos, \
		"lower right", largeX=True, largeY=True)
		
	## save figure
	if not cut:
		f.savefig(join(outputDirectory, "graphs", 'Input-vs-time-' + graphName + '.pdf'))
	else:
		for i, f in enumerate([p[0] for p in individualPlots]):
			f.savefig(join(outputDirectory, "graphs", 'Input-vs-time-' + graphName + '-' + str(i+1) + '.pdf'))
	pl.close()
	if verbose:
		print ""


## Creates plots showing the relation between the parameters used and the result files obtained by running Byblo (size and # of lines)
## No model as the parameters (x values) are not continuous data
def createPlotParametersVsFiles(paramList, statsDictionary, outputDirectory, graphName="", verbose=False, cut=False):
	## figure set up
	if not cut:
		f, (lines, sizes) = pl.subplots(2, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.45, bottom=0.15) ## add margins
		f.suptitle(graphName+'\nVariation of the parameters - Impact on result files', fontsize=14, fontweight='bold')
	else:
		individualPlots = [pl.subplots(1, 1) for x in xrange(2)]
		figures = [p[0] for p in individualPlots]
		lines, sizes = [p[1] for p in individualPlots]
		for f in figures:
			f.set_size_inches(8.3, 5.8) ## set figure size to A5
			f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2) ## add margins
			f.suptitle(graphName+'\nVariation of the parameters - Impact on result files', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5] 
	
	## draw LINES and SIZE plots for each result file
	suffixes = ['.entries.filtered', '.events.filtered', '.sims.neighbours']
	colors = ['RoyalBlue', 'red', 'gold', 'MediumVioletRed', 'DarkOrange', 'Chartreuse']
	ind = np.arange(len(paramList))
	width = 0.5 / len(paramList)
	## determine most appropriate unit
	array, sizeUnit = convertFileSize([s for suffix in suffixes for s in statsDictionary["Size_In_Bytes_Of_File_"+suffix]])
	
	for i, (suffix, color) in enumerate(zip(suffixes, colors)):
		linesResultFile  = convertFileSize(statsDictionary["Number_Of_Lines_In_File_"+suffix], sizeUnit)
		lines.bar(ind+width*i, linesResultFile, width=width, align='center', color=color, label=suffix)
		sizesResultFile = convertFileSize(statsDictionary["Size_In_Bytes_Of_File_"+suffix], sizeUnit)
		sizes.bar(ind+width*i, sizesResultFile, width=width, align='center', color=color, label=suffix)
		
	decorateGraph(lines, 'Lines in the result files', "parameter string for Byblo", "number of lines", \
		yLabelPos, "upper left", largeX=True, largeY=True)
	decorateGraph(sizes, 'Size of the result files', "parameter string for Byblo", "size ("+sizeUnit+")", \
		yLabelPos, "upper left", largeX=True)
	
	lines.set_xticks(ind+len(paramList)*width*0.5)
	sizes.set_xticks(ind+len(paramList)*width*0.5)
	lines.set_xticklabels(paramList, rotation=15, size='small')
	sizes.set_xticklabels(paramList, rotation=15, size='small')

	## save figure
	if not cut:
		f.savefig(join(outputDirectory, "graphs", 'Param-strings-vs-files-' + graphName + '.pdf'))
	else:
		for i, f in enumerate([p[0] for p in individualPlots]):
			f.savefig(join(outputDirectory, "graphs", 'Param-strings-vs-files-' + graphName + '-' + str(i+1) + '.pdf'))
	pl.close()
	if verbose:
		print ""


## Creates a plot showing the relation between the parameters used and the run time of Byblo
## No model as the parameters (x values) are not continuous data
def createPlotParametersVsTime(paramList, statsDictionary, outputDirectory, graphName="", verbose=False, cut=False):
	## figure set up
	f, (time) = pl.subplots(1, 1)
	f.set_size_inches(8.3, 5.8) ## set figure size to A5
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2) ## add margins
	f.suptitle(graphName+'\nVariation of the parameters - Impact on run time', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5] 
	
	## REPRESENT THE DATA
	ind = np.arange(len(paramList))
	width = 0.5 / len(paramList)
	times,  timeUnit = convertTimeRange(statsDictionary["Byblo_Run_Time"])
	time.bar(ind+width, times, width=width, align='center', color='aquamarine', label="run time on parameters")
		
	decorateGraph(time, "Run time based on Byblo parameters", "parameter string for Byblo", "time ("+timeUnit+")", \
		yLabelPos, largeY=True)
	
	time.set_xticks(ind+len(paramList)*width*0.5)
	time.set_xticklabels(paramList, rotation=15, size='small')

	## save figure
	f.savefig(join(outputDirectory, "graphs", 'Param-strings-vs-time-' + graphName + '.pdf'))
	pl.close()
	if verbose:
		print ""


## 
##
def createPlotSingleParamVsFiles(paramList, statsDictionary, outputDirectory, graphName="", \
	paramName="", confName="", verbose=False, cut=False):
	
	## figure set up
	if not cut:
		f, (lines, sizes) = pl.subplots(2, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.45, bottom=0.15) ## add margins
		f.suptitle(graphName+' ['+confName+']'+'\nVariation of '+paramName+' - impact on result files', fontsize=14, fontweight='bold')
	else:
		individualPlots = [pl.subplots(1, 1) for x in xrange(2)]
		figures = [p[0] for p in individualPlots]
		lines, sizes = [p[1] for p in individualPlots]
		for f in figures:
			f.set_size_inches(8.3, 5.8) ## set figure size to A5
			f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2) ## add margins
			f.suptitle(graphName+' ['+confName+']'+'\nVariation of '+paramName+' - impact on result files', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5] 
	
	## draw LINES and SIZE plots for each result file
	suffixes = ['.entries.filtered', '.events.filtered', '.sims.neighbours']
	colors = ['RoyalBlue', 'red', 'gold', 'MediumVioletRed', 'DarkOrange', 'Chartreuse']
	uselessArray, sizeUnit = convertFileSize([s for suffix in suffixes for s in statsDictionary["Size_In_Bytes_Of_File_"+suffix]]) # find most appropriate unit
	ind = np.arange(len(paramList))
	width = 0.5 / len(paramList)
	for suffix, color in zip(suffixes, colors):
		linesResultFile  = statsDictionary["Number_Of_Lines_In_File_"+suffix]
		lines.plot(paramList, linesResultFile, color=color, linestyle="dashed", marker='o', label=suffix, alpha=0.6)
		sizesResultFile = convertFileSize(statsDictionary["Size_In_Bytes_Of_File_"+suffix], sizeUnit)
		sizes.plot(paramList, sizesResultFile, color=color, linestyle="solid", marker='o', label=suffix, alpha=0.6)
		
	decorateGraph(lines, "Lines in the result files based on "+fullParamName(paramName), string.lower(paramName)+" value", "number of lines", \
		yLabelPos, "upper left", largeX=True, largeY=True)
	decorateGraph(sizes, "Size of the result files based on "+fullParamName(paramName), string.lower(paramName)+" value", "size ("+sizeUnit+")", \
		yLabelPos, "upper left", largeX=True)

	## save figure
	if not cut:
		f.savefig(join(outputDirectory, "graphs", paramName + '-vs-files'+'['+confName+']'+'-' + graphName + '.pdf'))
	else:
		for i, f in enumerate([p[0] for p in individualPlots]):
			f.savefig(join(outputDirectory, "graphs", paramName + '-vs-files'+'['+confName+']'+'-' + graphName + '-' + str(i+1) + '.pdf'))
	pl.close()
	if verbose:
		print ""


##
##
def createPlotSingleParamVsTime(paramList, statsDictionary, outputDirectory, graphName="", \
	paramName="", confName="", verbose=False, cut=False):
	
	## figure set up
	f, (time) = pl.subplots(1, 1)
	f.set_size_inches(8.3, 5.8) ## set figure size to A5
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2) ## add margins
	f.suptitle(graphName+' ['+confName+']'+'\nVariation of '+paramName+' - impact on run time', fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5] 
	
	## REPRESENT THE DATA
	times,  timeUnit = convertTimeRange(statsDictionary["Byblo_Run_Time"])
	time.plot(paramList, times, color='aquamarine', label="run time function of "+paramName, marker='o', linestyle='dashed')
		
	decorateGraph(time, "Run time based on "+fullParamName(paramName), string.lower(paramName)+" value", "time ("+timeUnit+")", \
		yLabelPos, largeY=True)

	## save figure
	f.savefig(join(outputDirectory, "graphs", paramName + '-vs-time'+'['+confName+']'+'-' + graphName + '.pdf'))
	pl.close()
	if verbose:
		print ""



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
	
	
## Deletes intermediary files needed for the statistics generation when they are in the list of item types (samples, 
## thesauri, event_stats, byblo_stats) that the user specified for deletion
def deleteOnExit(deleteList, outputDirectory, sampleFiles=[], originalInputFile="", verbose=False):
	print "\n>> start:deleteOnExit "
	if(deleteList != ["nothing"]):
		for directory in deleteList:
			os.system("rm -r " + join(outputDirectory, directory))
			if verbose:
				print "   Deleted directory " + directory
	print ">> end:deleteOnExit"


## Prints an array or a dictionary with customizable start, end, line length and title
def print_lines(list, min=0, max=0, line_max=0, title="List"):
	if max == 0:
		max = len(list)
	print "\n" + title
	if type(list) is dict:
		for i, key in enumerate(sorted(list.iterkeys())) :
			if i<max:
				if type(list[key]) is dict:
					print_lines(list[key], title=key)
				else:
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
		
	## PARSE COMMAND LINE
	parser = argparse.ArgumentParser(description="Generate statistics for Byblo.", prog="byblo_statistics")

	# data file name
	parser.add_argument('data', metavar='file', action='store', default="",
		help='data file containing  a feature / entry set')
	# percentages of entries to take into account
	parser.add_argument('-p', '--percentage', metavar='x', type=nonZeroPct, nargs='+', dest='pctList', 
		action='store', default=[100.],
		help='list of approximate percentages (]0.-100.]) of entries to keep  for a set of statistics '+\
		'\n(default: 100%%)')
	# thresholds for the entry frequency
	parser.add_argument('-fef', '--filter-entry-freq', metavar='x', type=positiveInt, nargs='*', dest='fef', 
		action='store', default=[],
		help='list of values consecutively used by Byblo as a minimum entry frequency threshold'+\
		'\n(default: 0.0)')
	# thresholds for the feature frequency
	parser.add_argument('-fff', '--filter-feature-freq', metavar='x', type=positiveInt, nargs='*', dest='fff', 
		action='store', default=[],
		help='list of values consecutively used by Byblo as a minimum feature frequency threshold'+\
		'\n(default: 0.0)')
	# thresholds for the event frequency
	parser.add_argument('-fvf', '--filter-event-freq', metavar='x', type=positiveInt, nargs='*', dest='fvf', 
		action='store', default=[],
		help='list of values consecutively used by Byblo as a minimum event frequency threshold'+\
		'\n(default: 0.0)')
	# thresholds for the similarity scores
	parser.add_argument('-Smn', '--similarity-min', metavar='x', type=zeroToOneFloat, nargs='*', dest='Smn', 
		action='store', default=[],
		help='list of values consecutively used by Byblo as a minimum similarity threshold'+\
		'\n(default: 0.0)')
	# Byblo parameters
	parser.add_argument('-P', '--Byblo-params', metavar='string', nargs='*', dest='paramStringsUser', 
		action='store', default=[],
		help='list of Byblo parameter strings consecutively passed to the script (only considered as a block, no  separate study of each parameter), '+\
		'use Byblo-help option for more info (default: "")')
	# Byblo help option
	parser.add_argument('-H', '--Byblo-help', dest='bybloHelp', action='store_true',
		help='display Byblo help, it is then possible to modify Byblo parameters (default: False)')
	# output directory
	parser.add_argument('-O', '--output-dir', metavar='path', dest='outputDir', 
		action='store', default=".",
		help='output directory where sub-directories "samples", "stats", "thesauri" and "graphs" can be found'+\
		'(default: ".")')
	# Byblo directory
	parser.add_argument('-B', '--Byblo-dir', metavar='path', dest='bybloDir', 
		action='store', default="../Byblo-2.0.1",
		help='Byblo directory, usually named "Byblo-x.x.x" '+\
		'(default: "../Byblo-2.0.1")')
	# reuse option (for files and directories created in the process of the graph generation)
	reuseOptions = ['samples', 'events_stats', 'byblo_stats', 'graphs']
	parser.add_argument('-r', '--reuse', metavar='string', nargs='*', dest='reuse',
		action='store', choices=reuseOptions,  default = ["nothing"],
		help='reuse specified items from "samples", "events_stats", "byblo_stats" and/or "graphs" during execution, '+\
		' everything when none specifed (default: reuse nothing)')
	# delete option (for files and directories created in the process of the graph generation)
	deleteOptions = ['samples', 'thesauri', 'stats']
	parser.add_argument('-d', '--delete', metavar='string', nargs='*', dest='delete',
		action='store', choices=deleteOptions, default = ["nothing"],
		help='delete specified items from "samples", "thesauri" and/or "stats" before exiting, '+\
		' everything when none specifed (default: delete nothing)')
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
		help='display Byblo output and statistics information for each step (default: False)')
	# cut option (for graphs)
	parser.add_argument('-c', '--cut', dest='cut', action='store_true',
		help='cut graphs so that all subgraphs are in distinct files (default: False)')
	
	## PROCESS PARAMETERS
	a = parser.parse_args()
	[a.pctList, a.fef, a.fff, a.fvf, a.Smn] = [sorted(set(list)) for list in [a.pctList, a.fef, a.fff, a.fvf, a.Smn] ]
	if not a.reuse:
		a.reuse = reuseOptions
	if not a.delete:
		a.delete = deleteOptions
	
	
	print "data\t", a.data
	print "p\t", a.pctList
	print "fef\t", a.fef
	print "fff\t", a.fff
	print "fvf\t", a.fvf
	print "Smn\t", a.Smn
	print "paramStringsUser\t", a.paramStringsUser
	print "bybloHelp\t", a.bybloHelp
	print "bybloDir\t", a.bybloDir
	print  "outputDir\t", a.outputDir
	print "reuse\t", a.reuse
	print "delete\t", a.delete
	print "verbose\t", a.verbose
	print "cut\t", a.cut
	print ""
	
	
	## start operations
	stime = datetime.datetime.now()
	print "***************************************************************************"
	print "BYBLO STATISTICS TOOL"
	print "***************************************************************************\n"

	
	## BYBLO HELP +  PARAMETER STRINGS
	if a.bybloHelp:
		a.paramStringsUser = displayBybloHelp(a.bybloDir, a.paramStringsUser, a.verbose)
	paramStringsAuto = generateParameterStrings(a.fef, a.fff, a.fvf, a.Smn, a.verbose)
	allParamStrings = sorted(set(a.paramStringsUser + paramStringsAuto))
	if not allParamStrings:
		allParamStrings.append("")
		
	## EVENTS STATS + SAMPLES
	sampleFileNames, statsFileNames = eventsStats(a.data, a.outputDir, a.pctList, a.reuse, a.verbose)
	if a.verbose:
		print_lines(statsFileNames, title="Statistics files after feature extraction")
	
	## BYBLO STATS (both auto and user patameter strings)
	statsFileNames += bybloStats(sampleFileNames, a.outputDir, a.bybloDir, allParamStrings, a.reuse, a.verbose)
	if a.verbose:
		print_lines(statsFileNames, title="Statistics files after byblo run")
	
	## HISTOGRAMS
	generateHistograms(sampleFileNames, allParamStrings, a.outputDir, a.bybloDir, a.reuse, a.verbose, a.cut)
	
	## PLOTS
	singleParameterValuesLists = [a.fef, a.fff, a.fvf, a.Smn]
	generatePlots(statsFileNames, allParamStrings, a.paramStringsUser, paramStringsAuto, singleParameterValuesLists, \
		a.pctList, a.data, a.outputDir, a.reuse, a.verbose, a.cut)
	
	## CLEAN UP
	deleteOnExit(a.delete, a.outputDir, sampleFileNames, a.data, a.verbose)
	
	etime = datetime.datetime.now()
	print "\n>Execution took", etime-stime, "hours"