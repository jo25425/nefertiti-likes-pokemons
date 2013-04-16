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
__version__ = "2.0"
__maintainer__ = "Joanne Robert"
__email__ = "jr317@sussex.ac.uk"
__status__ = "Development"
#~ -----------------------------------------------------------------------------

import sys
import argparse
import subprocess
import datetime
import random
import string
import numpy as np
import pylab as pl
import os
from os.path import *
from scipy import *
from scipy import optimize, stats, special
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


## Produces a string that gives information about parameters used and can be in a file name
## @return formatted parameters string
def paramSubstring(str):
	return '#' + string.replace(str, ' ', '_')
	
## Extracts from a file name the substrings produced by the 2 methods above
## @return string reflecting parameters chosen
def getParamSubstring(name):
	end = basename(name)
	index = string.find(end, '#')
	if index == -1: return index
	else: return "..." + end[index:]

## Transforms the abbreviated name of a Byblo parameter (e.g. "fef") into its full name
## (e.g. "entry frequency filter")
## @return complete name of the parameter
def fullParamName(initials):
	return { 'fef': 'entry frequency filter',
			'fff': 'feature frequency filter',
			'fvf': 'event frequency filter',
			'smn': 'minimum similarity '
			}[string.lower(initials)]


## Converts Byblo result files that use skip indexing so that the strings represented by the indexes 
## are restored (slower and heavier but better for result analysis, readability and adaptability)
def generateStringsFiles(sampleFileNames, paramList, thesauriDir, bybloDir, 
	reuse=[], verbose=False, filesToCreate=[]):
	savedFiles = []
	
	def vprint(s): print ("[INFO] " + s + '\n' if verbose else ""),
	
	for fileName in sampleFileNames:
		for paramStr in [paramSubstring(s) for s in paramList]:
			for typeSuffix in [".entries", ".features", ".events", ".sims"]:
				inputFileName = abspath(join(thesauriDir, basename(fileName)))
				
				## move to Byblo directory
				startDir=abspath(os.getcwd())
				os.chdir(bybloDir)
				vprint("Moved to " + os.getcwd())
				
				## convert both filtered (when it exists and hasn't already been converted) and 
				## unfiltered versions
				hasFiltered = isfile(inputFileName+paramStr+typeSuffix)
				for filterSuffix in ([''] if typeSuffix == ".sims" else 
					(['', '.filtered']) if hasFiltered else ['.filtered'] ):
					srcFile = inputFileName + paramStr + typeSuffix + filterSuffix
					
					## only convert wanted files
					if not filesToCreate or srcFile in filesToCreate:
						## ...but only if needed
						if "byblo_stats" in reuse and isfile(srcFile+".strings"):
							vprint("Reusing converted file " + getParamSubstring(srcFile))
						else:
							vprint("Restoring strings for file " + getParamSubstring(srcFile))
						
							logFile = open(os.devnull, 'w') if not verbose else None
							out = subprocess.call(abspath("./tools.sh") + " unindex-" + typeSuffix[1:] \
								+ " -i " + srcFile
								+ " -o " + srcFile + ".strings" 
								+ (" -Xe "+inputFileName+paramStr+".entry-index" 
									if typeSuffix != ".features" else "")
								+ (" -Xf "+inputFileName+paramStr+".feature-index" 
									if typeSuffix != ".entries" else "")
								+ " -et JDBM",
								shell = True, stdout = logFile, stderr = logFile)
							if  logFile: logFile.close()
								
						savedFiles.append(srcFile + ".strings")
				## move back to execution directory
				os.chdir(startDir)
				vprint("Moved back to " + os.getcwd() + '\n')
	return savedFiles


## Makes a graph nicer and clearer by adding a title, axes labels, a legend, limits on the axes (for an 
## exact fit or a "zoom"), a scientific notation (for large numbers)
def decorateGraph(subplot, title="", xLabel="", yLabel="", yLabelPos=None, legendPos=None,
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
	
	
## Offers a selection of fitting methods and model functions that allow to model some of the relations
## that exist within the data handled here
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
		label = "leastsq on log (for parameters) and powerlaw (for drawing)"
		logx, logy = log10(xdata), log10(ydata)
		pinit = [1.0, -1.0] if initialParameters == None else initialParameters
		out = optimize.leastsq(diffErr, pinit, args=(logx, logy, affine), full_output=1)
		amp, index = 10.0**out[0][0], out[0][1]
		fit = powerlaw(xdata, amp, index)
		
	## leastsq on log with affine + zipf
	elif method == 5:
		y = toProbabilities(ydata) ## convert to probabilities
		label = "leastsq on log and zipfian model"
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
	
	if verbose: print "[INFO] Modelling data using "+label+"[method "+str(method)+"]"
	return fit, label, colors[method]


## Generates histograms for Byblo result files (all sizes and parameter strings) containing counts of 
## entries, features and events, as well as thesauri containing similarity values, always combining 
## filtered and unfiltered versions of the same file
def generateHistograms(sampleFileNames, paramList, outputDir, bybloDir, reuse=[], 
	verbose=False, cut=False):
	thesauriDir = join(outputDir, "thesauri")
	graphsDir = join(outputDir, "graphs")
	if not exists(graphsDir): os.makedirs(graphsDir)
	
	def vprint(s): print ("[INFO] " + s + '\n' if verbose else ""),
	
	## file conversion (restore strings to ease comparison)
	generateStringsFiles(sampleFileNames, paramList, thesauriDir, bybloDir, reuse, verbose)
	
	savedFiles = []
	for fileName in sampleFileNames:
		for paramStr in [paramSubstring(s) for s in paramList]:
			for suffix in ['.entries', '.features', '.events', '.sims']:
				fileBaseName = basename(fileName) + paramStr + suffix
				
				if "graphs" in reuse and isfile(join(graphsDir, "Histogram-" + fileBaseName + ".pdf")):
					vprint("Reusing histogram for " + getParamSubstring(fileBaseName))
				else:
					vprint("Creating histogram for " + getParamSubstring(fileBaseName))
					
					## create histogram with a custom step value when required
					if suffix == '.sims':
						savedFiles += createSimilarityHistogram("Similarity", fileBaseName, 
							thesauriDir, graphsDir, verbose)
					else:
						savedFiles += createOccurenceHistogram("Occurence", fileBaseName, 
							thesauriDir, graphsDir, verbose, cut)
					if verbose: print ""
	return savedFiles


## Creates a histogram showing the distribution of frequencies for the chosen item (entry, feature, event)
## Both normal and logarithmic scale are created, together with a (for now (very) false) model attempt
## @return names of the saved files
def createOccurenceHistogram(label, fileName, thesauriDir, graphsDir, verbose=False, cut=False):
	## 0) Preparation of the data
	XBASE, YBASE = 2, 10
	LIMITS = [1, 2.0 ** 64]
	reducedFileSuffix = ".filtered"
	hasFiltered = isfile(join(thesauriDir, fileName))
	bins, hist = extractRowsValues(join(thesauriDir, fileName + (reducedFileSuffix if not hasFiltered else '')), 
		LIMITS, [XBASE, YBASE], verbose=verbose)
	if hasFiltered:
		reducedBins, reducedHist = extractRowsValues(join(thesauriDir, fileName + reducedFileSuffix), 	
			LIMITS, [XBASE, YBASE], verbose=verbose)
	
	## 1) Figure set up
	if not cut:
		f, (linearScale, logScale, fitLinearScale, fitLogScale) = pl.subplots(4, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.4) ## margins
		f.suptitle("Frequency distribution for " + fileName[string.rfind(fileName, '.')+1:], 
			fontsize=14, fontweight='bold')
	else:
		individualPlots = [pl.subplots(1, 1) for x in xrange(4)]
		figures = [p[0] for p in individualPlots]
		linearScale, logScale, fitLinearScale, fitLogScale = [p[1] for p in individualPlots]
		for f in figures:
			f.set_size_inches(8.3, 5.8) ## set figure size to A5
			f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2) ## margins
			f.suptitle("Frequency distribution for " + fileName[string.rfind(fileName, '.')+1:], 
				fontsize=14, fontweight='bold')
		
	yLabelPos = [-0.1, 0.5]
	
	## 2) Data representation
	## linear bar chart of 99% of the mass
	noThreshold = linearScale.bar(bins[:-1], hist, width=bins[1:] - bins[:-1], 
		color="orange", label="No threshold")
	if hasFiltered:
		withThreshold = linearScale.bar(reducedBins[:-1], reducedHist, 
			width=reducedBins[1:] - bins[:-1], color="red", label="With threshold")
	## log-log line
	noThreshold = logScale.loglog(bins[:-1], hist, 'x-', color="orange", label="No threshold")
	if hasFiltered:
		withThreshold = logScale.loglog(reducedBins[:-1], reducedHist, 'x-', 
			color="red", label="With threshold")
	
	pp = max([i for i,x in enumerate((np.cumsum(hist) / sum(hist)) < 0.99) if x])
	decorateGraph(linearScale, 'Linear scale', "number of occurences", "frequency",
		yLabelPos, "upper right" if hasFiltered else None, largeY=True)
	decorateGraph(logScale, 'Logarithmic scale', "number of occurences", "frequency", 
		yLabelPos, "upper right" if hasFiltered else None)
	linearScale.set_xlim(xmax=XBASE ** max(pp, XBASE ** 2))
	
	## 3) Model fitting
	## data
	removeZeroValues = lambda L, L2: [L[i] for i in xrange(len(L)) if L[i]>0 and L2[i]>0]
	x, y = removeZeroValues(bins[:-1], hist), removeZeroValues(hist, bins[:-1])
	fitLinearScale.fill_between(x, y, 0, color="lightgrey")
	fitLogScale.fill_between(x, 1e-50, y, color="lightgrey")
	
	## fit functions (if enough values)
	if len(x) > 4:
		#!>>MODIFY HERE THE METHODS TO USE<<!#
		methods = [1, 4]
		for m in methods:
			mFit, mLabel, mColor= fittingMethod(x, y, m, verbose=verbose)
			fitLinearScale.plot(x, mFit, label=mLabel, color=mColor)
			fitLogScale.loglog(x, mFit, label=mLabel, color=mColor)
	else:
		infoString = "(At least 5 values required to fit a model.)" 
		fitLinearScale.text(0.7, 0.1, infoString, horizontalalignment='center', verticalalignment='center',
			transform = fitLinearScale.transAxes, fontsize=10, color="firebrick")
		fitLogScale.text(0.7, 0.1, infoString, horizontalalignment='center', verticalalignment='center',
			transform = fitLogScale.transAxes, fontsize=10, color="firebrick")
		
	decorateGraph(fitLinearScale, 'Linear scale - zipfian model', "number of occurences", 
		"frequency", yLabelPos, "upper right", data=(x, y), largeY=True)
	decorateGraph(fitLogScale, 'Logarithmic scale - zipfian model', "number of occurences", 
		"frequency", yLabelPos, "upper right", data=(x, y))
	fitLinearScale.set_xlim(xmax=XBASE** max(pp, XBASE ** 2))
	
	## 4) Save file(s)
	savedFiles = []
	if not cut:
		f.savefig(join(graphsDir, 'Histogram-' + fileName + '.pdf'))
		savedFiles.append(join(graphsDir, 'Histogram-' + fileName + '.pdf'))
	else:
		for i, f in enumerate([p[0] for p in individualPlots]):
			f.savefig(join(graphsDir, 'Histogram-' + fileName + '-' + str(i+1) + '.pdf'))
			savedFiles.append(join(graphsDir, 'Histogram-' + fileName + '-' + str(i+1) + '.pdf'))
	pl.close()
	return savedFiles


## Creates a histogram showing the distribution of similarity scores that appear in a thesaurus
## Only normal scale (values between 0 and 1), no model
## @return names of saved files
def createSimilarityHistogram(label, fileName, thesauriDir, graphsDir, verbose=False):
	## 0)Preparation of the data
	XBASE, YBASE = 10, 10
	LIMITS = [0, 1]
	step = 0.01
	reducedFileSuffix = ".neighbours"
	hasFiltered = isfile(join(thesauriDir, fileName + reducedFileSuffix))
	bins, hist = extractRowsValues(join(thesauriDir, fileName + (reducedFileSuffix if not hasFiltered else '')),
		LIMITS, [XBASE, YBASE], step, verbose)
	if hasFiltered:
		reducedBins, reducedHist = extractRowsValues(join(thesauriDir, fileName + reducedFileSuffix),
			LIMITS, [XBASE, YBASE], step, verbose)
	
	## 1) Figure set up
	f, (sims) = pl.subplots()
	f.set_size_inches(8.3, 5.8) ## set figure size to A5
	f.subplots_adjust(left=0.15, right=0.85) ## add margins
	f.suptitle('Similarity score distribution for ' + fileName[string.rfind(fileName, '.')+1:], 
		fontsize=14, fontweight='bold')
	yLabelPos = [-0.1, 0.5]
	
	## 2) Representation of the data
	noThreshold = sims.bar(bins[:-1], hist, width=bins[1:] - bins[:-1], 
		color="orange", label="No threshold")
	if hasFiltered:
		withThreshold = sims.bar(reducedBins[:-1], reducedHist, 
			width=reducedBins[1:] - bins[:-1], color="red", label="With threshold")
	decorateGraph(sims, 'Linear scale', "similarity score", "frequency", 
		yLabelPos, "upper right" if hasFiltered else None, largeY=True)

	## 3) Save file
	f.savefig(join(graphsDir, 'Histogram-' + fileName + '.pdf'))
	pl.close()
	return [join(graphsDir, 'Histogram-' + fileName + '.pdf')]
	

## Creates bins of width "step" using all of the values from the specified file that are relevant 
## for the corresponding histogram
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
	
	if verbose: print "[INFO] " + str(int(sum(hist))) + " values extracted for the histogram"
	return bins, hist
	

## Generates plots for the impact of changes in the way Byblo is run (input file sizes and parameter strings)
## Representation with either linear plots or bar charts depending on the nature of the element varying
def generatePlots(statsDict, paramList, paramListUser, paramListAuto, paramValuesLists,\
	fileName, outputDir, reuse=[], verbose=False, cut=False):
	if len(paramList) < 2: return
	def vprint(s): print ("[INFO] " + s + '\n' if verbose else ""),
	
	graphsDir = join(outputDir, "graphs")
	if not exists(graphsDir): os.makedirs(graphsDir)
	fileBaseName = basename(fileName)
	savedFiles = []
	
	## create appropriate dictionary for each graph (sorted on parameter values)
	fullStringsDict, singleParamDicts, paramValuesLists = decomposeDictionary(statsDict, paramList, 
		paramListUser, paramListAuto, paramValuesLists)
	
	## all strings for chart on similarity between iterations
	if "graphs" in reuse and isfile(join(graphsDir, 'Param-strings-vs-iter-sim-'+ fileBaseName+'.pdf')):
		vprint("Reusing param-strings-vs-change plot")
	else:
		vprint("Creating param-strings-vs-change plot")
		savedFiles += createBarsAllParametersVsIterSim(
			paramList, statsDict, graphsDir, fileBaseName, verbose, cut)
	
	## parameter strings for bar charts (no linear study possible)
	if fullStringsDict:
		for plotType in ["files", "time", "sim"]:
			if "graphs" in reuse and isfile(join(graphsDir, 
				'Param-strings-vs-'+plotType+'-'+ fileBaseName+'.pdf')):
				vprint("Reusing param-strings-vs-" + plotType + " plot")
			else:
				vprint("Creating param-strings-vs-" + plotType + " plot")
				f = {"files"	: createBarsParametersVsFiles,
					"time"	: createBarsParametersVsTime,
					"sim"	: createBarsParametersVsWNSim}[plotType]
				savedFiles += f(paramListUser, fullStringsDict, graphsDir, fileBaseName, verbose, cut)
		
	## parameter strings for plots (linear study possible)
	if singleParamDicts:
		## each parameter varying...
		for paramKey, paramValues in \
			zip(singleParamDicts.iterkeys(), [vl for vl in paramValuesLists if vl]):
			## ... in each different configuration
			for confKey in singleParamDicts[paramKey].iterkeys():
				## plot for each interesting element
				for plotType in ["files", "time", "sim"]:
					if "graphs" in reuse and isfile(join(graphsDir, 
						paramKey+'-vs-'+plotType+'['+confKey+']'+'-'+ fileBaseName+'.pdf')):
						vprint("Reusing " + paramKey + "-vs-" + \
							plotType + " plot with " + confKey)
					else:
						vprint("Creating " + paramKey + "-vs-" + \
							plotType + " plot with " + confKey)
						
						f = {"files"	: createPlotSingleParamVsFiles,
							"time"	: createPlotSingleParamVsTime,
							"sim"	: createPlotSingleParamVsWNSim}[plotType]
						savedFiles += f(paramValues, singleParamDicts[paramKey][confKey], 
							graphsDir, fileBaseName, paramKey, confKey, verbose, cut)
	return savedFiles
	
	
## Produces dictionaries suitable to generate plots using...
##		1) the main dictionary (statistics for everything in the order in which the parameter strings 
##		were used), 
##		2) three separate lists for all of the parameter strings, those that can be studied in relation to others,
##		and those which are completely independent,
##		3) the lists of values for each of the studied parameters.
## The result is:
##		1) a dictionary for independent parameter strings
##			{ str1 : 	{ stat1 : 	[v1, v2, ...],
##					   stat2 :	... },
##			   str2 :	... }
##		2) a dictionary for the parameter strings studied (sorted following the ascending order of 
## 		parameter values)
##			{ param1 : { conf1 :	{ stat1 : 	[v1, v2, ...],
##					   			   stat2 :	... },
##					     conf2 : ... },
##			   param2 : ... }
##		3) the lists of values, but this time sorted in ascending order
## Therefore, as a result, everything is sorted so that it can be exploited in an efficient manner to 
## plot the statistics.
## @return userDict, singleParamVarDicts, paramValuesLists
def decomposeDictionary(mainDict, allParamStrings, paramStringsUser, paramStringsAuto, 
	paramValuesLists):
	
	## PARAMETER STRINGS SPECIFIED BY THE USER
	userDict = {}
	if paramStringsUser:
		for key in mainDict.iterkeys():
			if len(mainDict[key]) == 1:
				userDict[key] = mainDict[key]
			else :
				userDict[key] = (tuple(value for value, paramStr 
					in zip(mainDict[key], allParamStrings) if paramStr in paramStringsUser))
	
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
								confSettings.append([s for s 
									in [substring1, substring2, substring3] if s != ''])
				else:
					confSettings.append('')
				
				sortedValuesList = sorted([value for value in targetValuesList])
				
				## build a dictionary for each of these configurations
				singleParamDict = {}
				for settings in confSettings:
					confDict = {}
					for key in mainDict.iterkeys():
						if len(mainDict[key]) == 1:
							confDict[key] = mainDict[key]
						else :
							## select fields corresponding to automatically determined strings 
							## containing all of the fixed settings for this conf
							confDict[key] = [value for value, paramStr 
								in zip(mainDict[key], allParamStrings)
								if paramStr in paramStringsAuto and 
								sum([substring in paramStr for substring in settings]) == len(settings)]
							
							## sort dictionary so that it follows ascending order of parameter values
							confDict[key]= [confDict[key][j] 
								for j in [targetValuesList.index(paramValue) 
								for paramValue in sortedValuesList]]
					
					singleParamDict['"'+' '.join(settings)+'"'] = confDict
				singleParamVarDicts[string.upper(targetName)] = singleParamDict
	
	## finally, sort lists of parameter values
	paramValuesLists = [sorted(list) for list in paramValuesLists]
	
	return userDict, singleParamVarDicts, paramValuesLists


## Creates bar chart showing the relation between the parameters used and the result files obtained 
## by running Byblo (size and # of lines)
## No model: the parameters (x values) are not continuous data
## @return names of saved files
def createBarsParametersVsFiles(paramList, statsDictionary, graphsDir, graphName="", 
	verbose=False, cut=False):
	
	## 1) Figure set up
	f, (nbItems) = pl.subplots(1, 1)
	yLabelPos = [-0.1, 0.5] ## A5
	f.set_size_inches(8.3, 5.8)
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2)
	## title
	f.suptitle("Variation of the parameters\nImpact on item counts", 
		fontsize=14, fontweight='bold')
	## info on the file used
	f.text(0.02, 1, "Input file:" +  graphName, horizontalalignment='left', 
		verticalalignment='top', transform = f.transFigure, fontsize=10)
	
	## 2) Data representation
	itemTypes = ["Entries", "Features", "Events", "Thesaurus_Entries"]
	colors = ["RoyalBlue", "red", "gold", "Chartreuse"]
	ind = np.arange(len(paramList))
	width = 0.2
	for i, (item, color) in enumerate(zip(itemTypes, colors)):
		itemCount  = statsDictionary["Number_Of_"+item+"_After_Filtering"]
		nbItems.bar(ind+0.2*(i+1), itemCount, width=width, align="center", color=color, label=item)
	
	nbItems.set_xticks(ind+0.5)
	nbItems.set_xticklabels(paramList, rotation=15, size='small')
	nbItems.set_xlim(xmin=0, xmax=len(paramList))
	decorateGraph(nbItems, "Item counts", "full parameter string", "number of items", 
		yLabelPos, "upper left", largeY=True)
	
	## 3) Save file
	f.savefig(join(graphsDir, "Param-strings-vs-files-" + graphName + '.pdf'))
	pl.close()
	return [join(graphsDir, "Param-strings-vs-files-" + graphName + '.pdf')]


## Creates a bar chart showing the relation between the parameters used and the run time of Byblo
## No model: the parameters (x values) are not continuous data
## @return names of saved files
def createBarsParametersVsTime(paramList, statsDictionary, graphsDir, graphName="", 
	verbose=False, cut=False):
	
	## 1) Figure set up
	f, (time) = pl.subplots(1, 1)
	yLabelPos = [-0.1, 0.5]
	f.set_size_inches(8.3, 5.8) ## A5
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2)
	## title
	f.suptitle("Variation of the parameters\nImpact on running time", 
		fontsize=14, fontweight='bold')
	## info on the file used
	f.text(0.02, 1, "Input file:" +  graphName, horizontalalignment='left', 
		verticalalignment='top', transform = f.transFigure, fontsize=10)
	
	## 2) Data representation
	times, timeUnit = convertTimeRange(statsDictionary["Byblo_Run_Time"])
	ind = np.arange(len(paramList))
	width = 0.5
	time.bar(ind+0.5, times, width=width, align="center", color="aquamarine", 
		label="run time function of parameters")
	
	time.set_xticks(ind+0.5)
	time.set_xticklabels(paramList, rotation=15, size='small')
	time.set_xlim(xmin=0, xmax=len(paramList))
	decorateGraph(time, "Running time", "full parameter string", 
		"time ("+timeUnit+")", yLabelPos, largeY=True)

	## 3) Save file
	f.savefig(join(graphsDir, "Param-strings-vs-time-" + graphName + '.pdf'))
	pl.close()
	return [join(graphsDir, "Param-strings-vs-time-" + graphName + '.pdf')]


## Creates a bar chart showing the relation between the parameters used and the similarity of 
## the resultant thesaurus with WordNet
## No model: the parameters (x values) are not continuous data
## @return names of saved files
def createBarsParametersVsWNSim(paramList, statsDictionary, graphsDir, graphName="", 
	verbose=False, cut=False):
	
	## 1) Figure set up
	f, (sim) = pl.subplots(1, 1)
	yLabelPos = [-0.1, 0.5]
	f.set_size_inches(8.3, 5.8) ## A5
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2)
	## title
	f.suptitle("Variation of the parameters\nImpact on similarity with WordNet", 
		fontsize=14, fontweight='bold')
	## info on the file used
	f.text(0.02, 1, "Input file:" +  graphName, horizontalalignment='left', 
		verticalalignment='top', transform = f.transFigure, fontsize=10)
	
	## 2) Data representation
	simWithWn= statsDictionary["Similarity_With_WordNet"]
	ind = np.arange(len(paramList))
	width = 0.5
	sim.bar(ind+0.5, simWithWn, width=width, align="center", color="SpringGreen", 
		label="sim function of parameters")
	
	sim.set_xticks(ind+0.5)
	sim.set_xticklabels(paramList, rotation=15, size='small')
	sim.set_xlim(xmin=0, xmax=len(paramList))
	decorateGraph(sim, "Similarity with WordNet", "full parameter string", 
		"similarity score", yLabelPos)

	## 3) Save file
	f.savefig(join(graphsDir, "Param-strings-vs-WN-sim-" + graphName + '.pdf'))
	pl.close()
	return [join(graphsDir, "Param-strings-vs-WN-sim-" + graphName + '.pdf')]


## Creates a abar chart showing the relation between the parameters used and the similarity there is
## between one thesaurus and the next one (next iteration)
## No model: the parameters don't follow any rule and appear in simple chronological order
def createBarsAllParametersVsIterSim(paramList, statsDictionary, graphsDir, 
	graphName="", verbose=False, cut=False):
	
	
	## 1) Figure set up
	f, (sim) = pl.subplots(1, 1)
	yLabelPos = [-0.1, 0.5]
	f.set_size_inches(8.3, 5.8) ## A5
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2)
	## title
	f.suptitle("Variation of the parameters\nImpact on similarity between iterations", 
		fontsize=14, fontweight='bold')
	## info on the file used
	f.text(0.02, 1, "Input file:" +  graphName, horizontalalignment='left', 
		verticalalignment='top', transform = f.transFigure, fontsize=10)
	
	## 2) Data representation
	simWithPrev = statsDictionary["Similarity_With_Previous"][1:] + [0]
	ind = np.arange(len(paramList))
	width = 0.5
	sim.bar(ind+0.5, simWithPrev, width=width, align="center", color="goldenrod", 
		label="sim between iterations")
	
	sim.set_xticks(ind)
	sim.set_xticklabels(paramList, rotation=15, size='small')
	sim.set_xlim(xmin=0, xmax=len(paramList)-1)
	decorateGraph(sim, "Similarity between constructed thesauri", 
		"full parameter string", "similarity score", yLabelPos)

	## 3) Save file
	f.savefig(join(graphsDir, "Param-strings-vs-iter-sim-" + graphName + '.pdf'))
	pl.close()
	return [join(graphsDir, "Param-strings-vs-iter-sim-" + graphName + '.pdf')]


## Creates plots showing the relation between the parameters used and the result files obtained 
## by running Byblo (# of lines)
## Model: none (too many counts)
## @return names of saved files
def createPlotSingleParamVsFiles(paramList, statsDictionary, graphsDir, graphName="",
	paramName="", confName="", verbose=False, cut=False):
	
	## 1) Figure set up
	yLabelPos = [-0.1, 0.5]
	fileInfo, confInfo = "Input file:" +  graphName, "\nFixed settings:" +  confName
	f, (nbItems) = pl.subplots(1, 1)
	##size and margins
	f.set_size_inches(8.3, 5.8) ## A5
	f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2)
	## title
	f.suptitle("Variation of the " + fullParamName(paramName) +\
		"\nImpact on item counts", fontsize=14, fontweight='bold')
	## info on the file and configuration used
	f.text(0.02, 1, fileInfo, horizontalalignment='left', verticalalignment='top',
		transform = f.transFigure, fontsize=10)
	f.text(0.98, 1, confInfo, horizontalalignment='right', verticalalignment='top',
		transform = f.transFigure, fontsize=10)
	
	## 2) Data representation
	itemTypes = ["Entries", "Features", "Events", "Thesaurus_Entries"]
	colors = ["RoyalBlue", "red", "gold", "Chartreuse"]
	for item, color in zip(itemTypes, colors):
		itemCount  = statsDictionary["Number_Of_"+item+"_After_Filtering"]
		nbItems.plot(paramList, itemCount, color=color, linestyle="dashed", marker='o', 
			label=string.replace(string.lower(item), '_', ' '), alpha=0.6)

	decorateGraph(nbItems, "Item counts", "filter threshold value", "number of items", 
		yLabelPos, "upper left", largeX=True, largeY=True)
	
	## 3) Save file
	fileName = paramName + "-vs-counts["+confName+"]-" + graphName + ".pdf"
	f.savefig(join(graphsDir, fileName))
	pl.close()
	return join(graphsDir, fileName)


## Creates a plot showing the relation between the variation of ONE chosen parameter and
## the run time of Byblo
## Model: polynomial
## @return names of saved files
def createPlotSingleParamVsTime(paramList, statsDictionary, graphsDir, graphName="",
	paramName="", confName="", verbose=False, cut=False):
	
	## 1) Figure set up
	yLabelPos = [-0.1, 0.5]
	fileInfo, confInfo = "Input file:" +  graphName, "\nFixed settings:" +  confName
	if not cut:
		f, (time, fitTime) = pl.subplots(2, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.4) ## margins
		## title
		f.suptitle("Variation of the " + fullParamName(paramName) +\
			"\nImpact on running time", fontsize=14, fontweight='bold')
		## info on the file and configuration used
		f.text(0.02, 1, fileInfo, horizontalalignment='left', verticalalignment='top',
			transform = f.transFigure, fontsize=10)
		f.text(0.98, 1, confInfo, horizontalalignment='right', verticalalignment='top',
			transform = f.transFigure, fontsize=10)
	else:
		individualPlots = [pl.subplots(1, 1) for x in xrange(2)]
		figures = [p[0] for p in individualPlots]
		time, fitTime = [p[1] for p in individualPlots]
		for f in figures:
			## size and margins
			f.set_size_inches(8.3, 5.8)
			f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2)
			## title
			f.suptitle("Variation of the " + fullParamName(paramName) +\
				"\nImpact on running time", fontsize=14, fontweight='bold')
			## info on the file and configuration used
			f.text(0.02, 1, fileInfo, horizontalalignment='left', verticalalignment='top',
				transform = f.transFigure, fontsize=10)
			f.text(0.98, 1, confInfo, horizontalalignment='right', verticalalignment='top',
				transform = f.transFigure, fontsize=10)	
	
	## 2) Data representation
	times,  timeUnit = convertTimeRange(statsDictionary["Byblo_Run_Time"])
	time.plot(paramList, times, color="firebrick", marker='o', linestyle='dashed',
		label="run time function of "+paramName, alpha=0.6)
		
	decorateGraph(time, "Running time", "filter threshold value", 
		"time ("+timeUnit+")", yLabelPos, largeY=True)
	
	## 3) Model fitting
	fitTime.fill_between(paramList, times, color="lightgrey")
	## fit functions (if enough values)
	if len(paramList) >= 2:
		method = 6
		mFit, mLabel, mColor= fittingMethod(paramList, times, method, verbose=verbose)
		mLabel = "polynomial model proposed"
		fitTime.plot(paramList, mFit, label=mLabel, color=mColor)

	fitTime.set_xlim(xmax=paramList[-1])
	decorateGraph(fitTime, "Polynomial model", "filter threshold value", 
		"time ("+timeUnit+")", yLabelPos, "upper right")
	
	## 4) Save file(s)
	savedFiles = []
	if not cut:
		fileName = paramName + "-vs-time["+confName+"]-" + graphName + ".pdf"
		f.savefig(join(graphsDir, fileName))
		savedFiles.append(join(graphsDir, fileName))
	else:
		for i, f in enumerate([p[0] for p in individualPlots]):
			fileName = paramName + "-vs-time["+confName+"]-" + graphName +\
				'-' + str(i+1) + ".pdf"
			f.savefig(join(graphsDir, fileName))
			savedFiles.append(join(graphsDir, fileName))
	pl.close()
	return savedFiles


## Creates a plot showing the relation between the variation of ONE chosen parameter and similarity of 
## the resultant thesaurus with WordNet
## Model: polynomial
## @return names of saved files
def createPlotSingleParamVsWNSim(paramList, statsDictionary, graphsDir, graphName="",
	paramName="", confName="", verbose=False, cut=False):
	
	## 1) Figure set up
	yLabelPos = [-0.1, 0.5]
	fileInfo, confInfo = "Input file:" +  graphName, "\nFixed settings:" +  confName
	if not cut:
		f, (sim, fitSim) = pl.subplots(2, 1)
		f.set_size_inches(8.3, 11.7) ## set figure size to A4
		f.subplots_adjust(left=0.15, right=0.85, wspace=None, hspace=0.4) ## margins
		## title
		f.suptitle("Variation of the " + fullParamName(paramName) +\
			"\nImpact on similarity with WordNet", fontsize=14, fontweight='bold')
		## info on the file and configuration used
		f.text(0.02, 1, fileInfo, horizontalalignment='left', verticalalignment='top',
			transform = f.transFigure, fontsize=10)
		f.text(0.98, 1, confInfo, horizontalalignment='right', verticalalignment='top',
			transform = f.transFigure, fontsize=10)
	else:
		individualPlots = [pl.subplots(1, 1) for x in xrange(2)]
		figures = [p[0] for p in individualPlots]
		sim, fitSim = [p[1] for p in individualPlots]
		for f in figures:
			## size and margins
			f.set_size_inches(8.3, 5.8)
			f.subplots_adjust(left=0.15, right=0.85, wspace=None, top=0.8, bottom=0.2)
			## title
			f.suptitle("Variation of the " + fullParamName(paramName) +\
				"\nImpact on similarity with WordNet", fontsize=14, fontweight='bold')
			## info on the file and configuration used
			f.text(0.02, 1, fileInfo, horizontalalignment='left', verticalalignment='top',
				transform = f.transFigure, fontsize=10)
			f.text(0.98, 1, confInfo, horizontalalignment='right', verticalalignment='top',
				transform = f.transFigure, fontsize=10)	
	
	## 2) Data representation
	simWithWn = statsDictionary["Similarity_With_WordNet"]
	sim.plot(paramList, simWithWn, color="firebrick", marker='o', linestyle="dashed",
		label="sim function of "+paramName, alpha=0.6)
		
	decorateGraph(sim, "Similarity with WordNet", "filter threshold value", 
		"similarity score", yLabelPos)
	
	## 3) Model fitting
	fitSim.fill_between(paramList, simWithWn, color="lightgrey")
	## fit functions (if enough values)
	if len(paramList) >= 2:
		method = 6
		mFit, mLabel, mColor= fittingMethod(paramList, simWithWn, method, verbose=verbose)
		mLabel = "polynomial model proposed"
		fitSim.plot(paramList, mFit, label=mLabel, color=mColor)

	fitSim.set_xlim(xmax=paramList[-1])
	decorateGraph(fitSim, "Polynomial model", "filter threshold value", 
		"frequency", yLabelPos, "upper right")
	
	## 4) Save file(s)
	savedFiles = []
	if not cut:
		fileName = paramName + "-vs-WN-sim["+confName+"]-" + graphName + ".pdf"
		f.savefig(join(graphsDir, fileName))
		savedFiles.append(join(graphsDir, fileName))
	else:
		for i, f in enumerate([p[0] for p in individualPlots]):
			fileName = paramName + "-vs-WN-sim["+confName+"]-" + graphName +\
				'-' + str(i+1) + ".pdf"
			f.savefig(join(graphsDir, fileName))
			savedFiles.append(join(graphsDir, fileName))
	pl.close()
	return savedFiles
	
	
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
#~ def deleteOnExit(deleteList, outputDirectory, sampleFiles=[], originalInputFile="", verbose=False):
	#~ if(deleteList != ["nothing"]):
		#~ for directory in deleteList:
			#~ os.system("rm -r " + join(outputDirectory, directory))
			#~ if verbose:
				#~ print "   Deleted directory " + directory


