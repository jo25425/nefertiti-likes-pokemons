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
__credits__ = ["Joanne Robert"]
__license__ = "3-clause BSD"
__version__ = "1.0"
__maintainer__ = "Joanne Robert"
__email__ = "jr317@sussex.ac.uk"
__status__ = "Development"
#~ -----------------------------------------------------------------------------

import argparse
import time
import sys
import subprocess
import os
from os.path import *
from menusystem import *
import inputchecking
import outputformatting
import bybloeval
import cmpstats

## global variable holding extensions for output files of Byblo
BYBLO_OUTPUT_EXTENSIONS= [	".entries", ".entries.filtered", 
						".entry-index.d.0", ".entry-index.i.0",
						".enumerated",
						".events", ".events.filtered",
						".feature-index.d.0", ".feature-index.i.0",
						".features", ".features.filtered",
						".sims", ".sims.neighbours", ".sims.neighbours.strings"
					]

class BybloCmp:
	
	"""
	Compare methods for the construction of a disitributional thesaurus
	"""
	
	## Initialises the parameters of the module
	def __init__(self, inputFile, bybloDir, storageDir, outputFile, verbose):
		
		## input file for pre-processed data
		self.inputFile = abspath(inputFile)
		## location of the Byblo directory
		self.bybloDir = abspath(bybloDir)
		## storage directory for Byblo output
		self.storageDir = abspath(storageDir)
		## output file for  comparison results 
		self.outputFile = abspath(outputFile) if not isdir(outputFile)\
			else join(abspath(outputFile), "results.cmp")
		## verbose option
		self.verbose = verbose
		## sequence menu
		self.iterationEndMenu = self.initIterationEndMenu()
		## options menu
		self.optionsMenu = self.initOptionsMenu()
		## main menu
		self.mainMenu = self.initMainMenu()
		## record of iterations
		self.record = []
		## parser for input
		self.parser = inputchecking.Parser()
		## printer for output
		self.printer = outputformatting.Printer(self.verbose)
	
	## Run
	def run(self):
		## start operations.
		stime = time.time()
		self.printer.mainTitle("Byblo Configuration Helper")
		
		## direct operations from main menu
		self.mainMenu.waitForInput()
		
		etime = time.time()
		self.printer.info("Execution took "+str(etime-stime)+"seconds")
		
	
	## Creates the main menu.
	## It allows to navigate in the program, accessing all available functionalities."
	## @return menu
	def initMainMenu(self):
		mainChoices = [
			Choice(1, description="Execution of Byblo iterations",
				## execution function
				handler=self.execution),
			Choice(2, description="Iterations planning"),
			Choice(3, description="Comparison results"),
			Choice(4, description="Help (about Byblo usage)",
				## function displaying Byblo help
				handler=self.displayBybloHelp),
			Choice(5, description="Options", 
				## directly access the OPTIONS submenu
				subMenu=self.optionsMenu),
			Choice(6, description="Exit", 
				## function returning a 'False' to exit this menu
				handler=lambda anything:False)
		]
		return Menu("Main Menu", mainChoices, "What do you want to do?")
		
		
	## Creates the options menu.
	## It allows to change all of the settings that were defined by the command-line.
	## @return menu
	def initOptionsMenu(self):
		## temporary checking functions
		isFile = lambda str: True
		isDir = lambda str:True
		
		## Displays the current  value of an option.
		formatVal =  lambda str: "\n\t\t[" + str + "]"
		
		optionsChoices = [
			Choice(1, description="Input file for preprocessed data" + formatVal(self.inputFile),
				value=[isFile,1,  "inputFile"], handler=self.changeOption),
			Choice(2, description="Location of the Byblo directory" + formatVal(self.bybloDir),
				value=[isDir, 2, "bybloDir"], handler=self.changeOption),
			Choice(3, description="Storage directory for Byblo output" + formatVal(self.storageDir),
				value=[isDir, 3, "storageDir"], handler=self.changeOption),
			Choice(4, description="Output file for comparison results" + formatVal(self.outputFile),
				value=[isFile, 4, "outputFile"], handler=self.changeOption),
			Choice(5, description="Back to main menu", handler=lambda anything:False)
		]
		return Menu("Options Menu", optionsChoices, "What do you want to change?")
	
	
	## Creates the iteration end menu.
	## It allows to choose, at the end of an iteration, whether to continue working in the same sequence
	## or to end it and consequently plot the appropriate statistics over this entire sequence.
	## (Sequence = series of iterations sharing common non-studied parameters)
	## @return menu
	def initIterationEndMenu(self):
		sequenceChoices = [
			Choice(1, description="Continue working in the same sequence of iterations",
				handler=lambda anything:False),
			Choice(2, description="End this sequence and plot some statistics over this entire sequence",
				handler=self.plotSequence)
		]
		return Menu("Iteration End Menu", sequenceChoices, "Iteration complete. Now...")
	
	
	## Changes one of the options.
	## Only produces correct values that are accepted by external verifiers (type/file/parameters)
	def changeOption(self, args):
		checkingFunction, selector, attr = args # separate function and attribute name
		oldVal = getattr(self, attr) # store old value
		choiceIndex = selector-1 # determine option index in list of choices
			
		while True:
			newVal = raw_input("\nNew value? ")
			if checkingFunction(newVal):
				## modify corresponding option
				setattr(self, attr, newVal)
				
				## update options menu
				setattr(self.optionsMenu.choices[choiceIndex], "description", 
					self.optionsMenu.choices[choiceIndex].description.replace(oldVal, newVal))
				break
			else:
				print "Invalid value."
	
	
	## Displays Byblo help to give information about the possible parameter strings, how to build them
	## and what they mean
	def displayBybloHelp(self, args):
		## move to Byblo directory, and run it in a subprocess
		startDir=abspath(os.getcwd())
		os.chdir(self.bybloDir)
		self.printer.info("Moved to " + os.getcwd())
		out = subprocess.call(abspath("./byblo.sh ") + " --help", shell = True)
		## fail?
		if(out != 0 and out != 255):
			self.printer.info("Byblo help failed.\n       Fail Code: " + str(out))
		os.chdir(startDir)
		self.printer.info("Moved back to " + os.getcwd())
		
		return True # to stay in this level's menu
	
	
	## Executes Byblo once, going through the following steps:
	## 1. determine the parameters to use (checked by parser)
	## 2. run Byblo with these parameters
	## 3. evaluate the resultant thesaurus against the one from the previous iteration
	## 4. generate histograms showing the distributions observed (of entries, features, events, and
	##     similarity scores)
	## FInally, the user is offered the possibility to end the current sequence.
	## @return  True
	def execution(self, args):
		## determine parameters to use
		self.printer.stage(1, 5, "Determining parameters")
		"""
			for now directly, later from planned stuff... or not? (could fairly easily just read from a file)
		"""
		while True:
			parameters = raw_input("Parameters for this iteration of Byblo? ")
			if self.parser.checkBybloSettings(parameters, type='studied'): break
			else: print ""
		
		## run Byblo
		self.printer.stage(2, 5, "Running Byblo")
		returnCode = self.runByblo(parameters)
		if returnCode != 0: return True # if Byblo failed, end iteration
		
		## evaluate built thesaurus
		self.printer.stage(3, 5, "Comparing with previous iteration")
		self.evalIteration()
		
		## generate histograms
		self.printer.stage(4, 5, "Generating histograms")
		self.plotIteration()
		
		## write to output file
		self.printer.stage(5, 5, "Writing iteration summary")
		self.writeIteration()
		
		## either continue or end current sequence and plot
		self.iterationEndMenu.waitForInput()
		
		return True # to stay in this level's menu
		
		
	## Runs Byblo in a subprocess, renames files produced so that they reflect the settings used,
	## and consequently adds an entry to the record of iterations
	## @return code indicating success/failure of the iteration
	def runByblo(self, parameters):
		## create output directory when required
		if not exists(self.storageDir):
			os.makedirs(self.storageDir)
			
		## create thesauri directory when required
		thesauriDir = join (self.storageDir, "thesauri")
		if not exists(thesauriDir):
			os.makedirs(thesauriDir)
		
		## move to Byblo directory
		startDir=abspath(os.getcwd())
		os.chdir(self.bybloDir)
		self.printer.info("Moved to " + os.getcwd())
		
		## execute Byblo in a subprocess
		logFile = open(devnull, 'w') if not self.verbose else None
		stime = time.time()
		out = subprocess.call(abspath("./byblo.sh ") + " -i " + self.inputFile + " -o " + thesauriDir +\
			" "+ parameters, shell=True, stdout=logFile, stderr=logFile)
		etime = time.time()
		if logFile: logFile.close()
		
		## move back to initial directory
		os.chdir(startDir)
		self.printer.info("Moved back to " + os.getcwd())
		
		## stop here in case of fail
		if out != 0:
			self.printer.info("Byblo failed with settings: " + parameters + "\n       Fail Code: " + str(out))
			return -1
		
		## rename files to include the parameter string in their name
		base = join(thesauriDir, basename(self.inputFile))
		for ext in BYBLO_OUTPUT_EXTENSIONS:
			filename = base + ext
			if exists(filename):
				newname = base + cmpstats.paramSubstring(parameters) + ext
				os.rename(filename, newname)
		
		## update record
		aboutIteration = {
			"input" 		: self.inputFile,
			"settings" 	: parameters if parameters else 'None',
			"output" 		: self.storageDir,
			"runtime"	: etime-stime,
			"sim-with-WN"	: 'None',
			"change-from-prev": 'None'
			}
		self.record.append(aboutIteration)
		return 0


	## Evaluates the result of this last iteration of Byblo against...
	## 1) WordNet, in order to have a fixed reference for comparison
	## 2) the result of the previous iteration, in order to "measure" change
	def evalIteration(self):
		findThesaurus = lambda rec: join(rec["output"], "thesauri",
					basename(rec["input"])\
					+ cmpstats.paramSubstring(rec["settings"]
						if rec['settings'] != 'None' else "") 
					+ ".sims.neighbours.strings")
		evalOutput = "eval-iter"+str(len(self.record))+".tmp"
		
		
		
		## compute similarity with WordNet
		with open(findThesaurus(self.record[-1]), 'r') as currentTh:
			evalTask1 = bybloeval.BybloEval([currentTh], evalOutput, method="rank", 
				testIndex=0, maxRank=None, maxIndex=None, verbose=self.verbose)
			evalTask1.run()
			
			## write evaluation result to record
			self.record[-1]["sim-with-WN"] =  \
				100.0 * float(open(evalOutput, 'r').readlines()[0].split('\t')[0])
			
		## compute similarity with previously built thesaurus
		if len(self.record) < 2:
			self.printer.info("Comparative evaluation only possible from the second iteration. Skipped.")
		else: 
			with open(findThesaurus(self.record[-2]), 'r') as previousTh,\
				open(findThesaurus(self.record[-1]), 'r') as currentTh:
				evalTask2 = bybloeval.BybloEval([currentTh, previousTh], evalOutput, method="rank", 
					testIndex=0, maxRank=None, maxIndex=None, verbose=self.verbose)
				evalTask2.run()
				
				## write evaluation result to record
				self.record[-1]["change-from-prev"] = \
					100.0 * (1.0 - float(open(evalOutput, 'r').readlines()[0].split('\t')[0]))
				
		## delete temporary file
		os.remove(evalOutput)
		self.printer.lines(self.record[-1], title="Record for the current interation")
	
	
	## Plots some statistics over a single iteration of Byblo using the cmpstats module.
	##
	## The graphs generated are histograms representing the distributions of entries, features,
	## and events in the corpus, as well as the distribution of similarity scores in the resultant
	## thesaurus.
	##
	## The histograms will be put in the subdirectory "graphs" of the main output directory.
	def plotIteration(self):
		## retrieve settings strings
		findSettings = lambda dict: dict["settings"] if dict["settings"] != 'None' else ''
		paramString = findSettings(self.record[-1])
		## generate histograms
		self.record[-1]["histograms"] = cmpstats.generateHistograms(
			[self.inputFile], 	# sampleFileNames
			[paramString],	# paramList
			self.storageDir, 	# outputDir
			self.bybloDir, 		# bybloDir
			[], 				# reuse
			False, 			# verbose
			False)			# cut
			
		self.printer.info("Histograms created in subdirectory \"graphs\" of main output directory.")
	
	
	##
	##
	def writeIteration(self):
		findThesaurus = lambda rec: join(rec["output"], "thesauri",
					basename(rec["input"])\
					+ cmpstats.paramSubstring(rec["settings"]
						if rec["settings"] != "None" else "") 
					+ ".sims.neighbours.strings")
		
		with open(self.outputFile, ('w' if not isfile(self.outputFile) else 'a')) as output:
			## sequence start if needed
			out = lambda s : output.write(s + '\n')
			if len(self.record) == 1:
				out( ''.join(['-' for i in range(100)]))
				out("%-20s" %  "Local current time: " +
					time.asctime( time.localtime(time.time()) ))
				out("%-20s" % "Data file: " + self.inputFile)
				out("%-20s" % "Storage directory: " + self.storageDir)
				out('') # empty line
			
			## about last iteration
			out = lambda s : output.write('\t' + s + '\n')
			title = "Iteration " + str(len(self.record))
			out(title)
			out(''.join(['-' for i in range(len(title))]))
			out("%-25s" % "Parameters used: " + self.record[-1]["settings"])
			out("%-25s" % "Running time: " + 
				"%.2f" % self.record[-1]["runtime"] + " seconds")
			out("%-25s" % "Similarity with WordNet: " +
				"%.2f" % self.record[-1]["sim-with-WN"] + " %")
			out("%-25s" % "Change from previous: " + 
				(("%.2f" % self.record[-1]["change-from-prev"] + " %") 
				if type(self.record[-1]["change-from-prev"]) == float
				else str(self.record[-1]["change-from-prev"])))
			out("%-25s" % "Resultant thesaurus: " + findThesaurus(self.record[-1]))
			
			out("Corresponding graphs:")
			for histogram in self.record[-1]["histograms"]:
				out(' * ' + histogram)
			out('') # empty line
	
	
	## Plots some statistics over the complete iteration sequence using the cmpstats module.
	##
	## The graphs generated plot the parameters used...
	## 		- either individually when they take on several values, e.g. -fef = {10, 20, 50},
	## 		- or as a whole string when they don't relate, e.g. "-Smn 0.5" and "-fvf 10"
	## ...against one of the characteristics of the Byblo run and result, namely
	##		- the running time,
	##		- the number of entries, features, events kept,
	##		- the similarity with WordNet.
	## There is also an additional graph showing the amount of change from one iteration to the next one.
	##
	## The graphs will be put in the subdirectory "graphs" of the main output directory.
	## @return False
	def plotSequence(self, args=None):
		self.printer.info("Computing stats over the entire sequence")
		
		## build a dictionary containing all the wanted statistics
		statsDict = self.recordToDictionary()
		self.printer.lines( statsDict )
		allStrings, stringsForCharts, stringsForPlots, paramDict = self.sortParameterStrings()
		
		## generate plots
		self.record[-1]["plots"] = cmpstats.generatePlots(
			statsDict,			# statsDict
			allStrings,			# paramList 
			stringsForCharts,		# paramListUser 
			stringsForPlots,		# paramListAuto
			paramDict,			# paramValuesLists
			self.inputFile,			# fileName
			self.storageDir,		# outputDir
								# reuse=[]
			verbose=True)		# verbose=False
								# cut=False
	
		self.printer.info("Plots created in subdirectory \"graphs\" of main output directory.")
		self.writeSequenceEnd()
		
		return False # to go back to main menu (after question)
	
	
	## Produces from the record a dictionary containing some statistics for the sequence ending:
	## 1) from the input: number of distinct entries and events, average number of events per entry
	## 2) from the output: number of filtered entries, features, events and neighbour sets, running 
	## time of Byblo, results of the evaluation against WordNet and the previous iteration
	## @return dictionary from the record
	def recordToDictionary(self):
		findSettings = lambda dict: dict["settings"] if dict["settings"] != 'None' else ''
		countLines = lambda file: sum([1 for entry in open(file)])
		statsDict = {}
		
		## compute statistics regarding input (one per sequence)
		nbEntries, nbEvents = 0, 0
		for entry in open(self.inputFile, 'r'):
			nbEntries += 1
			nbEvents += len(entry.split()) -1
		
		## add them to thr dictionary
		statsDict["Total_Number_Of_Distinct_Entries"] = [nbEntries]
		statsDict["Total_Number_Of_Distinct_Events"] = [nbEvents]
		statsDict["Average_Number_Of_Events_By_Entry"] = [1.0 * nbEntries / nbEvents]
		
		## compute statistics regarding output (once per iteration)
		nbFilteredEntries, nbFilteredFeatures, nbFilteredEvents, nbFilteredNeighbourSets = [], [], [], []
		runTime, simWithWN, changeFromPrev = [], [], []
		for iteration in self.record:
			startFileName = join(self.storageDir, "thesauri", basename(self.inputFile))\
				+ cmpstats.paramSubstring(findSettings(iteration))
				
			nbFilteredEntries.append( countLines(startFileName + ".entries.filtered") )
			nbFilteredFeatures.append( countLines(startFileName + ".features.filtered") )
			nbFilteredEvents.append( countLines(startFileName + ".events.filtered") )
			nbFilteredNeighbourSets.append( countLines(startFileName + ".sims.neighbours") )
			runTime.append( iteration["runtime"] )
			simWithWN.append( iteration["sim-with-WN"] )
			changeFromPrev.append( iteration["change-from-prev"] )
		
		## add them to thr dictionary
		statsDict["Number_Of_Lines_In_File_.entries.filtered"] = nbFilteredEntries
		statsDict["Number_Of_Lines_In_File_.features.filtered"] = nbFilteredFeatures
		statsDict["Number_Of_Lines_In_File_.events.filtered"] = nbFilteredEvents
		statsDict["Number_Of_Lines_In_File_.sims.neighbours"] = nbFilteredNeighbourSets
		statsDict["Byblo_Run_Time"] = runTime
		statsDict["Similarity_With_WordNet"] = simWithWN
		statsDict["Change_From_Previous"] = changeFromPrev

		return statsDict
	
	
	## Analyses parameter strings used throughout this sequence of Byblo iteration to...
	## 1) Produce, for each studied parameter (fef, fff, fvf, Smn), a list of the distinct values it takes 
	## (if there are several of them)
	## 2) Based on these lists, determine wether each string can be studied in a linear way (line plot), 
	## or will have to be studied separately (bar chart)
	## @return allStrings, stringsForCharts, stringsForPlots, parameterValuesLists
	def sortParameterStrings(self):
		findSettings = lambda dict: dict["settings"] if dict["settings"] != 'None' else ''
		
		## initialise lists and dictionary
		allStrings = map(findSettings, self.record)
		stringsForPlots, stringsForCharts = [], []
		paramNames = ["fef", "fff", "fvf", "Smn"]
		paramDict = {}
		for p in paramNames: paramDict[p] = []
		
		## create lists of distinct values used for each parameter
		for settings in allStrings:
			things = settings.split(' ')
			## look for each parameter in the string...
			for param in [p for p in paramNames if p in settings]:
				## ... and add its value to the right list (if it isn't already there)
				value = things[ things.index('-'+param) + 1 ]
				if value not in paramDict[param]: paramDict[param].append(value)
		
		## create lists of values in this dictionary# and sort strings
		for settings in allStrings:
			## look for parameters that are 1) in the string, and 2) varying (several values)
			if [p for p in paramNames if p in settings and len(paramDict[p]) > 1]: 
				stringsForPlots.append(settings)
			else: 
				stringsForCharts.append(settings)
		
		## lists of values for parameters that will be studied "several values"
		paramValuesLists = [paramDict[p] if len(paramDict[p]) > 1 else [] for p in paramNames]
		return allStrings, stringsForCharts, stringsForPlots, paramValuesLists
	
	
	##
	##
	def writeSequenceEnd(self):
		with open(self.outputFile, 'a') as output:
			out = lambda s : output.write(s + '\n')
			out("Corresponding graphs:")
			for plot in self.record[-1]["plots"]:
				out(' * ' + plot)
			out('') # empty line
	
## Parses a command-line.
if __name__=='__main__':
	argParser = argparse.ArgumentParser(description='Compare different parameterisations for Byblo.')
	
	## input file for preprocessed data
	argParser.add_argument(metavar='file', dest='inputFile', action='store', 
		help='input file for preprocessed data')
	## location of the Byblo directory
	argParser.add_argument('-b', '--byblo', metavar='dir', dest='bybloDir',
		action='store', default="../Byblo-2.1.0",
		help='directory for the location Byblo (default: "../Byblo-2.1.0")')
	## storage directory for Byblo output
	argParser.add_argument('-s', '--storage', metavar='dir', dest='storageDir',
		action='store', default="./byblo-output",
		help='storage directory for Byblo output (default: "./byblo-output")')
	## output file for comparison results
	argParser.add_argument('-o', '--output', metavar='file', dest='outputFile',
		action='store', default="./results.cmp",
		help='output file for comparison results (default: "./results.cmp")')
	## verbose option
	argParser.add_argument('-v', '--verbose', dest='verbose', 
		action='store_true', default=False,
		help='display information about operations (default: False)')
		
	a = argParser.parse_args()
	bybloCmp = BybloCmp(a.inputFile, a.bybloDir, a.storageDir, a.outputFile, a.verbose)
	bybloCmp.run()
	

	