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
import byblostats


class BybloCmp:
	
	"""
	Compare methods for the construction of a disitrbutional thesaurus
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
		"""
		This will eventually be modified so that iterations can be launched using scripts.
		1) a script used to run Byblo each time, but menu still used
		2) a general script that includes a configuration file holding all of the information for iteration planning
		At least ideally....
		"""
		etime = time.time()
		self.printer.info( "Execution took "+str(etime-stime)+"seconds")
		
	
	## Creates the main menu.
	## It allows to navigate in the program, accessing all available functionalities."
	def initMainMenu(self):
		mainChoices = [
			Choice(1, description="Execution of Byblo iterations",
				## execution function
				handler=self.execution),
			Choice(2, description="Iterations planning"),
			Choice(3, description="Comparison results"),
			Choice(4, description="Options", 
				## directly access the OPTIONS submenu
				subMenu=self.optionsMenu),
			Choice(5, description="Exit", 
				## function returing a 'False' to exit this menu
				handler=lambda anything:False)
		]
		return Menu("Main Menu", mainChoices, "What do you want to do?")
		
		
	## Creates the options menu.
	## It allows to change all of the settings that were defined by the command-line.
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
	
	
	def execution(self, args):
		
		## determine parameters to use
		self.printer.stage(1, 4, "Determining parameters")
		# for now directly, later from planned stuff
		while True:
			parameters = raw_input("Parameters for this iteration of Byblo? ")
			if self.parser.checkBybloSettings(parameters, type='studied'): break
			else: print ""
		
		## run Byblo
		self.printer.stage(2, 4, "Running Byblo")
		self.runByblo(parameters)
		
		## evaluate built thesaurus
		self.printer.stage(3, 4, "Comparing with previous iteration")
		self.evalIteration()
		
		## generate histograms
		self.printer.stage(4, 4, "Generating histograms")
		self.plotIteration()
		
		print "Woooooh... Not yet, let's all calm down."
		""" NEED TO:
			X 1) compute similarity between resultant thesaurus and gold-standard
			   2) produce some graphs (distributions only)
		"""
		
		## results over an entire sequence
		#~ self.printer.stage(5, 5, "Computing stats over the entire sequence")
		"""
		Compute rest of the graphs (those working on several outputs)
		"""
		
		return True # to stay in this level's menu
		
		
	## Runs Byblo in a subprocess and consequently adds an entry to the record of iterations 
	def runByblo(self, parameters):
		## create output directory when required
		if not exists(self.storageDir):
			os.makedirs(self.storageDir)
		
		## move to Byblo directory
		startDir=abspath(os.getcwd())
		os.chdir(self.bybloDir)
		self.printer.info("moved to " + os.getcwd())
		
		## execute Byblo in a subprocess
		logFile = open(devnull, 'w') if not self.verbose else None
		stime = time.time()
		out = subprocess.call(abspath("./byblo.sh ") + " -i " + self.inputFile + " -o " + self.storageDir +\
			" "+ parameters, shell=True, stdout=logFile, stderr=logFile)
		etime = time.time()
		if logFile != None:
			logFile.close()
		if out != 0:
			self.printer.info("Byblo failed with settings:\n" + parameters + "\n   Fail Code: " + str(out))
			
		## move back to initial directory
		os.chdir(startDir)
		self.printer.info("moved back to " + os.getcwd())
		
		## rename files to include the parameter string in their name
		EXTENSIONS= [	".entries", ".entries.filtered", 
						".entry-index.d.0", "entry-index.i.0",
						".enumerated",
						".events", ".events.filtered",
						".feature-index.d.0", "feature-index.i.0",
						".features", ".features.filtered",
						".sims", ".sims.neighbours", ".sims.neighbours.strings"
					]
		#~ for file in [self.inputFile + ext for ext in EXTENSIONS]:
			#~ os.
		
		## update record
		aboutIteration = {
			"input" 		: self.inputFile,
			"settings" 	: parameters,
			"output" 		: self.storageDir,
			"runtime"	: etime-stime,
			"eval-vs-WN"	: '',
			"eval-vs-prev": ''
			}
		self.record.append(aboutIteration)
	


	## Evaluates the result of this last iteration of Byblo against...
	## 1) WordNet, in order to have a fixed reference for comparison
	## 2) the result of the previous iteration, in order to "measure" change
	def evalIteration(self):
		findThesaurus = lambda rec: join(rec["output"], 
			basename(rec["input"])+".sims.neighbours.strings")
		evalOutput = "eval-iter"+str(len(self.record))+".tmp"
		
		## compute similarity with WordNet
		with open(findThesaurus(self.record[-1]), 'r') as currentTh:
			evalTask1 = bybloeval.BybloEval([currentTh], evalOutput, method="rank", 
				testIndex=0, maxRank=None, maxIndex=None, verbose=True)
			evalTask1.run()
			
			## write evaluation result to record
			self.record[-1]["eval-vs-WN"] = open(evalOutput, 'r').readlines()[0].split('\t')[0]
			
		## compute similarity with previously built thesaurus
		if len(self.record) < 2:
			self.printer.info("Only possible from the second iteration. Skipped.")
		else: 
			with open(findThesaurus(self.record[-2]), 'r') as previousTh,\
			open(findThesaurus(self.record[-1]), 'r') as currentTh:
				evalTask2 = bybloeval.BybloEval([currentTh, previousTh], evalOutput, method="rank", 
					testIndex=0, maxRank=None, maxIndex=None, verbose=True)
				evalTask2.run()
				
				## write evaluation result to record
				self.record[-1]["eval-vs-prev"] = open(evalOutput, 'r').readlines()[0].split('\t')[0]
				
		## delete temporary file
		os.remove(evalOutput)
		print self.record[-1]
	
	
	def plotIteration(self):
		## retrieve settings strings
		findSettings = lambda dict: dict["settings"]
		allSettingsStrings = map(findSettings, self.record)
		print allSettingsStrings
		
		#~ byblostats.generateHistograms(self.inputFile, allSettingsStrings, a.outputDir, a.bybloDir, a.reuse, a.verbose, a.cut)
		print "something"
	
	
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
	

	