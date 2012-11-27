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
from menusystem import Menu, Choice


class BybloCmp:
	
	"""
	Compare methods for the construction of a disitrbutional thesaurus
	"""
	
	## Initialises the parameters of the module
	def __init__(self, inputFile, bybloDir, storageDir, outputFile, verbose):
		
		## input file for pre-processed data
		self.inputFile = inputFile 
		## location of the Byblo directory
		self.bybloDir = bybloDir
		## storage directory for Byblo output
		self.storageDir = storageDir
		## output file for  comparison results 
		self.outputFile = outputFile
		## verbose option
		self.verbose = verbose
		## options menu
		self.optionsMenu = self.initOptionsMenu()
		## main menu
		self.mainMenu = self.initMainMenu()
	
	
	## Creates the main menu.
	## It allows to navigate in the program, accessing all available functionalities."
	def initMainMenu(self):
		mainChoices = [
			Choice(1, description="Execution of Byblo iterations"),
			Choice(2, description="Iterations planning"),
			Choice(3, description="Comparison results"),
			Choice(4, description="Options", subMenu=self.optionsMenu),
			Choice(5, description="Exit", value='Bye!', handler=lambda anything:False)
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
	
	
	def	execution(self):
		
		## get settings
		
		## run
		
		## compare
		
		## display results of comparison
	
	
	## Run
	def run(self):
		## Start operations.
		stime = time.time()
		print "***************************************************************************"
		print "PARAMETERISATION COMPARISON TOOL"
		print "***************************************************************************\n"
		
		## direct operations from main menu
		self.mainMenu.waitForInput()
		
		"""
		This will eventually be modified so that iterations can be launched using scripts.
		1) a script used to run Byblo each time, but menu still used
		2) a general script that includes a configuration file holding all of the information for iteration planning
		At least ideally....
		"""
		
		etime = time.time()
		print "\n>Execution took", etime-stime, "seconds"   
	

## Parses a command-line.
if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Compare different parameterisations for Byblo.')
	
	## input file for preprocessed data
	parser.add_argument(metavar='file', dest='inputFile', action='store', 
		help='input file for preprocessed data')
	## location of the Byblo directory
	parser.add_argument('-b', '--byblo', metavar='dir', dest='bybloDir',
		action='store', default="./Byblo-2.1.0",
		help='directory for the location Byblo (default: "./Byblo-2.1.0")')
	## storage directory for Byblo output
	parser.add_argument('-s', '--storage', metavar='dir', dest='storageDir',
		action='store', default="./byblo-output",
		help='storage directory for Byblo output (default: "./byblo-output")')
	## output file for comparison results
	parser.add_argument('-o', '--output', metavar='file', dest='outputFile',
		action='store', default="./results.cmp",
		help='output file for comparison results (default: "./results.cmp")')
	## verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', 
		action='store_true', default=False,
		help='display information about operations (default: False)')
		
	a = parser.parse_args()
	bybloCmp = BybloCmp(a.inputFile, a.bybloDir, a.storageDir, a.outputFile, a.verbose)
	bybloCmp.run()
	

	