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
import string
import sys
import os
from os.path import *

## global variable holding extensions for output files of Byblo
BYBLO_OUTPUT_EXTENSIONS= [	".entries", ".entries.filtered", 
						".entry-index.d.0", ".entry-index.i.0",
						".enumerated",
						".events", ".events.filtered",
						".feature-index.d.0", ".feature-index.i.0",
						".features", ".features.filtered",
						".sims", ".sims.neighbours", ".sims.neighbours.strings",
						".runtime"
					]

## Produces a string that gives information about parameters used and can be in a file name
## @return formatted parameters string
def paramSubstring(str):
	return '#' + string.replace(str, ' ', '_')


def extractText(fileName, lineIndex, start, end, mod):
	text = list(open(fileName, 'r'))[lineIndex]
	if text.startswith(start) and text.endswith(end):
		return mod(text[ len(start) : text.index(end) ])
	else:
		return -1

def retrieveInputFile(logFileName):
	start, end = " * Input instances file: ", "\n"
	identity = lambda x: x
	inputFile = extractText(logFileName, 7, start, end, identity)
	return inputFile
	
def retrieveRunningTime(logFileName):
	start, end = " * Elapsed time: ", "\n"
	toSeconds = lambda timeString: sum( [float(t)*60**(2-i) for i, t in enumerate( timeString.split(':') )] )
	time = extractText(logFileName, -4, start, end, toSeconds)
	return time 

## Renames ALL of the output files of Byblo so that their name includes the parameter string used
## by Byblo to create them.
def renameBybloOutputFiles(parameters, base):
	for ext in BYBLO_OUTPUT_EXTENSIONS:
		filename = base + ext
		if exists(filename):
			newname = base + paramSubstring(parameters) + ext
			os.rename(filename, newname)

## Moves ALL of the output files of Byblo to put them in a new directory
def moveBybloOutputFiles(parameters, base, newdir):
	newbase = join(newdir, basename(base))
	for ext in BYBLO_OUTPUT_EXTENSIONS:
		filename = base + paramSubstring(parameters) + ext
		if exists(filename):
			newname = newbase + paramSubstring(parameters) + ext
			os.rename(filename, newname)

if __name__=='__main__':
	argParser = argparse.ArgumentParser(description='Prepare Byblo output for use by BybloCmp.')
	
	## log file for Byblo log
	argParser.add_argument(metavar='file', dest='logFile', action='store', 
		help='log file for Byblo execution')
	## thesauri dir for Byblo output
	argParser.add_argument(metavar='dir', dest='thesauriDir', action='store', 
		help='thesauri directory for Byblo output')
	## parameters for Byblo log
	argParser.add_argument(metavar='string', dest='parameters', action='store', 
		help='parameters that were used by Byblo to produce this log')
		
	a = argParser.parse_args()
	## ilog file for Byblo log
	logFile = abspath(a.logFile)
	## thesauri dir for Byblo output
	thesauriDir = abspath(a.thesauriDir)
	## parameters for Byblo log
	parameters = a.parameters
	
	print "** Preparing Byblo output for use by BybloCmp **"
	
	## 1/3)  Create a file containing the running time
	## extract input file name
	inputFileName = retrieveInputFile(logFile)
	print "Input file:", inputFileName
	## extract running time
	runningTime = retrieveRunningTime(logFile)
	print "Running time:", runningTime
	
	## create the file (correctly named)
	if inputFileName == -1 or runningTime == -1:
		print "Error reading log file."
		sys.exit()
	base = join(thesauriDir, basename(inputFileName))
	open(base + paramSubstring(parameters) + ".runtime", 'w').write(str(runningTime))
	print "Runtime file written!"

	## 2/3)  Rename the output files created by Byblo
	renameBybloOutputFiles(parameters, base)
	print "Byblo output files renamed!"
	
	## 3/3)  Move all of these files (including the runtime one) to a "thesauri" subdirectory
	## create thesauri directory when required
	newThesauriDir = join(thesauriDir, "thesauri")
	if not exists(newThesauriDir): os.makedirs(newThesauriDir)
	moveBybloOutputFiles(parameters, base, newThesauriDir)
	print "All files moved!"
	print "Everything done! :)"
	