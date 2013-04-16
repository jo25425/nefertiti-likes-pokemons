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

class Printer:
	"""
	Very basic object whose function is to simplifying output formatting within an application.
	Can print several types of elements: a main (module/class/task) title, some standard information,
	the current stage in a process, an array or dictionary, and the content of a file.
	Takes a verbose option used initially into account, which avoids having to pass it as an argument
	each type.
	"""
	
	## Sets verbose option
	def __init__(self, verbose):
		self.verbose = verbose
	
	## Important title, e.g. class or module
	def mainTitle(self, string):
		if self.verbose:
			for i in range(80/2): print '*',
			print '\n*' + (string.upper()).center(77) + '*'
			for i in range(80/2): print '*',
			print '\n'
	
	## Standard information
	def info(self, string):
		if self.verbose:
			print "[INFO] " + string
	
	## Current stage in a process
	def stage(self, number, total, title):
		if self.verbose:
			print  "\n+++ Stage " + str(number) + " of " + str(total) + ": " + title + " +++\n"

	## Array or dictionary with customisable start, end, line length and title
	def lines(self, list, min=0, max=None, line_max=None, title="List"):
		if self.verbose:
			stop = lambda x, L: x if x and x < len(L) else len(L)
			stopList = stop(max, list)
			
			print title
			
			## dictionary
			if type(list) is dict:
				for index, key in enumerate(sorted(list.iterkeys())):
					if index in range(min, stopList):
						if type(list[key]) is dict:
							print_lines(list[key], title=key) # recursive print of nested dictionaries
						else:
							lineText = key + '\t' + str(list[key])
							stopText = stop(line_max, lineText)
							print "   " + lineText[:stopText] + ("..." if stopText < len(lineText) else "")
			## list
			else:
				for index in range(min, stopList):
					lineText = str(list[index])
					stopText = stop(line_max, lineText)
					print lineText[:stopText] + ("..." if stopText < len(lineText) else "")
				
			print "..." if stopList != len(list) else ""
	
	## File with a customisable start and title
	def file(self, fileName, min=0, title="File content"):
		print title
		fileLines = [ line[:line.index('\n')] for line in open(fileName, 'r') ]
		self.lines(fileLines, min=min, title="From " + fileName)
