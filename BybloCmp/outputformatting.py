"""
Loads of useless stuff here...
'cause that's what academics like.
"""
	
class Printer:
	
	def __init__(self, verbose):
		self.verbose = verbose
		 
	def mainTitle(self, string):
		if self.verbose:
			for i in range(80/2): print '*',
			print '\n*' + (string.upper()).center(77) + '*'
			for i in range(80/2): print '*',
			print '\n'
	
	
	def info(self, string):
		if self.verbose:
			print "[INFO] " + string
	

	def stage(self, number, total, title):
		if self.verbose:
			print  "\n+++ Stage " + str(number) + " of " + str(total) + ": " + title + " +++\n"


	## Prints an array with customizable start, end, line length and title
	def lines(self, list, min=0, max=None, line_max=None, title="List"):
		if self.verbose:
			stop = lambda x, L: x if x and x < len(L) else len(L)
			
			print title
			
			stopList = stop(max, list)
			for index in range(min, stopList):
				lineText = str(list[index])
				stopText = stop(line_max, lineText)
				print lineText[:stopText] + ("..." if stopText < len(lineText) else "")
				
			print "..." if stopList != len(list) else ""
	
	def record_iter(self, out="standard"):
		if out is "standard":
			print "to standard output"
		else:
			print "to file"
		

""""
Let's make a class out of this
"""

