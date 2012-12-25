"""
Loads of useless stuff here...
'cause that's what academics like.
"""
	

def mainTitle(string, verbose=True):
	if verbose:
		for i in range(80/2): print '*',
		print '*' + (string.upper()).center(77) + '*',
		for i in range(80/2): print '*',
		print ''
	
	
def info(string, verbose=True):
	if verbose:
		print "[INFO] " + string
	

def stage(number, total, title):
	print  "\n+++ Stage " + str(number) + " of " + str(total) + ": " + title + " +++\n"


## Prints an array with customizable start, end, line length and title
def lines(list, min=0, max=None, line_max=None, title="List", verbose=True):
	if verbose:
		stop = lambda x, L: x if x and x < len(L) else len(L)
		
		print title
		
		stopList = stop(max, list)
		for index in range(min, stopList):
			lineText = str(list[index])
			stopText = stop(line_max, lineText)
			print lineText[:stopText] + ("..." if stopText < len(lineText) else "")
			
		print "..." if stopList != len(list) else ""