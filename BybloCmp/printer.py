"""
Loads of useless stuff here...
'cause that's what academics like.
"""

def info(string, verbose=True):
	if verbose:
		print [INFO]
		print string
	

def stage(number, total, title):
	print  "\n+++ Stage " + str(number) + " of " + str(total) + ": " + title + " +++\n"
	

def mainTitle(string, verbose=True):
	if verbose:
		for i in range(80/2): print '*',
		print '*' + (string.upper()).center(77) + '*',
		for i in range(80/2): print '*',
		print ''