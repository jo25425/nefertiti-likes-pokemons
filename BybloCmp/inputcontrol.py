"""
Blah blah blah.
...
...
Blah.
"""

import argparse

def checkBybloSettings(input, type="studied and fixed"):
	description = "Choose " + type + " settings for Byblo."
	
	bybloParser = argparse.ArgumentParser(description=description)
	
	if "studied" in type:
		## minimum entry frequency threshold
		bybloParser.add_argument('-fef', '--filter-entry-freq',  dest='fef', 
			action='store', default=0.0,
			help='Minimum entry frequency threshold (default: 0.0)')
		## minimum feature frequency threshold
		bybloParser.add_argument('-fff', '--filter-feature-freq',  dest='fff', 
			action='store', default=0.0,
			help='Minimum feature frequency threshold (default: 0.0)')
		## minimum event frequency threshold
		bybloParser.add_argument('-fvf', '--filter-event-freq',  dest='fvf', 
			action='store', default=0.0,
			help='Minimum event frequency threshold (default: 0.0)')
		## minimum similarity threshold
		bybloParser.add_argument('-Smn', '--similarity-min', dest='Smn',
			action='store', default=None,
			help='Minimum similarity threshold (default: none)')
			
			
	if "fixed" in type:
		## caracter set encoding to use for both reading input and writing output files
		bybloParser.add_argument('-c', '--charset',  dest='charset', 
			action='store', default="UTF-8",
			help='Caracter set encoding to use for both reading input and  writing '+
			'output files (default: "UTF-8")')
		## beta parameter to CRMI measure
		bybloParser.add_argument('--crmi-beta',  dest='crmi_beta', 
			action='store', default=0.5,
			help='Beta parameter to CRMI measure (default: 0.5)')
		## gamma parameter to CRMI measure
		bybloParser.add_argument('--crmi-gamma',  dest='crmi_gamma', 
			action='store', default=0.5,
			help='Gamma parameter to CRMI measure (default: 0.5)')
		## regular expression that accepted entries must match
		bybloParser.add_argument('fep', '--filrer-entry-pattern', dest='fep',
			action='store',
			help='Regular expression that accepted entries must match')
		## whitelist file containing entries of interest
		bybloParser.add_argument('-few', '--filter-entry-whitelist', dest='few',
			action='store',
			help='Whitelist file containing entries of interest (all others will be ignored)')
		## regular expression that accepted features must match
		bybloParser.add_argument('ffp', '--filrer-feature-pattern', dest='ffp',
			action='store',
			help='Regular expression that accepted features must match')
		## whitelist file containing features of interest
		bybloParser.add_argument('-ffw', '--filter-feature-whitelist', dest='ffw',
			action='store',
			help='Whitelist file containing features of interest (all others will be ignored)')
		## produce similarity between pair of identical entries
		bybloParser.add_argument('-ip', '--identity-pairs', dest='ip',
			action='store_true', default=False,
			help='Produce similarity between pair of identical entries (default: False)')
		## alpha parameter to Lee alpha-skew divergence measure
		bybloParser.add_argument('--lee-alpha', dest='lee_alpha',
			action='store', default=0.99,
			help='Alpha parameter to Lee alpha-skew measure (default: 0.99)')
		## similarity measure to use
		bybloParser.add_argument('-m', '--measure', dest='measure',
			action='store', default='Lin',
			help='SImilarity measure to use (default: "Lin")')
		## swap similarity measure inputs
		bybloParser.add_argument('--measure-reversed', dest='m_reversed',
			action='store_true', default=False,
			help='Swap similarity measure inputs (default: False')
		## p parameter to MInkowski distance measure
		bybloParser.add_argument('-mink-p', dest='mink_p',
			action='store', default=2.0,
			help='P parameter to Minkowski distance measure (default: 2.0)')
		## maximum similarity threshold
		bybloParser.add_argument('-Smx', '--similarity-max', dest='Smx',
			action='store', default=None,
			help='Maximum similarity threshold (default: none)')
		## temporary directory
		bybloParser.add_argument('-T', '--temp-dir', dest='temp_dir',
			action='store', default='subdirectory',
			help='Temporary directory, used during processing '+
			'(default: subdirectory created inside the outptu directory)')
		## number of concurrent processing threads
		bybloParser.add_argument('-t', '--threads', dest='threads',
			action='store', default=4,
			help='NUmber of concurreny processing threads (default: 4)')
		## maximum number of neighbours to produce per word
		bybloParser.add_argument('-k', dest='k',
			action='store', default=100,
			help='Maximum number of neighbours to produce per word (default:100)')
	
	try:
		return bybloParser.parse_args(input.split())
	except:
		return False
	