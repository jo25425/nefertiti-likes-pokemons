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
__version__ = "1.1.0"
__maintainer__ = "Joanne Robert"
__email__ = "jr317@sussex.ac.uk"
__status__ = "Development"
#~ -----------------------------------------------------------------------------

import sys
import argparse
import time
import string
import stat
import multiprocessing
import math
import os
from os.path import *
import thesauto
import outputformatting


## global variable holding names of special entries in thesauri
SPECIAL_ENTRIES = ["___FILTERED___"]

## Lin's WordNet similarity function
def linWordnetSimilarity(word1, word2):
		return max([wn.lin_similarity(synset1, synset2, ic) \
				for synset1 in wn.synsets(word1) \
				for synset2 in wn.synsets(word2) \
				if synset1.pos == synset2.pos \
				and synset1.pos in ic and synset2.pos in ic])


## custom type for positive integer (maximum neighbour rank)
def positiveInt(text):
	value = int(text)
	if not 0 < value:
		raise argparse.ArgumentTypeError("stricly positive integer expected")
	return value


class BybloEval:
	
	"""
	Compute the percentage of similarity between two thesauri
	"""
	
	##Initialises the parameters of the module
	def __init__(self, inputFiles, baseThesaurus, outputFile, method, testIndex, 
		maxRank, maxIndex, verbose):
		## input files containing the thesauri
		self.inputFiles =  inputFiles
		## base WordNet thesaurus for evaluation against WordNet 
		inputName = self.inputFiles[0].name
		self.baseThesaurus = inputName + ".WN" if not baseThesaurus\
			else (abspath(baseThesaurus) if not isdir(baseThesaurus)
			else join(abspath(baseThesaurus), basename(inputName) + ".WN"))
		## output file for comparison results
		self.outputFile = abspath(outputFile) if not isdir(outputFile)\
			else join(abspath(outputFile), "result.thsim")
		## measure used for neighbour set (i.e. thesaurus entry) comparison
		self.method = method
		## index for a single test entry
		self.testIndex = testIndex
		## maximum neighbour rank
		self.maxRank = maxRank
		## maximum entry index
		self.maxIndex = maxIndex
		## verbose option
		self.verbose = verbose
		## self.printer to output
		self.printer = outputformatting.Printer(self.verbose)
	
	
	## Runs the comparison between two thesauri
	def run(self):
		
		## start operations
		stime = time.time()
		self.printer.mainTitle("Thesaurus comparison tool")
		
		## [1/3] read files
		self.printer.stage(1, 3, "Extracting words from input files")
		thesauri = [self.extractTerms(f) for f in self.inputFiles]
		names = [f.name for f in self.inputFiles]
		
		if len(self.inputFiles) == 1:
			thesaurusFile = names[0] + ".WN" # corresponds to default name in thesauto
			self.createWordnetThesaurus(names[0], thesaurusFile, self.baseThesaurus, self.maxRank)
			thesauri.append(self.extractTerms( open(thesaurusFile, 'r') ))
			names.append("WordNet")
		
		## verify extracted information
		for th, name in zip(thesauri, names):
			self.printer.lines(th, max=10, line_max=75, title="\""+name+"\" after sort:")
		
		## [2/3] compare one entry as a test
		self.printer.stage(2, 3, "Comparing neighbour sets for test index " + str(self.testIndex))
		set1, set2 = [th[self.testIndex] for th in thesauri]
		sim_score1 = {
			  'Lin': self. neighbourSetSimLin,
			  'rank': self.neighbourSetSimRank
		}[self.method](set1, set2, verbose=self.verbose)
		self.printer.info("Similarity score: " + str(sim_score1) + "\n")
		
		## [3/3]compare thesauri
		self.printer.stage(3, 3, "Comparing thesauri") 
		with open(self.outputFile, 'w') as output:
			sim_score2 = self.thesaurusSim(thesauri, self.method, 
				n=self.maxIndex, verbose=self.verbose)
			output.write(str(sim_score2) + '\t' + '\t'.join(names))
		self.printer.info("Similarity score: " + str(sim_score2) + "\n")

		etime = time.time()
		self.printer.info("Execution took " + str(etime-stime) + " seconds\n")
	
	
	## Extracts thesaurus entries from a file
	## The entries are formatted as follows:
	##	term 	neighbour1	sim_score1	neighbour2	sim_score2	...	...
	##
	## When the similarity measure chosen for comparison between neighbour sets uses ranks,
	## similarity scores are replaced by their rank in the resultant array.
	## @return array representing a thesaurus
	def extractTerms(self, file):
		terms = []
		for line in file:
			fields = string.split(line, "\t")
			
			if fields[0] not in SPECIAL_ENTRIES:
				if self.maxRank: fields = fields[:self.maxRank*2+1]
				
				## build tuples depending on method chosen
				if self.method == "Lin":
					line_terms = [(fields[i], float(fields[i+1])) for i in xrange(1, len(fields)) \
						if i%2 == 1 and fields[i] not in SPECIAL_ENTRIES]
				elif self.method == "rank":
					line_terms = [(fields[i], i/2) for i in xrange(1, len(fields)) \
						if i%2 == 1 and fields[i] not in SPECIAL_ENTRIES]
					
				line_terms.sort()
				line_terms = [fields[0]] + line_terms
				terms.append(line_terms)
		terms.sort()
		return terms


	## Creates a thesaurus from WordNet, using the terms already present in an existing thesaurus
	## @return array representing a neighbour set
	def createWordnetThesaurus(self, inputFile, thesaurusFile, baseThesaurus, maxRank=None):
		thesautoTask = thesauto.ThesAuto(inputFile, thesaurusFile, baseThesaurus,
			maxRank, discard=True, verbose=self.verbose)
		thesautoTask.run()
		
		self.printer.info("Thesaurus written in " + thesaurusFile)
	
	
	## Computes the similarity score between two versions of a neighbour set (coming from different
	## thesauri) using Lin's set similarity measure.
	## @return float in range [0, 1]
	def  neighbourSetSimLin(self, set1, set2, verbose=False):
		sim_score, pt1, pt2 = 0, 1, 1
		square_score = lambda x: x[1] * x[1]
		
		## determine maximum possible similarity score (set specific)
		max_sim_score = sum([square_score(neighbour) for neighbour in set1[1:]]) \
			* sum([square_score(neighbour) for neighbour in set2[1:]])
		max_sim_score = math.sqrt(max_sim_score)
		
		## avoid division by zero
		if max_sim_score == 0:
			return 0
			
		while pt1 < len(set1) and pt2 < len(set2) :
			if set1[pt1][0] < set2[pt2][0]:
				pt1 += 1
				
			elif set1[pt1][0] > set2[pt2][0]:
				pt2 += 1
			
			else: # potential neighbour present in both sets
				score_product = set1[pt1][1]*set2[pt2][1]
				sim_score += score_product
				
				if verbose: # what is happening?
					print str(set1[pt1][1]) + ' vs. ' + str(set2[pt2][1]) +' \t'\
						+ set1[pt1][0] + ' \t=> ' + str(score_product)
				
				pt1 += 1
				pt2 += 1
				
		return 1.0 * sim_score / max_sim_score
	
	
	## Computes the similarity score between two versions of a neighbour set (coming from different
	## thesauri) using an adaptation of Lin's set similarity measure based on ranks and not values. This
	## allows to compare thesauri featuring scores obtained with different similarity measures.
	## @return float in range [0, 1]
	def neighbourSetSimRank(self, set1, set2, verbose=False):
		sim_score, pt1, pt2 = 0, 1, 1
		
		## determine maximum rank distance to compute maximum similarity score
		k = max(len(set1), len(set2)) - 1 # target word doesn't count
			
		## compute maximum similarity score 
		max_sim_score = sum([i*i for i in xrange(1, k+1)])
			
		## avoid division by zero
		if max_sim_score == 0:
			return 0
		
		while pt1 < len(set1) and pt2 < len(set2) :
			if set1[pt1][0] < set2[pt2][0]:
				pt1 += 1
					
			elif set1[pt1][0] > set2[pt2][0]:
				pt2 += 1
				
			else: # potential neighbour present in both sets
				score_product = (k - set1[pt1][1])*(k - set2[pt2][1])
				sim_score += score_product
				
				if verbose: # what is happening?
					print str(set1[pt1][1]) + ' vs. ' + str(set2[pt2][1]) +' \t'\
						+ set1[pt1][0] + '\t => ' + str(score_product)
				
				pt1 += 1
				pt2 += 1
				
		return 1.0 * sim_score / max_sim_score
	
	
	## Computes the similarity between two thesauri, using one of the two measures available (based on
	## either values or ranks), and with the possibility to limit the comparison to a given number of entries
	## (i.e. neighbour sets).
	## @return float in range [0, 1]
	def thesaurusSim(self, thesauri, method="rank", n=None, verbose=False):
		th1, th2 = thesauri[0], thesauri[1]
		global_sim_score, pt1, pt2 = 0, 0, 0
		neighbour_set_sim = {
			'Lin': self.neighbourSetSimLin,
			'rank': self.neighbourSetSimRank
		}[method]
		
		## determine maximum thesaurus size if not passed (n allows comparison on a subset)
		if n == None:
			n = max(len(th1), len(th2))
			
		## avoid division by zero
		if n == 0:
			return 0
		
		while pt1 < n and pt2 < n and pt1 < len(th1) and pt2 < len(th2) :
			
			if th1[pt1][0] == th2[pt2][0]: # target word present in both sets
				set_sim = neighbour_set_sim(th1[pt1], th2[pt2])
				global_sim_score += set_sim
				
				if verbose: # what is happening?
					print th1[pt1][0] + ' \t=> ' + str(set_sim)
				
				pt1 += 1
				pt2 += 1
				
			elif th1[pt1][0] < th2[pt2][0]:
				pt1 += 1
			else:
				pt2 += 1
		return 1.0 * global_sim_score / n



## Parses a command line
if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Compare two thesauri.')
	
	## Forces one or two values exactly for an argument
	def required_length(nmin,nmax):
		class RequiredLength(argparse.Action):
			def __call__(self, parser, args, values, option_string=None):
				if not nmin<=len(values)<=nmax:
					msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format( \
						f=self.dest,nmin=nmin,nmax=nmax)
					raise argparse.ArgumentTypeError(msg)
				setattr(args, self.dest, values)
				return RequiredLength

	## input files containing the thesauri
	parser.add_argument('inputFiles', metavar='file', type=file, nargs='+',
		action=required_length(1, 2),
		help='files containing the thesauri to compare (against WordNet if there is only one)')
	## base WordNet thesaurus for evaluation against WordNet 
	parser.add_argument('-th', '--thesaurus', metavar='file', dest='baseThesaurus',
		action='store',
		help='base WordNet thesaurus for evaluation against WordNet (specific to input) ' +
			'(default: input file name + \".WN\")')
	## output file for comparison results
	parser.add_argument('-o', '--output-file', metavar='file', dest='outputFile',
		action='store', default="./result.thsim",
		help='output file in which similarity scores and file names will be written'+\
		'(default: "./result.thsim")')
	## measure used for neighbour set (i.e. thesaurus entry) comparison
	parser.add_argument('-l', '--Lin', metavar='string', dest='method', action='store_const',
		const='Lin', default='rank',
		help='method used to calculate the similarity between two neighbour sets'+\
		'(default: rank)')
	## index for a single test entry
	parser.add_argument('-i', '--test-index', metavar='n', type=int, dest='testIndex',
		action='store', default=0,
		help="index of a single word to compare as a test (default: 0)")
	## maximum neighbour rank
	parser.add_argument('-k', '--max-rank', metavar='n', type=positiveInt, dest='maxRank', 
		action='store', default=None,
		help="maximum rank for neighbours to compare")
	## maximum entry index
	parser.add_argument('-n', '--max-index', metavar='n', type=positiveInt, dest='maxIndex',
		action='store', default=None,
		help="maximum index for entries to compare")
	## verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_const',
		const=True, default=False,
		help="display information about each comparison (default: False)")
		
	a = parser.parse_args()
	bybloEval = BybloEval(a.inputFiles, a.baseThesaurus, a.outputFile, a.method, a.testIndex, 
		a.maxRank, a.maxIndex, a.verbose)
	bybloEval.run()