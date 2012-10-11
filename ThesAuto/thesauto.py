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

import sys
import argparse
import time
import string
import stat
import os
from os.path import *
from nltk.corpus import wordnet as wn, wordnet_ic as wn_ic
from operator import itemgetter


class ThesAuto:
	
	"""
	Build a thesaurus for a set of words based on similarity scores determined using WordNet.
	
	The set of words can be specified via a file - the first field of each line being a word 
	- or will be the entire set of words appearing in WordNet.
	This tool comes with a database containing this entire set.
	It is possible to choose the number of neighbours kept for each word.
	"""
	
	## Initialises the parameters of the module
	def __init__(self, inputFile, outputFile, k, discard, verbose):
		## input file for words
		self.inputFile = inputFile 
		## output file for thesaurus
		self.outputFile = outputFile
		## fixed maximum rank distance
		self.k = k
		## discard option for words
		self.discard = discard
		## verbose option
		self.verbose = verbose
		## WordNet words database
		self.database = join(dirname(realpath(__file__)), ".wordnet.words")
	
	
	## Runs the automatic construction of a thesaurus
	def run(self):
		## start operations
		stime = time.time()
		print "***************************************************************************"
		print "WORDNET THESAURUS CREATION TOOL"
		print "***************************************************************************\n"
		
		## prepare WordNet IC
		#~ brown_ic = wn_ic.ic('ic-brown.dat')
		semcor_ic = wn_ic.ic('ic-semcor.dat')
		
		## read file
		print "[1] Extracting words from input file\n" + \
		      "------------------------------------"
		terms = self.extract_terms(self.inputFile, self.database, self.discard, self.verbose)
		if self.verbose:
			self.print_lines(terms, max=20, title="-- Extracted " + str(len(terms)) + " terms --")
		
		## create thesaurus
		print "\n[2] Creating thesaurus in output file\n" + \
		      "-------------------------------------"
		with open(self.outputFile, 'w') as output:
			for t in terms:
				newSet = self.create_wordnet_neighbour_set(terms, t, semcor_ic, self.k, self.verbose)
				output.write( self.set_to_string(newSet) )
		
		etime = time.time()
		print "\n>Execution took", etime-stime, "seconds"   
	
	
	## Extracts from the specified files the term that are to appear in the thesaurus.
	## Each term is the first field of a line.
	## @return list of extracted terms
	def extract_terms(self, file, database="./.wordnet.words", discard=True, verbose=False):
		terms = []
		tryDb, makeDb = file is None, False
		
		## read terms from file
		if file:
			try:
				for line in open(file, 'r'):
					## get 1st field
					t = string.split(line)[0]
					## keep if in WordNet
					if wn.synsets(t):
						terms.append(t)
					elif verbose:
						print t, "not found in WordNet"
			except IOError as e:
				tryDb = True
				
		## read terms from database
		if tryDb:
			try:
				terms = [string.split(line)[0] for line in open(database, 'r')]
			except IOError as e:
				makeDb = True
		
		## find terms in WordNet
		if makeDb:
			for synset in list(wn.all_synsets()):
				terms += [lemma.name for lemma in synset.lemmas \
					if lemma.name not in terms]
			## write them in a file
			if not discard:
				open(database, 'w').write('\n'.join(terms))
				os.chmod(database, stat.SF_IMMUTABLE)
		return terms


	## Lin's WordNet similarity function
	def lin_wordnet_sim(self, word1, word2, ic):
		possibleValues = [wn.lin_similarity(synset1, synset2, ic) \
					for synset1 in wn.synsets(word1) \
					for synset2 in wn.synsets(word2) \
					if synset1.pos == synset2.pos \
					and synset1.pos in ic and synset2.pos in ic]
		return max(possibleValues) if possibleValues else -1


	## Creates a neighbour set for a target word using a list of words available and an ic (Information Content
	## of a corpus). This set can be limited to a chosen number of neighbours k.
	## Synsets whose PoS doesn't appear in the ic corpus are ignored.
	## @return neighbour set [type: list of string and (string, float) tuples]
	def create_wordnet_neighbour_set(self, words, targetWord, ic, k=None, verbose=False):
		## neighbour set creation
		if verbose:
			print "Creating set for " + targetWord

		set = [(w, self.lin_wordnet_sim(targetWord, w, ic)) for w in words if w != targetWord]
		set.sort(key=itemgetter(1), reverse=True)
		if k:
			set = set[:k]
		
		return [targetWord] + set


	## Converts a neighbour set (containing tuples) to a printable string
	## @return string for the neighbour set
	def set_to_string(self, set):
		neighbour_to_string = lambda tuple: tuple[0] + '\t' + str(tuple[1])
		return '\t'.join([set[0]] + [neighbour_to_string(w) for w in set[1:]]) + '\n'


	## Prints an array with customizable start, end, line length and title
	def print_lines(self, list, min=0, max=None, line_max=None, title="List"):
		stop = lambda x, L: x if x and x < len(L) else len(L)
		
		print title
		
		stopList = stop(max, list)
		for index in range(min, stopList):
			lineText = str(list[index])
			stopText = stop(line_max, lineText)
			print lineText[:stopText] + ("..." if stopText < len(lineText) else "")
			
		print "..." if stopList != len(list) else ""



if __name__=='__main__':
	## parse command line
	parser = argparse.ArgumentParser(description='Compare two thesauri.')
	
	# input file for words
	parser.add_argument('-i', '--input', metavar='file', dest='inputFile',
		action='store', default="...ThesAuto/wordnet.words",
		help='input file for terms to look for in WordNet (default: "...ThesAuto/wordnet.words")')
	# output file for thesaurus
	parser.add_argument('-o', '--output', metavar='file', dest='outputFile',
		action='store', default="./wordnet.thesaurus",
		help='output file where the thesaurus will be written (default: "./wordnet.thesaurus")')
	# fixed maximum rank distance
	parser.add_argument('-k', '--max-rank', type=int, dest='k',
		action='store', default=100,
		help='number of neighbours to keep (default: 100)')
	# discard option for words
	parser.add_argument('-d', '--discard', dest='discard',
		action='store_true', default=False,
		help=' discard the record of words that appear in WordNet (default: False)')
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', 
		action='store_true', default=False,
		help='display information about operations (default: False)')
		
	thesAuto = ThesAuto(a.inputFile, a.outputFile, a.k, a.discard, a.verbose)
	thesAuto.run()
	

	