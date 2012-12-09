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
import multiprocessing
import itertools
import os
from os.path import *
from math import ceil
from nltk.corpus import wordnet, brown
from operator import itemgetter

	
## Replaces an 's' (for ADJ SAT) PoS tag by a simple 'a' (for ADJ) in order to add some smoothing
## to the processing of PoS tags in the IC (Information Content) of a corpus
## @return same synset, possibly with a smoothed PoS tag
def smoothed(synset):
	if synset.pos == 's':
		synset.pos = 'a'
	return synset
		
		
## Lin's WordNet similarity function
## @return similarity score, or None if there is no valid score
def linWordnetSimilarity(word1, word2, ic):
	
	possibleValues = [wordnet.lin_similarity(synset1, synset2, ic) \
				for synset1 in [smoothed(s) for s in wordnet.synsets(word1)] \
				for synset2 in [smoothed(s) for s in wordnet.synsets(word2)] \
				if synset1.pos == synset2.pos and synset1.pos in ic]
	possibleValues = [v for v in possibleValues if v >= 0]
	return max(possibleValues) if possibleValues else None


class ThesAuto:
	
	"""
	Build a thesaurus for a set of words based on similarity scores determined using WordNet.
	
	The set of words can be specified via a file - the first field of each line being a word 
	- or will be the entire set of words appearing in WordNet.
	This tool comes with a database containing this entire set.
	It is possible to choose the number of neighbours kept for each word.
	"""
	
	## Initialises the parameters of the module
	def __init__(self, inputFile, outputFile, maxRank, discard, verbose, nbCPU):
		## input file for words
		self.inputFile = inputFile 
		## output file for thesaurus
		self.outputFile = outputFile
		## fixed maximum rank distance
		self.maxRank = maxRank
		## discard option for words
		self.discard = discard
		## verbose option
		self.verbose = verbose
		## number of processors
		self.nbCPU = nbCPU if nbCPU else multiprocessing.cpu_count()
		## WordNet words database
		self.database = join(dirname(realpath(__file__)), ".wordnet.words")
		## WordNet thesaurus
		self.thesaurus = []
	

	
	## Extracts from the specified files the term that are to appear in the thesaurus.
	## Each term is the first field of a line.
	## @return list of extracted terms
	def extractTerms(self, file, database="./.wordnet.words", discard=True, verbose=False):
		terms = []
		tryDb, makeDb = file is None, False
		
		## read terms from file
		if file:
			try:
				for line in open(file, 'r'):
					## get 1st field
					t = string.split(line)[0]
					## keep if in WordNet
					if wordnet.synsets(t):
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
			for synset in list(wordnet.all_synsets()):
				terms += [lemma.name for lemma in synset.lemmas \
					if lemma.name not in terms]
			## write them in a file
			if not discard:
				open(database, 'w').write('\n'.join(terms))
				os.chmod(database, stat.SF_IMMUTABLE)
		return terms		
	
	
	##
	##
	def createWordnetFullThesaurus(self, words):
		## prepare WordNet IC and other tools
		brownIC = wordnet.ic(brown, False, 1.0)
		processes, fullThesaurus = [], []
		chunkSize = int( ceil( 1.0*len(words)/self.nbCPU ) )
		resultQueue = multiprocessing.Queue()
		
		## compute thesaurus chunks in different processes
		for counter in range(self.nbCPU):
			chunkWords = words[counter*chunkSize:min( (counter+1)*chunkSize, len(words) )]
			p = multiprocessing.Process(target=self.createWordnetPartialThesaurus, 
				args=(resultQueue, chunkWords, words, brownIC, counter))
			processes.append(p)
			p.start()
		
		return sorted(itertools.chain.from_iterable(resultQueue.get() for p in processes))
		
	
	
	##
	##
	def createWordnetPartialThesaurus(self, queue, selectedWords, allWords, ic, nb):
		print "-> starting  for part", nb, "with", len(selectedWords), "words"
		partialThesaurus =  [set for set in [self.createWordnetNeighbourSet( \
					allWords, targetWord, ic, self.maxRank) \
					for targetWord in selectedWords] if set]
		queue.put(partialThesaurus)
		
		print "<- ending process", nb


	## Creates a neighbour set for a target word using a list of words available and an ic 
	## (Information Content of a corpus). This set can be limited to a chosen number of 
	## neighbours k.
	## Synsets whose PoS doesn't appear in the ic corpus are ignored.
	## @return neighbour set (type: [string, (string, float), (string, float), ...]
	def createWordnetNeighbourSet(self, words, targetWord, ic, k=None, verbose=False):
		## neighbour set creation
		if verbose:
			print "Creating set for " + targetWord

		set = [(w, linWordnetSimilarity(targetWord, w, ic)) 
			for w in words if w != targetWord]
		set.sort(key=itemgetter(1), reverse=True)
		
		return ([targetWord] + set[:k]) if set else None


	## Converts a neighbour set (containing tuples) to a printable string
	## @return string for the neighbour set
	def setToString(self, set):
		neighbourToString = lambda tuple: tuple[0] + '\t' + str(tuple[1])
		return '\t'.join([set[0]] + [neighbourToString(w) for w in set[1:]]) + '\n'

	
	## Runs the automatic construction of a thesaurus
	def run(self):
		## start operations
		stime = time.time()
		print "***************************************************************************"
		print "WORDNET THESAURUS CREATION TOOL"
		print "***************************************************************************\n"
		
		
		
		## read file
		print "[1] Extracting words from input file\n" + \
		      "------------------------------------"
		terms = self.extractTerms(self.inputFile, self.database, self.discard, self.verbose)
		#~ if self.verbose:
			#~ self.print_lines(terms, max=20, title="-- Extracted " + str(len(terms)) + " terms --")
		self.print_lines(terms, max=20, 
			title="-- Extracted " + str(len(terms)) + " terms --")
		
		
		
		## create thesaurus
		print "\n[2] Creating thesaurus in output file\n" + \
		      "-------------------------------------"
		      
		
		thesaurus = self.createWordnetFullThesaurus(terms)
		self.print_lines(thesaurus, max=20, line_max=75, \
			title="\n\n-- Constructed " + str(len(thesaurus)) + " neighbour sets --")
			
		with open(self.outputFile, 'w') as output:
			for set in thesaurus:
				output.write( self.setToString(set) )
		
		
		
		etime = time.time()
		print "\n>Execution took", etime-stime, "seconds"

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
	parser.add_argument('-k', '--max-rank', metavar='n', type=int, dest='maxRank',
		action='store', default=100,
		help='number of neighbours to keep (default: 100)')
	# discard option for words
	parser.add_argument('-d', '--discard', dest='discard',
		action='store_true', default=False,
		help='discard the record of words that appear in WordNet (default: False)')
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', 
		action='store_true', default=False,
		help='display information about operations (default: False)')
	# number of processors
	parser.add_argument('-p', '--processors', metavar='n', dest='processors',
		action='store',
		help='number of processors to use (default: available)')
		
	a = parser.parse_args()
	thesAuto = ThesAuto(a.inputFile, a.outputFile, a.maxRank, a.discard, a.verbose, a.processors)
	thesAuto.run()
	

	