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
from nltk.corpus import brown, wordnet, wordnet_ic
from operator import itemgetter

import outputformatting


## global variable holding names of special entries in thesauri
SPECIAL_ENTRIES = ["___FILTERED___"]

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
	

## custom type for positive integer (maximum neighbour rank)
def positiveInt(text):
	value = int(text)
	if not 0 < value:
		raise argparse.ArgumentTypeError("stricly positive integer expected")
	return value
	
	
class ThesAuto:
	
	"""
	Build a thesaurus for a set of words based on similarity scores determined using WordNet.
	
	The set of words can be specified via a file - the first field of each line being a word 
	- or will be the entire set of words appearing in WordNet.
	This tool comes with a database containing this entire set.
	It is possible to choose the number of neighbours kept for each word.
	"""
	
	## Initialises the parameters of the module
	def __init__(self, inputFile, outputFile, database, maxRank, discard, verbose):
		## input file for words
		self.inputFile = abspath(inputFile)
		## output file for thesaurus
		self.outputFile = self.inputFile + ".WN" if not outputFile\
			else (abspath(outputFile) if not isdir(outputFile)
			else join(abspath(outputFile), basename(self.inputFile) + ".WN"))
		## fixed maximum rank distance
		self.maxRank = maxRank
		## discard option for words
		self.discard = discard
		## verbose option
		self.verbose = verbose
		## printer to output
		self.printer = outputformatting.Printer(self.verbose)
		## database for similarity scores
		self.database = self.setDatabase(database)

	
	## Determines from the database file potentially chosen by the user whether a
	## database is to be used or not, and if so, which one.
	## Tries the user's choice first, then the default  location, otherwise doesn't use
	## a database.
	## @return name of the database file or None
	def setDatabase(self, chosen):
		final = None
		if chosen and chosen != "...ThesAuto/database.sims":
			try:
				test = open(chosen, 'r')
				test.close()
				final = chosen
			except IOError as e:
				self.printer.info("Can't read specified database file.")
				
		if not final:
			default = join(dirname(realpath(__file__)), ".database.sims")
			try:
				test = open(default, 'r')
				final = default
			except IOError as e:
				self.printer.info("Can't read default database file.")
		return final
	
	
	## Extracts from the input file the words that are to appear in the thesaurus.
	## Each term is the first field of a line.
	## @return list of extracted words
	def extractwords(self):
		words = []
		try:
			for line in open(self.inputFile, 'r'):
				## get 1st field
				t = string.split(line)[0]
				
				if t not in SPECIAL_ENTRIES:
					## keep if in WordNet
					if self.database or wordnet.synsets(t):
						words.append(t)
					else:
						self.printer.info(t + " not found in WordNet")
		except IOError as e:
			self.printer.info("Fatal error: can't read input file.")
			sys.exit()
		return words	
	
	
	## Builds a thesaurus, set by set, using a list of words to find neighbours for and an
	## ic (Information Content of a corpus).
	## Empty sets are removed.
	def buildWordnetThesaurus(self, words, IC):
		## find wanted scores in the database (base thesaurus)
		if self.database:
			baseTh = sorted([line for line in open(self.database)])
			basePt, wordsPt = 0, 0
			
			while basePt < len(baseTh) and wordsPt < len(words):
				tmpSet = baseTh[basePt].split('\t')
				## wanted term found in base: 1) format it (create tuples), 2) move on
				if tmpSet[0] == words[wordsPt]:
					baseTh[basePt] = self.buildNeighbourSetFromDatabase(
						[t for t in words if t != words[wordsPt]], tmpSet[0], tmpSet, k=self.maxRank)
					wordsPt += 1
					basePt += 1
				## unwanted term in the base: delete it
				elif tmpSet[0] < words[wordsPt]: baseTh.remove(baseTh[basePt])
				## term not found in the database: ignore it
				else: wordsPt += 1
			
			return baseTh[:basePt] ## remove unwanted words left in base thesaurus
		
		## compute scores using WordNet
		else: 
			return [set for set in [self.buildNeighbourSetFromWordnet\
				(words, targetTerm, ic=IC, k=self.maxRank) for targetTerm in words] if set]


	## Builds a neighbour set from Wordnet for a target word, using a list of words available
	## and an ic (Information Content of a corpus). This set can be limited to a chosen number
	## of neighbours k.
	## Synsets whose PoS doesn't appear in the ic corpus are ignored.
	## @return neighbour set (type: [string, (string, float), (string, float), ...])
	def buildNeighbourSetFromWordnet(self, words, targetWord, ic, k=None):
		set = []
		for otherWord in [w for w in words if w != targetWord]:
			similarityScore = linWordnetSimilarity(targetWord, otherWord, ic)
			if similarityScore:
				set.append((otherWord, similarityScore))
			
		set.sort(key=itemgetter(1), reverse=True)
		return ([targetWord] + set[:k]) if set[:k] else None


	## Built a neighbour set from the database file for a target word, using a list of words available. 
	## This set can be limited to a chosen number of neighbours k.
	## @return neighbour set (type: [string, (string, float), (string, float), ...])
	def buildNeighbourSetFromDatabase(self, words, targetWord, fields, k=None):
		baseSet = sorted([(fields[i], float(fields[i+1])) for i in xrange(1, len(fields)) if i%2 == 1])
		basePt, wordsPt = 0, 0
		
		while basePt < len(baseSet) and wordsPt < len(words):
			## wanted term found in base: keep it and move on
			if baseSet[basePt][0] == words[wordsPt]:
				wordsPt += 1
				basePt += 1
			## unwanted term in the base: delete it
			elif baseSet[basePt][0] < words[wordsPt]: baseSet.remove(baseSet[basePt])
			## term not found in the database: ignore it
			else: wordsPt += 1
		
		## remove unwanted words left in base set
		baseSet = baseSet[:basePt]
		baseSet.sort(key=itemgetter(1), reverse=True)
		return [fields[0]] + baseSet[:k] if baseSet[:k] else None


	## Converts a neighbour set (containing tuples) to a printable string
	## @return string for the neighbour set
	def setToString(self, set):
		neighbourToString = lambda tuple: tuple[0] + '\t' + str(tuple[1])
		return '\t'.join([set[0]] + [neighbourToString(w) for w in set[1:]]) + '\n'

	
	## Runs the automatic construction of a thesaurus
	def run(self):
		## start operations
		stime = time.time()
		self.printer.mainTitle("Thesauto - automated creation of a thesaurus using WordNet")
		
		## read file
		self.printer.stage(1, 4, "Extracting words from input file")
		self.printer.info("File: " + basename(self.inputFile))
		words = self.extractwords()
		self.printer.lines(words, max=20, 
			title="-- Extracted " + str(len(words)) + " words --")
		
		## prepare WordNet IC
		self.printer.stage(2, 4, "Preparing WordNet IC (Information Content)")
		if not self.database:
			IC = wordnet.ic(brown, False, 1.0)
		else:
			IC = None
			self.printer.info("Using base WordNet thesaurus as a database instead. Skipped.")
			self.printer.info("Base: " + basename(self.database))
			
		## create a thesaurus for each set of words
		self.printer.stage(3, 4, "Building thesaurus")
		thesaurus = self.buildWordnetThesaurus(words, IC)
		self.printer.lines(thesaurus, max=10, line_max=75,
			title="-- Built thesaurus --")
			
		## save the final thesaurus
		self.printer.stage(4, 4, "Saving full thesaurus")
		open(self.outputFile, 'w').write( \
			''.join([self.setToString(set) for set in thesaurus]))
		self.printer.info(self.outputFile + " written.")
		
		etime = time.time()
		self.printer.info("Execution took " + str(etime-stime) + " seconds\n")
	

if __name__=='__main__':
	## parse command line
	parser = argparse.ArgumentParser(description='Compare two thesauri.')
	
	# input file for words
	parser.add_argument(metavar='file', dest='inputFile', action='store',
		help='input file for words to build a thesaurus for')
	# output file for thesaurus
	parser.add_argument('-o', '--output', metavar='file', dest='outputFile',
		action='store',
		help='output file where the thesaurus will be written ' +
			'(default: input file name + \".WN\")')
	# database for similarity scores (base thesaurus)
	parser.add_argument('-db', '--database', metavar='file', dest='database',
		action='store', default="...ThesAuto/database.sims",
		help='base thesaurus acting as a database for similarity scores ' +
			'(default: "...ThesAuto/database.sims")')
	# fixed maximum rank distance
	parser.add_argument('-k', '--max-rank', metavar='n', type=positiveInt, dest='maxRank',
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
		
	a = parser.parse_args()
	thesAuto = ThesAuto(a.inputFile, a.outputFile, a.database, a.maxRank, a.discard, a.verbose)
	thesAuto.run()

	