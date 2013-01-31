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
from heapq import *
from math import ceil
from nltk.corpus import wordnet, brown
from operator import itemgetter

import printer
	
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
		self.inputFile = abspath(inputFile)
		## output file for thesaurus
		self.outputFile = abspath(outputFile) if not isdir(outputFile)\
			else join(abspath(outputFile), "wordnet.thesaurus")
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
	def extractTerms(self, file, database=".wordnet.words", discard=True, verbose=False):
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
					else:
						printer.info(t+" not found in WordNet", verbose)
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
	def buildPosTermsDictionary(self, terms):
		posTerms = {'a':[], 'v':[], 'n':[], 'r':[]}
		for t in terms:
			for p in posTerms.iterkeys():
				if wordnet.synsets(t, pos=p):
					posTerms[p].append(t)
		posTerms = dict((p, terms) for p, terms in posTerms.items() if len(terms) >= 2)
		return posTerms
	
	##
	##
	def buildWordnetFullThesaurus(self, words, ic):
		## prepare process handling
		one=time.time()
		
		processes, fullThesaurus = [], []
		nbProcesses = min(self.nbCPU, len(words))
		nbFullChunks = len(words) % nbProcesses
		if nbFullChunks == 0: nbFullChunks = nbProcesses
		chunkSize = int( ceil( 1.0*len(words)/nbProcesses ) )
		resultQueue = multiprocessing.Queue()
		#~ print "nbProcess", nbProcesses
		#~ print "nbFullChunks", nbFullChunks
		#~ print "chunkSize", chunkSize
		
		## compute thesaurus chunks in different processes
		for counter in range(nbProcesses):
			start, end = self.determineChunkSlice(nbFullChunks, counter, chunkSize)
			chunkWords = words[start:end]
			p = multiprocessing.Process(name="p"+str(counter), target=self.buildWordnetPartialThesaurus, 
				args=(resultQueue, chunkWords, words, ic, counter))
			processes.append(p)
			p.start()
			#~ print p
			
		for p in processes:
			fullThesaurus += resultQueue.get()
			p.join()
		fullThesaurus.sort()
		
		two=time.time()
		#~ print "Finished creating thesaurus in", two-one, "seconds!"
			
		return fullThesaurus
	
	
	## faster than using % ?
	##
	def determineChunkSlice(self, nbFullChunks, counter, chunkSize):
		start = min(counter, nbFullChunks) * chunkSize
		end = start + chunkSize
		if counter >= nbFullChunks:
			start += (counter - nbFullChunks) * (chunkSize - 1)
			end = start + (chunkSize - 1)
		return (start, end)
		
	
	##
	##
	def buildWordnetPartialThesaurus(self, queue, selectedWords, allWords, ic, nb):
		#~ print "-> starting  for part", nb, "with", len(selectedWords), "words"
		partialThesaurus =  [set for set in [self.buildWordnetNeighbourSet( \
					allWords, targetWord, ic, self.maxRank) \
					for targetWord in selectedWords] if set]
		queue.put(partialThesaurus)
		
		#~ print "<- ending process", nb


	## Builds a neighbour set for a target word using a list of words available and an ic 
	## (Information Content of a corpus). This set can be limited to a chosen number of 
	## neighbours k.
	## Synsets whose PoS doesn't appear in the ic corpus are ignored.
	## @return neighbour set (type: [string, (string, float), (string, float), ...]
	def buildWordnetNeighbourSet(self, words, targetWord, ic, k=None, verbose=False):
		## neighbour set creation
		#~ if verbose:
			#~ print "Creating set for " + targetWord

		set = []
		for otherWord in [w for w in words if w != targetWord]:
			similarityScore = linWordnetSimilarity(targetWord, otherWord, ic)
			if similarityScore:
				set.append((otherWord, similarityScore))
			
		set.sort(key=itemgetter(1), reverse=True)
		
		return ([targetWord] + set[:k]) if set else None
	

	def lookForItem(self, what, attr, where):
		for i, thing in enumerate(where):
			if thing[attr] == what[attr]:
				return i
		return -1
	
	
	## @return set resulting from the merge
	def mergeNeighbourSets(self, neighbours1, neighbours2, maxRank=None):
		pt1, pt2 = 0, 0
		target, resultSet = neighbours1[0], []
		neighbours1, neighbours2 = neighbours1[1:], neighbours2[1:]
		len1, len2 = len(neighbours1), len(neighbours2)
		
		# produce set of duplicated items
		l1, l2 = [n[0] for n in neighbours1], [n[0] for n in neighbours2]
		#~ print "\n** LIST L1 **\n", l1
		#~ print "** LIST L2 **\n", l2
		commonNeighbours = set(l1).intersection(l2)
		#~ print "** COMMON NEIGHBOURS **\n", commonNeighbours, "\n"
		
		while (not maxRank or len(resultSet) < maxRank) and (pt1 < len1 or pt2 < len2) :
			tuple1 = neighbours1[pt1] if pt1 < len1 else ("", -1)
			tuple2 = neighbours2[pt2] if pt2 < len2 else ("", -1)
			
			if tuple1[1] >= tuple2[1]:
				greater = tuple1
				pt1 += 1
			else:
				greater = tuple2
				pt2 += 1
			
			if greater[0] in commonNeighbours: # if is one of the duplicated terms
				i = self.lookForItem(greater, 0, resultSet)
				
				# no duplicate yet
				if i == -1:
					#~ print greater[0], "\t: adding first instance of duplicate"
					resultSet.append(greater)
				# duplicate found
				else:
					commonNeighbours.remove(greater[0])
					# with lower score
					if resultSet[i][1] < greater[1]:
						#~ print greater[0], ": changing", resultSet[i][1] , "for", greater[1]
						resultSet[i] = greater
					# with higher score
					#~ else :
						#~ print greater[0], ": keeping", resultSet[i][1] , "over", greater[1]
			else:
				#~ print greater[0], "\t: adding unique term"
				resultSet.append(greater)
		#~ print len(neighbours1), len(neighbours2), len(resultSet)
		return [target] + resultSet
	

	##
	## @return thesaurus resulting from merge
	def mergeThesauri(self, posThesauri, maxRank=None):
		
		# initialisation
		thHeads = []
		for p in posThesauri.iterkeys():
			thHeads.append( (posThesauri[p][0], p) )
			posThesauri[p] = posThesauri[p][1:]
		heapify(thHeads)
		
		# processing
		resultTh = []
		while thHeads:
			next = heappop(thHeads)
			set, p = next[0], next[1]
			
			# add set to thesaurus
			# if same term, merge them
			if resultTh != [] and resultTh[-1][0] == set[0]:
				set = self.mergeNeighbourSets(resultTh[-1], set, maxRank)
				resultTh[-1] = set
			else:
				resultTh.append(set)
			
			# move next to heap if any left
			if posThesauri[p]:
				heappush(thHeads, (posThesauri[p][0], p) )
				posThesauri[p] = posThesauri[p][1:]
		
		return resultTh


	## Converts a neighbour set (containing tuples) to a printable string
	## @return string for the neighbour set
	def setToString(self, set):
		neighbourToString = lambda tuple: tuple[0] + '\t' + str(tuple[1])
		return '\t'.join([set[0]] + [neighbourToString(w) for w in set[1:]]) + '\n'

	
	## Runs the automatic construction of a thesaurus
	def run(self):
		## start operations
		stime = time.time()
		printer.mainTitle("Thesauto - automated creation of a thesaurus using WordNet", self.verbose)
		
		## read file
		printer.stage(1, 5, "Extracting words from input file")
		terms = self.extractTerms(self.inputFile, self.database, self.discard, self.verbose)
		printer.lines(terms, max=20, 
			title="-- Extracted " + str(len(terms)) + " terms --")
		
		###############################################################################################
		# NOW, DISTRIBUTE TERMS IN POS DICTIONARY AND EXECUTE STEP 2 FOR ALL DICTIONARIES SEPARATELY
		# THEN, USE MERGING OPERATIONS DEVELOPED DURING EXPERIMENT 2 TO REbuild THE COMPLETE THESAURUS
		###############################################################################################
		TESTCUT  =  -1
		###############################################################################################
		
		## sort terms depending on PoS (with duplicates)
		printer.stage(2, 5, "Creating PoS dictionary")
		posTerms = self.buildPosTermsDictionary(terms)		
		for p in posTerms.iterkeys():
			###########################################################
			#~ posTerms[p] = posTerms[p][:TESTCUT]
			###########################################################
			printer.lines(posTerms[p], max=10, line_max=75, \
				title="\n\n** DICTIONARY for PoS "+p+" with "+str(len(posTerms[p]))+" items **")
		
		
		## create a thesaurus for each set of terms
		printer.stage(3, 5, "Creating PoS thesauri")
		brownIC = wordnet.ic(brown, False, 1.0)
		posThesauri = {}
		for p in posTerms.iterkeys():
			start = time.time()
			thesaurus = self.buildWordnetFullThesaurus(posTerms[p], brownIC)
			if thesaurus: posThesauri[p] = thesaurus
			end = time.time()
			
			# display
			printer.lines(posThesauri[p], max=10, line_max=75, \
				title="\n\n-- THESAURUS for PoS "+p+" in "+str(end-start)+" seconds --",
				verbose=self.verbose)
		
		
		printer.stage(4, 5, "Merging PoS thesauri")
		start = time.time()
		resultTh = self.mergeThesauri(posThesauri, self.maxRank)
		end = time.time()
		printer.lines(resultTh, max=10, line_max=75, \
			title="\n\n-- THESAURUS for ALL PoS rebuilt in "+str(end-start)+" seconds --",
			verbose=self.verbose)
		
		
		## save the final thesaurus
		printer.stage(5, 5, "Saving full thesaurus")
		open(self.outputFile, 'w').write( \
			''.join([self.setToString(set) for set in resultTh]))
		printer.info(self.outputFile + " written.", self.verbose)
		
		etime = time.time()
		printer.info("Execution took " + str(etime-stime) + " seconds", self.verbose)



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
	

	