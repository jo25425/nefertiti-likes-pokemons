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
__version__ = "1.0.0"
__maintainer__ = "Joanne Robert"
__email__ = "jr317@sussex.ac.uk"
__status__ = "Development"
#~ -----------------------------------------------------------------------------

import nltk, os, sys
from nltk.corpus import wordnet as wn
from nltk.corpus import wordnet_ic as wn_ic
from operator import itemgetter
import time, string, argparse


## Extracts from the specified files the term that are to appear in the thesaurus.
## Each term is the first field of a line.
## @return list of extracted terms
def extract_terms(file, verbose=False):
	terms = []
	for line in file:
		## get 1st field
		t = string.split(line, "\t")[0]
		## keep if in WordNet and in IC
		if wn.synsets(t):
			terms.append(t)
		elif verbose:
			print t, "not found in WordNet"
			
	return terms


## Creates a neughbour set for a target word using a list of words available and an ic (Information Content
## of a corpus). This set can be limited to a chosen number of neighbours k.
## Synsets whose PoS doesn't appear in the ic corpus are ignored.
## @return neighbour set [type: list of string and (string, float) tuples]
def create_wordnet_neighbour_set(words, targetWord, ic, k=None, verbose=False):
	
	## Lin's WordNet similarity function
	lin_wordnet_sim = lambda word1, word2 : max([wn.lin_similarity(synset1, synset2, ic) \
				for synset1 in wn.synsets(word1) \
				for synset2 in wn.synsets(word2) \
				if synset1.pos == synset2.pos \
				and synset1.pos in ic and synset2.pos in ic])
	
	## neighbour set creation
	if verbose:
		print "Creating set for " + targetWord
	words.remove(targetWord)
	set = [(w, lin_wordnet_sim(targetWord, w)) for w in words]
	set.sort(key=itemgetter(1), reverse=True)
	if k:
		set = set[:k]
	
	return [targetWord] + set


## Converts a neighbour set (containing tuples) to a printable string
## @return string for the neighbour set
def set_to_string(set):
	neighbour_to_string = lambda tuple: tuple[0] + '\t' + str(tuple[1])
	return '\t'.join([set[0]] + [neighbour_to_string(w) for w in set[1:]]) + '\n'


## Prints an array with customizable start, end, line length and title
def print_lines(list, min=0, max=None, line_max=None, title="List"):
	if max is None:
		max = len(list)
	print title
	for index in range(min, max):
		lineText = str(list[index])
		print (lineText[:line_max] + "...") if line_max and len(lineText) > line_max \
			else lineText
	print "\n"


if __name__=='__main__':

	## PARSE COMMAND LINE
	parser = argparse.ArgumentParser(description='Compare two thesauri.')
	
	# input file for words
	parser.add_argument('inputFile', metavar='file', type=file, \
		action='store', help='file for the terms to look for in WordNet')
	# output file for thesaurus
	parser.add_argument('-o', '--output-file', metavar='file', dest='outputFile', 
		action='store', default="./wordnet.thesaurus",
		help='output file in which the thesaurus will be written'+\
		'(default: "./wordnet.thesaurus")')
	# fixed maximum rank distance
	parser.add_argument('-k', '--max-rank', type=int, dest='k', \
		action='store', default=None, \
		help="number of neighbours to keep")
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_const', \
		const=True, default=False, \
		help="display information about operations")
	
	a = parser.parse_args()
	
	## start operations
	stime = time.time()
	print "***************************************************************************"
	print "WORDNET THESAURUS CREATION TOOL"
	print "***************************************************************************\n"

	## prepare WordNet IC
	#~ brown_ic = wn_ic.ic('ic-brown.dat')
	semcor_ic = wn_ic.ic('ic-semcor.dat')
	
	## read file
	print "[1] Extracting words from input file",\
	      "\n------------------------------------"
	terms = extract_terms(a.inputFile,  a.verbose)
	if a.verbose:
		print_lines(terms, title="-- Extracted terms --")
	
	## create thesaurus
	print "\n[2] Creating thesaurus in output file",\
	      "\n-------------------------------------"
	with open(a.outputFile, 'w') as output:
		for t in terms:
			newSet = create_wordnet_neighbour_set(terms, t, semcor_ic, a.k, a.verbose)
			output.write( set_to_string(newSet) )
	
	etime = time.time()
	print "\n>Execution took", etime-stime, "seconds"   
	