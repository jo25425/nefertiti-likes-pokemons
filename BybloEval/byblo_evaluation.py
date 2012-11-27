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

import time, string, argparse, math
import nltk, os, sys
from nltk.corpus import wordnet as wn
from nltk.corpus import wordnet_ic as wn_ic
from operator import itemgetter, attrgetter

##
##
def extract_terms(file, method="Lin", k=None):
	terms = []
	for line in file:
		fields = string.split(line, "\t")
		
		# reduce neighbour set if needed
		if k !=None:
			fields = fields[:k*2+1]
		
		# build tuples depending on method specified
		if method == "Lin":
			line_terms = [(fields[i], float(fields[i+1])) for i in xrange(1, len(fields)) if i%2 == 1]
		else:
			line_terms = [(fields[i], i/2) for i in xrange(1, len(fields)) if i%2 == 1]
			
		line_terms.sort()
		line_terms = [fields[0]] + line_terms
		terms.append(line_terms)
	terms.sort()
	return terms


##
##
def create_wordnet_neighbour_set(originalSet, ic, method="Lin"):
	
	## Lin's WordNet similarity function
	lin_wordnet_sim = lambda word1, word2 : max([wn.lin_similarity(synset1, synset2, ic) \
				for synset1 in wn.synsets(word1) \
				for synset2 in wn.synsets(word2) \
				if synset1.pos == synset2.pos \
				and synset1.pos in ic and synset2.pos in ic])
	
	## neighbour set creation
	word = originalSet[0]
	set = [(neighbour, lin_wordnet_sim(word, neighbour)) 
		for (neighbour, val) in originalSet[1:] 
		if wn.synsets(neighbour)]
	
	## conversion to ranks when needed
	if method != "Lin":
		set.sort(key=itemgetter(1), reverse=True)
		set = [(neighbour, rank) for rank, (neighbour, val) in enumerate(set)]
	
	return [word] + sorted(set)


##
##
def create_wordnet_thesaurus(originalThesaurus, method="Lin"):
	
	ic = wn_ic.ic('ic-brown.dat')
	#~ semcor_ic = wn_ic.ic('ic-semcor.dat')
	return [create_wordnet_neighbour_set(set, ic, method) 
		for set in originalThesaurus if wn.synsets(set[0])]


##
##
def neighbour_set_sim_lin(set1, set2, k=None, max_sim_score=None, verbose=False):
	sim_score, pt1, pt2 = 0, 1, 1
	square_score = lambda x: x[1] * x[1]
	
	## determine maximum possible similarity score
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
			
			# what is happening?
			if verbose:
				print "["+str(pt1)+"]"+str(set1[pt1])
				print "["+str(pt2)+"]"+str(set2[pt2])
				print "\t=> "+str(score_product)
			
			pt1 += 1
			pt2 += 1
			
	return 1.0 * sim_score / max_sim_score

##
##
def neighbour_set_sim_rank(set1, set2, k=None, max_sim_score=None, verbose=False):
	sim_score, pt1, pt2 = 0, 1, 1
	
	## determine maximum rank distance if not passed
	if k == None:
		k = max(len(set1), len(set2)) - 1 # target word doesn't count
		
	## determine maximum similarity score if not passed
	if max_sim_score == None:
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
				
			# what is happening?
			if verbose:
				print "["+str(pt1)+"]"+str(set1[pt1])
				print "["+str(pt2)+"]"+str(set2[pt2])
				print "\t=> "+str(score_product)
			
			pt1 += 1
			pt2 += 1
			
	return 1.0 * sim_score / max_sim_score


##
##
def thesaurus_sim(th, method="rank", n=None, k=None, verbose=False):
	th1, th2 = th[0], th[1]
	global_sim_score, pt1, pt2 = 0, 0, 0
	fixed_k = (k != None)
	neighbour_set_sim = {
		'Lin'	: neighbour_set_sim_lin,
		'rank': neighbour_set_sim_rank,
	}[method]
	
	## determine maximum thesaurus size if not passed (n allows comparison on a subset)
	if n == None:
		n = max(len(th1), len(th2))
		
	## avoid division by zero
	if n == 0:
		return 0
	
	## determine a possible maximum rank distance if not passed
	if k == None:
		k = len(th1[0]) - 1
	
	## and a possible maximum neighbour set similarity score
	## (to avoid unneccessary computations of unchanged values)
	max_sim_score = sum([i*i for i in xrange(1, k+1)])
	
	while pt1 < n and pt2 < n and pt1 < len(th1) and pt2 < len(th2) :
		
		if th1[pt1][0] == th2[pt2][0]: # target word present in both sets
			l1, l2 = len(th1[pt1])-1, len(th2[pt2])-1
			
			if fixed_k and k <= l1 and k <= l2 or k == l1 == l2 :
				set_sim = neighbour_set_sim(th1[pt1], th2[pt2], k, max_sim_score)
			else:
				set_sim = neighbour_set_sim(th1[pt1], th2[pt2])
			global_sim_score += set_sim
			
			# what is happening?
			if verbose:
				print "th1[", pt1, "]", th1[pt1][0], ": th2[", pt2, "]", th2[pt2][0], \
				"=> ", set_sim
			
			pt1 += 1
			pt2 += 1
			
		elif th1[pt1][0] < th2[pt2][0]:
			pt1 += 1
		else:
			pt2 += 1
	return 1.0 * global_sim_score / n
	
	
## Prints an array with customizable start, end, line length and title
def print_lines(list, min=0, max=None, line_max=None, title="List"):
	if max is None:
		max = len(list)
	print title + "\n"
	for index in range(min, max):
		print str(list[index])[:line_max], "..."
	print "\n"
	


if __name__=='__main__':

	## PARSE COMMAND LINE
	parser = argparse.ArgumentParser(description='Compare two thesauri.')
	subparsers = parser.add_subparsers()
	
	# subparser for one thesaurus
	oneFileParser = subparsers.add_parser('single')
	oneFileParser.add_argument('files', metavar='file', type=file, nargs=1, \
		action='store', help='file for the thesaurus')
	oneFileParser.add_argument('-w', '--Wordnet', dest='method', action='store_const', \
		const='Wordnet', default='Wordnet', \
		help='method used to compare a thesaurus to WordNet')
	
	#subparser for two thesauri
	twoFilesParser = subparsers.add_parser('pair')
	twoFilesParser.add_argument('files', metavar='file', type=file, nargs=2, \
		action='store', help='files for the thesauri')
	
	
	# output file to write the results in
	parser.add_argument('-o', '--output-file', metavar='file', dest='outputFile', 
		action='store', default="./result.thsim",
		help='output file in which similarity scores and file names will be written'+\
		'(default: "./result.thsim")')
	# neighbour sets similarity measure
	parser.add_argument('-l', '--Lin', dest='method', action='store_const', \
		const='Lin', default='rank', \
		help='method used to calculate the similarity between two neighbour sets')
	# test index
	parser.add_argument('-i', '--test-index', type=int, dest='i', \
		action='store', default=0, \
		help="index of a single word to compare as a test")
	# fixed maximum rank distance
	parser.add_argument('-k', '--max-rank', type=int, dest='k', \
		action='store', default=None, \
		help="maximum number of neighbours to compare")
	# number of lines to compare
	parser.add_argument('-n', '--nb-lines', type=int, dest='n', \
		action='store', default=None, \
		help="maximum number of lines to compare")
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_const', \
		const=True, default=False, \
		help="display information about each comparison")
		

	## PROCESS PARAMETERS
	a = parser.parse_args()
	for item in vars(a):
		print item, ":", vars(a)[item]
	
	## start operations
	stime = time.time()
	print "***************************************************************************"
	print "THESAURI COMPARISON TOOL"
	print "***************************************************************************\n"

	## read files
	thesauri = [extract_terms(file, a.method, a.k) for file in a.files]
	names = [f.name for f in a.files]
	if len(a.files) < 2:
		thesauri.append(create_wordnet_thesaurus(thesauri[0], a.method))
		names.append("WordNet")
	
	## verify extracted information
	for th, name in zip(thesauri, names):
		print_lines(th, max=10, line_max=75, title="\""+name+"\" after sort:")
		print '\n'
	
	## NEIGHBOUR SET COMPARISON
	print "\nComparing sets at index", a.i
	print "*******************************"
	set1, set2 = [th[a.i] for th in thesauri]
	sim_score1 = {
	  'Lin'	: 	neighbour_set_sim_lin(set1, set2, verbose=a.verbose),
	  'rank': 	neighbour_set_sim_rank(set1, set2, verbose=a.verbose)
	}[a.method]
	print "\nSimilarity score = ", sim_score1, "\n"
	
	
	## THESAURI COMPARISON
	print "\nComparing thesauri"
	print "******************"
	with open(a.outputFile, 'w') as output:
		sim_score2 = thesaurus_sim(thesauri, a.method, n=a.n, verbose=a.verbose)
		output.write(str(sim_score2) + '\t' + '\t'.join(names))
	print "\nSimilarity score = ", sim_score2, "\n"

	etime = time.time()
	print "\n>Execution took", etime-stime, "seconds"   
	