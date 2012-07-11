
import os, sys, nltk, time, string, argparse, math

#################################################################
def extract_terms(file, method="rank", k=None):
	terms = []
	for line in file:
		fields = string.split(line, "\t")
		
		# reduce neighbor set if needed
		if k !=None:
			fields = fields[:k*2+1]
		
		# build tuples depending on method specified
		if method == "rank":
			line_terms = [(fields[i], i/2) for i in xrange(1, len(fields)) if i%2 == 1]
		else:
			line_terms = [(fields[i], float(fields[i+1])) for i in xrange(1, len(fields)) if i%2 == 1]
			
		line_terms.sort()
		line_terms = [fields[0]] + line_terms
		terms.append(line_terms)
	terms.sort()
	return terms


def neighbor_set_similarity(set1, set2, method="rank", k=None, max_sim_score=None, verbose=False):
	sim_score, pt1, pt2 = 0, 1, 1
	
	if method == "rank":
		# determine maximum rank distance if not passed
		if k == None:
			k = max(len(set1), len(set2)) - 1 # target word doesn't count	
		# determine maximum similarity score if not passed
		if max_sim_score == None:
			max_sim_score = sum([i*i for i in xrange(1, k+1)])
	else:
		# determine maximum similarity score
		max_sim_score = sum([set1[i][1]*set1[i][1] for i in xrange(1, len(set1))]) \
			* sum([set2[i][1]*set2[i][1] for i in xrange(1, len(set2))])
		max_sim_score = math.sqrt(max_sim_score)
		
	while pt1 < len(set1) and pt2 < len(set2) :
		if set1[pt1][0] < set2[pt2][0]:
				pt1 += 1
		else:
			if set1[pt1][0] > set2[pt2][0]:
				pt2 += 1
			
			else:	
				# potential neighbor present in both sets
				if method == "rank":
					score_product = (k - set1[pt1][1])*(k - set2[pt2][1])
				else:
					score_product = set1[pt1][1]*set2[pt2][1]
				sim_score += score_product
				
				# what is happening?
				if verbose:
					print "set1[", pt1, "]", set1[pt1], ": set2[", pt2, "]", set2[pt2],  \
					"=> ", score_product
				
				pt1 += 1
				pt2 += 1
			
	return 1.0 * sim_score / max_sim_score
	
def thesaurus_similarity(th1, th2, n=None, k=None, verbose=False):
	global_sim_score, pt1, pt2 = 0, 0, 0
	fixed_k = (k != None)
	
	# determine maximum thesaurus size if not passed
	# (n allows comparison on a subset)
	if n == None:
		n = max(len(th1), len(th2))
	print "\n# of entries to compare = ", n
	
	# determine a possible maximum rank distance if not passed
	if k == None:
		k = len(th1[0]) - 1
	
	# and a possible maximum neighbor set similarity score
	# (to avoid unneccessary computations of unchanged values)
	max_sim_score = sum([i*i for i in xrange(1, k+1)])
	
	while pt1 < n and pt2 < n and pt1 < len(th1) and pt2 < len(th2) :
		
		if th1[pt1][0] == th2[pt2][0]:
			# target word present in both sets
			l1, l2 = len(th1[pt1])-1, len(th2[pt2])-1
			if fixed_k and k <= l1 and k <= l2 or k == l1 == l2 :
				set_sim = neighbor_set_similarity(th1[pt1], th2[pt2], k, max_sim_score)
			else:
				set_sim = neighbor_set_similarity(th1[pt1], th2[pt2])
			global_sim_score += set_sim
			
			# what is happening?
			if verbose:
				print "th1[", pt1, "]", th1[pt1][0], ": th2[", pt2, "]", th2[pt2][0], \
				"=> ", set_sim
			
			pt1 += 1
			pt2 += 1
		else:
			if th1[pt1][0] < th2[pt2][0]:
				pt1 += 1
			else:
				pt2 += 1
	return 1.0 * global_sim_score / n
	
	
def print_lines(list, min=0, max=None, line_max=None, title="List"):
	if max is None:
		max = len(list)
	print title + "\n"
	for index in range(min, max):
		print str(list[index])[:line_max], "..."
	print "\n"
	
#################################################################

if __name__=='__main__':

	parser = argparse.ArgumentParser(description='Compare two thesauri.')
	# comparison method
	parser.add_argument('-l', '--Lin', dest='method', action='store_const', \
		const='Lin', default='rank', \
		help='method used to calculate the similarity between two neighbor sets')
	# test index
	parser.add_argument('-i', '--test-index', type=int, nargs=1, dest='i', \
		action='store', default=None)
	# fixed maximum rank distance
	parser.add_argument('-k', '--max-rank', type=int, nargs=1, dest='k', \
		action='store', default=None)
	# number of lines to compare
	parser.add_argument('-n', '--nb-lines', type=int, nargs=1, dest='n', \
		action='store', default=None)
	# verbose option
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_const', \
		const=True, default=False)
	# file names
	parser.add_argument('files', metavar='file', type=file, nargs=2, \
		action='store', help='files for the thesauri')
	
	args = parser.parse_args()
	#############################################################
	i = args.i[0] if args.i !=None else 0
	k = args.k[0] if args.k !=None else None
	n = args.n[0] if args.n !=None else None
	
	stime = time.time()
	print "***************************************************************************"
	print "THESAURI COMPARISON TOOL"
	print "***************************************************************************\n"

	# read files
	th1 = extract_terms(args.files[0], args.method, k)
	th2 = extract_terms(args.files[1], args.method, k)
	
	# verify extracted information
	print_lines(th1, max=10, line_max=75, title="\""+args.files[0].name+"\" after sort:")
	print_lines(th2, max=10, line_max=75, title="\""+args.files[1].name+"\" after sort:")
	
	# test for a single target word
	print "\nComparing sets at index", i
	print "*******************************"
	sim_score1 = neighbor_set_similarity(th1[i], th2[i], method=args.method, verbose=args.verbose)
	print "\nSimilarity score = ", sim_score1, "\n"
	
	# and finally compare the pair of thesauri
	print "\nComparing thesauri"
	print "******************"
	sim_score2 = thesaurus_similarity(th1, th2, n, args.verbose)
	print "\nSimilarity score = ", sim_score2, "\n"

	etime = time.time()
	print ">Execution took", etime-stime, "seconds"
	
	# pause at the end
	raw_input()

	    