# test using ngrams paper
import sys
import itertools
import random
import numpy

import ngrams.lib.NGramSet
import ngrams.Sanitizer

def read_file(fname):
    docs = []
    with open(fname, 'r') as f:
        n, e, t = map(int, f.readline().split())
        readingDocs = False
        while True:
            line = f.readline()
            if not line: break
            line = line[:-1]
            if not line: continue
            if line == '--': readingDocs = True
            elif readingDocs:
                docs.append(map(int, line.split()))
    lmax = max(map(len, docs))
    return n, e, t, lmax, docs

def answer_query(ngramset, query):
    return ngramset[u''.join(map(chr, map(lambda x:x-1, query)))]

def dump_histogram(rules):
    cs = map(lambda x: x/10.0, range(10))
    h = {}
    for c in cs:
        h[c] = 0

    for r in rules:
        rc = rules[r][0]
        for c in cs:
            if rc >= c:
                h[c] += 1

    for c in sorted(cs, reverse=True):
        print '{:5.2f}\t{}'.format(c, h[c])

def main(fname, epsilon, rlmax, k, seed=42):
    print 'Running ngrams model on {} with budget {:5.2f}.'.format(fname, epsilon)
    numpy.random.seed(seed)
    n, e, t, lmax, dataset = read_file(fname)
    ngramset = ngrams.lib.NGramSet.NGramSet(lmax, n)
    ngramset.parse_sequences_memory(dataset)
    ngramset_np = ngrams.lib.NGramSet.NGramSet(lmax, n)
    ngramset_np.parse_sequences_memory(dataset)
    tree, ngramset = ngrams.Sanitizer.ngram(ngramset, n, budget=epsilon,
            sensitivity=lmax)

    top_rules = {}
    for rl in xrange(2, rlmax + 1):
        for ab in itertools.product(xrange(1, n+1), repeat=rl):
            for al in xrange(1, rl):
                a = ab[:al]
                x = answer_query(ngramset, a)
                if  not x: continue
                y = answer_query(ngramset, ab)
                x1 = answer_query(ngramset_np, a)
                y1 = answer_query(ngramset_np, ab)
                top_rules[(a, ab)] = (y / x, (y1 + 0.0)/ x1 if x1 else float('nan'))
                if len(top_rules) > k:
                    del top_rules[min(top_rules, key=lambda x: top_rules[x][0])]
    dump_histogram(top_rules)

if __name__ == '__main__':
    fname = sys.argv[1]
    epsilon = float(sys.argv[2])
    rlmax = int(sys.argv[3])
    k = int(sys.argv[4])
    if len(sys.argv) > 5:
        seed = int(sys.argv[5])
    else:
        seed = 42
    main(fname, epsilon, rlmax, k, seed)
