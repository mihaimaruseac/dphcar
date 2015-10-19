# test using ngrams paper
import sys
import itertools
import random
import numpy
import progressbar

import ngrams.lib.NGramSet
import ngrams.lib.ProgressBar
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
                docs.append(filter(lambda x: x > 0, map(int, line.split())))
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

def main(fname):
    non_private_path = fname + '_non_private'
    n, e, t, lmax, dataset = read_file(fname)
    ngramset = ngrams.lib.NGramSet.NGramSet(lmax, n)
    ngramset.parse_sequences_memory(dataset)
    ngramset.dump(non_private_path)
    print n, lmax

def mine(fname, n, lmax, rlmax, k, epsilon, seed):
    print 'Running ngrams model on {} with budget {:5.2f}.'.format(fname, epsilon)
    numpy.random.seed(seed)
    non_private_path = fname + '_non_private'
    ngramset = ngrams.lib.NGramSet.NGramSet(lmax, n)
    ngramset.load_dump(non_private_path)
    ngramset_np = ngrams.lib.NGramSet.NGramSet(lmax, n)
    ngramset_np.load_dump(non_private_path)
    tree, ngramset = ngrams.Sanitizer.ngram(ngramset, n, budget=epsilon,
            sensitivity=lmax)

    top_rules = {}
    for rl in xrange(2, rlmax + 1):
        pb = ngrams.lib.ProgressBar.MyProgressBar(
                'Mining length {}'.format(rl), pow((n + 1), rl))
        cnt = 0
        for ab in itertools.product(xrange(1, n+1), repeat=rl):
            for al in xrange(1, rl):
                try:
                    a = ab[:al]
                    x = answer_query(ngramset, a)
                    if  not x: continue
                    y = answer_query(ngramset, ab)
                    x1 = answer_query(ngramset_np, a)
                    y1 = answer_query(ngramset_np, ab)
                    top_rules[(a, ab)] = (y / x, (y1 + 0.0)/ x1 if x1 else float('nan'))
                    if len(top_rules) > k:
                        del top_rules[min(top_rules, key=lambda x: top_rules[x][0])]
                except UnicodeDecodeError:
                    continue
                except ValueError:
                    continue
            cnt += 1
            pb.update(cnt)
        pb.finish()

    dump_histogram(top_rules)

building = 0
if __name__ == '__main__':
    if building:
        main(sys.argv[1])
    else:
        fname = sys.argv[1]
        n = int(sys.argv[2])
        lmax = int(sys.argv[3])
        rlmax = int(sys.argv[4])
        k = int(sys.argv[5])
        epsilon = float(sys.argv[6])
        if len(sys.argv) > 7:
            seed = int(sys.argv[7])
        else:
            seed = 42
        mine(fname, n, lmax, rlmax, k, epsilon, seed)
