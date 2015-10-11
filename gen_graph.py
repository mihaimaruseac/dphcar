# params:
#   seed        - seed for rng
#   n           - number of nodes
#   fr          - fill ratio - number of edges 1 = complete graph ($K_n$), 0 = empty graph
#   num_docs    - number of docs to generate
#   min_doc_len - min doc len
#   max_doc_len - max doc len
import itertools
import random
import sys

def update(edges, x, y):
    l = edges.get(x, [])
    l.append(y)
    edges[x]=l

def init_args(seed, n, fr, num_docs, min_doc_len, max_doc_len):
    # type conversions
    seed=int(seed)
    n = int(n)
    fr = float(fr)
    num_docs = int(num_docs)
    min_doc_len=int(min_doc_len)
    max_doc_len=int(max_doc_len)

    # init random seed
    random.seed(seed)
    return seed, n, fr, num_docs, min_doc_len, max_doc_len

def save(n, e, num_docs, edges, docs):
    print n, e, num_docs
    for x in edges:
        print x, ' '.join(map(str, edges[x]))
    print '--'
    for d in docs:
        print ' '.join(map(str, d))

def main_random_original(seed=42, n=10, fr=0.5, num_docs=10,
        min_doc_len=3, max_doc_len=6):
    seed, n, fr, num_docs, min_doc_len, max_doc_len = init_args(seed, n, fr,
            num_docs, min_doc_len, max_doc_len)

    # generate graph
    e = int(fr * (n * (n -1) / 2))
    tuples = random.sample([(1+x, 1+y) for x in xrange(n)
                                  for y in xrange(n)
                                  if x < y], e)

    edges = {}
    for x, y in sorted(tuples):
        update(edges, x, y)
        update(edges, y, x)

    # generate next transitions (2gram model for now)
    nexts = {}
    for cn in xrange(1, n + 1):
        e1 = edges.get(cn, [])
        if not e1:
            nexts[cn] = []
        else:
            cnt = min(4, len(e1))
            nexts[cn] = [(random.uniform(0, 0.2), random.choice(e1))
                    for i in xrange(cnt)]

    # generate docs
    docs = []
    for i in xrange(num_docs):
        doclen = random.randint(min_doc_len, max_doc_len)
        cn = random.randint(1, n)
        pth = []
        for j in xrange(doclen):
            pth.append(cn)
            if not nexts[cn]: break
            s = random.random()
            l = len(nexts[cn])
            cnnx = None
            for i in xrange(l):
                s -= nexts[cn][i][0]
                if s < 0:
                    cnnx = nexts[cn][i][1]
                    break
            if not cnnx:
                cn = random.choice(edges[cn])
            else:
                cn = cnnx
        docs.append(pth)

    save(n, e, num_docs, edges, docs)

def main_ring(seed=42, n=10, fr=0.5, num_docs=10,
        min_doc_len=3, max_doc_len=6):
    seed, n, fr, num_docs, min_doc_len, max_doc_len = init_args(seed, n, fr,
            num_docs, min_doc_len, max_doc_len)

    # generate graph
    e = n
    edges = {}
    for x, y in [(1 + x, 1 + (x + 1) % n) for x in xrange(n)]:
        update(edges, x, y)
        update(edges, y, x)
    # generate docs
    docs = []
    for i in xrange(num_docs):
        doclen = random.randint(min_doc_len, max_doc_len)
        pth = range(1, doclen + 1)
        docs.append(pth)
    save(n, e, num_docs, edges, docs)

if __name__ == '__main__':
    #main=main_random_original
    main=main_ring
    main(**dict(map(lambda x: x.split('='), sys.argv[1:])))
