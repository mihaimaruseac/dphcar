import itertools
import random
import sys

def update(edges, x, y):
    l = edges.get(x, [])
    l.append(y)
    edges[x]=l

# params:
#   seed        - seed for rng
#   n           - number of nodes
#   fr          - fill ratio - number of edges 1 = complete graph ($K_n$), 0 = empty graph
#   num_docs    - number of docs to generate
#   min_doc_len - min doc len
#   max_doc_len - max doc len
def main(seed=42, n=10, fr=0.5, num_docs=10,
        min_doc_len=3, max_doc_len=6):
    # type conversions
    seed=int(seed)
    n = int(n)
    fr = float(fr)

    # init random seed
    random.seed(seed)

    # generate graph
    e = int(fr * (n * (n -1) / 2))
    tuples = random.sample([(1+x, 1+y) for x in xrange(n)
                                  for y in xrange(n)
                                  if x < y], e)

    edges = {}
    for x, y in sorted(tuples):
        update(edges, x, y)
        update(edges, y, x)

    # generate docs
    docs = []
    for i in xrange(num_docs):
        doclen = random.randint(min_doc_len, max_doc_len)
        cn = random.randint(1, n)
        pth = []
        for j in xrange(doclen):
            pth.append(cn)
            if not edges.get(cn, []): break
            cn = random.choice(edges[cn])
        docs.append(pth)

    # output datafile
    print n, e, num_docs
    for x in edges:
        print x, ' '.join(map(str, edges[x]))
    print '--'
    for d in docs:
        print ' '.join(map(str, d))

if __name__ == '__main__':
    main(**dict(map(lambda x: x.split('='), sys.argv[1:])))
