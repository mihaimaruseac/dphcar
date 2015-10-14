import sys

def read_file(fname):
    docs = []
    with open(fname, 'r') as f:
        for line in f:
            line = map(int, line[:-1].split())
            docs.append(line)
    return docs

def update(edges, x, y):
    l = edges.get(x, [])
    if y not in l:
        l.append(y)
    edges[x]=l

def main(fname):
    docs = read_file(fname)
    n = max(map(max, docs))
    t = len(docs)
    graph = {}
    for d in docs:
        for x, y in zip(d, d[1:]):
            update(graph, x, y)
            update(graph, y, x)
    e = sum(map(lambda x: len(graph[x]), graph)) / 2
    print '{} {} {}'.format(n, e, t)
    for x in xrange(1, n+1):
        l = graph.get(x, [])
        print '{} {}'.format(x, ' '.join(map(str, sorted(l))))
    print '--'
    for d in docs:
        print ' '.join(map(str, d))

if __name__ == '__main__':
    main(sys.argv[1])
