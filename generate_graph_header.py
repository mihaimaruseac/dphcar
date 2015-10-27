import sys

def read_file(fname):
    docs = []
    with open(fname, 'r') as f:
        for line in f:
            line = map(int, line[:-1].split())
            docs.append(filter(lambda x: x > 0, line))
    return docs

def update(edges, x, y):
    l = edges.get(x, [])
    if y not in l:
        l.append(y)
    edges[x]=l

def fix_graph(n, graph, docs):
    actual = min(len(graph), 255)
    mapping = {}
    for i in sorted(graph):
        mapping[i] = len(mapping) + 1
        if len(mapping) == actual:
            break

    docs = [[mapping[x] for x in d if mapping.has_key(x)] for d in docs]
    graph = dict([(mapping[x], [mapping[z] for z in y if mapping.has_key(z)])
            for x, y in graph.iteritems() if mapping.has_key(x)])

    return actual, graph, docs

def main(fname):
    docs = read_file(fname)
    n = max(map(max, docs))
    t = len(docs)
    graph = {}
    for d in docs:
        for x, y in zip(d, d[1:]):
            update(graph, x, y)
            update(graph, y, x)
    #n, graph, docs = fix_graph(n, graph, docs)
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
