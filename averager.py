import sys

d = {}
for fname in sys.argv[1:]:
    with open(fname, "r") as f:
        f.readline()
        exp = '\t'.join(f.readline().split()[2:-1])
        d.setdefault(exp, {})
        for i in xrange(8): f.readline()
        for i in xrange(10):
            k, v = f.readline().split()[:2]
            d[exp][k] = d[exp].get(k, 0) + int(v)

for exp in sorted(d.keys()):
    print "Experiment", exp
    s = sum(d[exp].values()) + 0.0
    v = 0
    for k in sorted(d[exp].keys(), reverse=True):
        v = v + d[exp][k]
        print '\t{}\t{:>10}\t{:3.2f}\t{:>12}\t{:3.2f}'.format(k, d[exp][k], d[exp][k] / s, v, v / s)
    print
