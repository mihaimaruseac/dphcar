import sys

d = {}
npd = {}
c = {}
for fname in sys.argv[1:]:
    notDone = False
    with open(fname, "r") as f:
        f.readline()
        exp = '\t'.join(f.readline().split()[2:-1])
        d.setdefault(exp, {})
        npd.setdefault(exp, {})
        c[exp] = c.get(exp, 0) + 1
        for i in xrange(10): f.readline()
        for i in xrange(10):
            line = f.readline()
            if not line:
                notDone = True
                break
            k, v = line.split()[:2]
            d[exp][k] = d[exp].get(k, 0) + int(v)
        if notDone:
            del d[exp]
            del npd[exp]
            continue
        for i in xrange(1): f.readline()
        for i in xrange(10):
            k, v = f.readline().split()[:2]
            npd[exp][k] = npd[exp].get(k, 0) + int(v)

for exp in sorted(d.keys()):
    print "Experiment", exp
    s = sum(d[exp].values()) + 0.0
    v = 0.0
    npv = 0.0
    cnt = 0.0 + c[exp]
    s = s/c[exp]
    expected = int(exp.split()[-1])
    for k in sorted(d[exp].keys(), reverse=True):
        v = v + d[exp][k]
        npv = npv + npd[exp][k]
        vv = v/cnt
        #npvv = npv/cnt
        #vnpv = vv / npvv if npv else float('nan')
        vs = vv / s if s else float('nan')
        ve = vv / expected if expected else float('nan')
        #print '\t{}\t{:>6.2f}\t{:>10.2f}\t{:3.2f}\t{:3.2f}'.format(k, vv, npvv, vs, vnpv)
        print '\t{}\t{:>6.2f}\t{:3.2f}\t{:3.2f}'.format(k, vv, vs, ve)
    print
