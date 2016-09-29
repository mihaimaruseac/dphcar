import sys

def add_to_dict(d, k, v):
    l = d.get(k, [])
    l.append(v)
    d[k] = l

def mean(l):
    s = sum(l)
    n = len(l)
    return (s + 0.0) / n

d_r = {}
d_r50 = {}
d_rr = {}

with open(sys.argv[1], 'r') as f:
    f.readline()
    for line in f:
        db,eps,epsh,c0,rmax,ni,bf,seed,_,_,_,r,r50,_,rr,_ = line.split(",")
        key = db,eps,c0,rmax,ni,bf
        add_to_dict(d_r, key, int(r))
        add_to_dict(d_r50, key, int(r50))
        add_to_dict(d_rr, key, int(rr))

for k in sorted(d_r):
    r = mean(d_r[k])
    r50 = mean(d_r50[k])
    rr = mean(d_rr[k])
    p = r50 / r
    rc = r50 / rr
    print "{}\t{:>11.1f}\t{:>11.1f}\t{:>11.1f}\t{:>5.2f}\t{:>5.2f}".format(
            '\t'.join(k), r, r50, rr, p, rc)
