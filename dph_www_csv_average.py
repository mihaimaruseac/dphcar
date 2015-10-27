import sys

d = {}
fname = sys.argv[1]
with open(fname, 'r') as f:
    next(f)
    for line in f:
        line = line[:-1]
        parts = line.split(',')
        db, lmax, k, eps, seed = parts[:5]
        data = [1] + map(float, parts[5:])
        old_line = d.get((db, lmax, k, eps), [0] * 11)
        new_line = map(lambda x,y: x+y, data, old_line)
        d[(db, lmax, k, eps)] = new_line
        print line

for k in sorted(d):
    if d[k][0] != 10:
        continue
    print '{},avg,{}'.format(','.join(k),
            ','.join(map(lambda x: '{:.2f}'.format(x / (d[k][0] + 0.0)), d[k])[1:]))
