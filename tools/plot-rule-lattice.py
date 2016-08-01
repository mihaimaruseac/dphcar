import sys
import itertools
import matplotlib.pyplot as plt

gScale = 10000
assert(len(sys.argv) == 2)

supports = {}
icnt = 0
with open(sys.argv[1], "r") as f:
    for line in f:
        parts = line.split()
        items = map(int, parts[:-1])
        support = int(parts[-1][1:-1])
        if len(items) == 1: icnt += 1
        items.sort()
        supports[tuple(items)] = support
assert(2**icnt == len(supports) + 1)

def valid_subsets(s):
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in
            range(1, len(s)))

def mark(la, lab):
    s = 'osph*xd'
    return s[(la-1)%len(s)]

def markcol(la, lab):
    s = 'bgrcmyk'
    return s[(lab-1)%len(s)]

rcnt = 0
for k in supports:
    if len(k) == 1:
        continue # no rule to generate
    for a in valid_subsets(k):
        rcnt += 1
        plt.plot(supports[a], supports[k], hold=True,
                marker=mark(len(a), len(k)),
                markerfacecolor=markcol(len(a), len(k)))
assert(rcnt == 3**icnt - 2**(icnt+1) + 1)

mx = (max(supports.values())/gScale + 1) * gScale
plt.plot([0, mx], [0, mx], "b-", hold=True)
plt.plot([0, mx/0.5], [0, mx], "r-", hold=True)
plt.axis([0, mx, 0, mx])
plt.title(sys.argv[1].split('/')[-1].split('.')[0])
plt.show()
