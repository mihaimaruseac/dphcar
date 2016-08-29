import sys
import itertools
import matplotlib.pyplot as plt

gScale = 1000
c0 = 0.5

assert(len(sys.argv) == 5)
n = int(sys.argv[1])
rmax = int(sys.argv[2])
fname = sys.argv[3]
mode = sys.argv[4]

def valid_subsets(s, m):
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in
            range(1, m + 1))


supports = {}
with open(fname, "r") as f:
    for line in f:
        parts = line.split()
        items = map(int, parts[:-1])
        support = int(parts[-1][1:-1])
        if len(items) > rmax: continue
        items.sort()
        supports[tuple(items)] = support

top_its = []
for k in supports:
    if len(k) == 1:
        top_its.append((k[0], supports[k]))
top_its.sort(cmp=lambda x, y:y[1]-x[1])
top_its = map(lambda x:x[0], top_its)[:n]

def plot_rules_from(k, mf='o', mfc='w', mfca='w'):
    for a in valid_subsets(k, rmax):
        plt.plot(supports[a], supports[k], hold=True,
                marker=mf,
                markerfacecolor=mfc,
                markerfacecoloralt=mfca,
                fillstyle="top" if mfc != mfca else "full")

# plot all rules on background
for k in supports:
    if len(k) == 1:
        continue # no rules to generate
    if any([x not in top_its for x in k]):
        continue # not to be considered now
    plot_rules_from(k)

# abstract method for plot_sigma and plot_delta (strategy pattern)
def plot_delta_sigma(rmax, it, f):
    for i in range(len(it) + 1, rmax + 1):
        print it
        l = [(list(k), f(k))
                for k in supports
                if len(k) == i and all([x in k for x in it])]
        l.sort(cmp=lambda x, y: 1 if y[1] - x[1] > 0 else -1 if y[1] - x[1] < 0 else 0)
        print l[:5], '..'
        it = tuple(l[0][0])
        plot_rules_from(it, mfc='b')
    print ">", it
    return it

# plot rules sampled by q_sigma all the way
def plot_sigma(rmax, it):
    return plot_delta_sigma(rmax, it, lambda k: supports[k])

def plot_delta(rmax, it):
    return plot_delta_sigma(rmax, it,
            lambda k: supports[k] - min([supports[i]
                            for i in valid_subsets(k, rmax) if i != k]))

def plot_d(rmax, it):
    return plot_delta_sigma(rmax, it,
            lambda k: supports[k]/c0 - min([supports[i]
                            for i in valid_subsets(k, rmax) if i != k]))

it = (top_its[0],)
if mode == "sigma":
    plot_sigma(rmax, it)
elif mode == "delta":
    plot_delta(rmax, it)
elif mode == "d":
    plot_d(rmax, it)
elif mode == "sigma-d":
    plot_d(rmax, plot_sigma(rmax, it))
elif mode == "delta-d":
    plot_d(rmax, plot_delta(rmax, it))

mx = (max(supports.values())/gScale + 1) * gScale
plt.plot([0, mx], [0, mx], "b-", hold=True)
plt.plot([0, mx/c0], [0, mx], "r-", hold=True)
plt.axis([0, mx, 0, mx])
fname = '{}_{}_{}_{}'.format(fname.split('/')[-1].split('.')[0], n, rmax, mode)
plt.savefig('{}.png'.format(fname))
