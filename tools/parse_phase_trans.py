import os
import sys
import matplotlib.pyplot as plt

times = {}
rules = {}
precs = {}

def add_to_dict(d, k, v):
    l = d.get(k, [])
    l.append(v)
    d[k] = l

for r, _, files in os.walk(sys.argv[1]):
    for fname in files:
        with open(os.path.join(r, fname), "r") as f:
            k = None
            for line in f:
                if line.startswith("./dph"):
                    _, _, ds, eps, es, c0, rmax, k, nits, rf, _ = line.split()
                    k = ds.split('/')[-1].split('.')[0],\
                        float(eps), float(es), float(c0),\
                        int(rmax), int(k), int(rf), int(nits)
                elif line.startswith("Round"):
                    add_to_dict(times, k, float(line.split()[-1]))
                elif line.startswith("\t0.50"):
                    _, _, _, rs, prec = line.split()
                    add_to_dict(rules, k, int(rs))
                    add_to_dict(precs, k, float(prec))

exps = {}
for k in times:
    k1, ni = k[:-1], k[-1]
    add_to_dict(exps, k1, (ni, times[k], rules[k], precs[k]))

def plot_exp(exps, ix, title, ylabel):
    for k in exps:
        fname = "{}_{}.png".format('_'.join(title.split()), '_'.join(map(str, k)))
        #f = plt.figure()
        plt.title("{} for {}".format(title, k))
        plt.xlabel("num items")
        plt.ylabel(ylabel)
        for it in exps[k]:
            x = it[0]
            ys = it[ix]
            plt.plot([x]*len(ys), ys, 'bo')
        xmax = max([it[0] for it in exps[k]])
        ymax = max([max(it[ix]) for it in exps[k]])
        plt.axis([-0.1, 1.1*xmax, -0.1, 1.1*ymax])
        plt.savefig(fname)
        plt.clf()

plot_exp(exps, 1, "Time", "time (s)")
plot_exp(exps, 2, "Rule count", "# rules")
plot_exp(exps, 3, "Precision", "precision")
