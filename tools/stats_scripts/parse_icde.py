import os
import sys
import matplotlib.pyplot as plt

class Experiment:
    def __init__(self, db, eps, eps_share, c0, rmax, nits, bf, seed):
        self.db = db
        self.eps = eps
        self.eps_share = eps_share
        self.c0 = c0
        self.rmax = rmax
        self.nits = nits
        self.bf = bf
        self.seed = seed
        self.leaves = 0
        self.time = 0
        self.tpl = 0
        self.rules50 = 0
        self.rules = 0
        self.prec = 0

    def recordLeaves(self, leaves):
        self.leaves = leaves

    def recordTime(self, time):
        self.time = time
        self.tpl = self.time / self.leaves

    def recordRules50(self, rules50):
        self.rules50 = rules50

    def recordPrecision(self, prec):
        self.prec = prec

    def recordRules(self, rules):
        self.rules = rules

    def __str__(self):
        s = "{0.db:10}\t{0.eps}\t{0.eps_share}\t{0.c0}".format(self)
        s = "{1}\t{0.rmax}\t{0.nits}\t{0.bf}\t{0.seed}".format(self, s)
        s = "{1}\t{0.leaves}\t{0.time}\t{0.tpl:6.3f}".format(self, s)
        s = "{1}\t{0.rules}\t{0.rules50}\t{0.prec:5.2f}".format(self, s)
        return s

    @staticmethod
    def legend():
        s = "db\t\teps\tepsh\tc0\trmax\tnits\tbf\tseed\tleaves\ttime"
        s += "\ttpl\trules\trules50\tprec50"
        return s

experiments = []
for r, _, files in os.walk(sys.argv[1]):
    for fname in files:
        with open(os.path.join(r, fname), "r") as f:
            exp = None
            for line in f:
                if line.startswith("./dph"):
                    _, ds, eps, es, c0, rmax, nits, bf, seed = line.split()
                    exp = Experiment(ds.split('/')[-1].split('.')[0],\
                            float(eps), float(es), float(c0), int(rmax),\
                            int(nits), int(bf), int(seed))
                elif line.startswith("Total leaves"):
                    exp.recordLeaves(int((line.split()[-1])))
                elif line.startswith("Total time"):
                    exp.recordTime(float((line.split()[-1])))
                elif line.startswith("\t0.50"):
                    _, _, _, rs, prec = line.split()
                    exp.recordRules50(int(rs))
                    exp.recordPrecision(float(prec))
                elif line.startswith("\t0.00"):
                    _, _, _, rs, _ = line.split()
                    exp.recordRules(int(rs))
            experiments.append(exp)
    break # disable recursion

def print_experiments():
    print Experiment.legend()
    for exp in experiments:
        print exp

def plot_exp(exps, xfun, yfun, selecfuns=[], outname=None, title=None,
        xlabel=None, ylabel=None, xrng=None, yrng=None):
    xmax = ymax = None
    for exp in experiments:
        if not all([sf(exp) for sf in selecfuns]):
            continue
        xmax = max(xmax, xfun(exp))
        ymax = max(ymax, yfun(exp))
        plt.plot(xfun(exp), yfun(exp), 'bo')
    if xmax is None or ymax is None:
        return
    if xrng: plt.xlim(xrng)
    else: plt.xlim([0, 1.1*xmax])
    if yrng: plt.ylim(yrng)
    else: plt.ylim([0, 1.1*ymax])
    if title: plt.title(title)
    if xlabel: plt.xlabel(xlabel)
    if ylabel: plt.ylabel(ylabel)
    if outname: plt.savefig(outname)
    else: plt.show()
    plt.clf()

getDB = lambda xp: xp.db
getRMax = lambda xp: xp.rmax
getNI = lambda xp: xp.nits
getBF = lambda xp: xp.bf
getR = lambda xp: xp.rules
getR50 = lambda xp: xp.rules50
getP50 = lambda xp: xp.prec
getTime = lambda xp: xp.time
getTimeLeaf = lambda xp: xp.tpl

# have the x part and the selector filled in
def plot_against(exps, xfun, selecfuns, xlabel, xrng, xtitle):
    title = 'Total time vs {}'.format(xtitle)
    outname = '{}.png'.format(title.replace('/','_').replace(' ','_'))
    plot_exp(exps, xfun, getTime, selecfuns, xlabel=xlabel, xrng=xrng,
            ylabel="time (s)", yrng=None, title=title, outname=outname)
    title = 'Time per leaf vs {}'.format(xtitle)
    outname = '{}.png'.format(title.replace('/','_').replace(' ','_'))
    plot_exp(exps, xfun, getTimeLeaf, selecfuns, xlabel=xlabel, xrng=xrng,
            ylabel="time (s)", yrng=None, title=title, outname=outname)
    title = 'Total rules vs {}'.format(xtitle)
    outname = '{}.png'.format(title.replace('/','_').replace(' ','_'))
    plot_exp(exps, xfun, getR, selecfuns, xlabel=xlabel, xrng=xrng,
            ylabel="# rules", yrng=None, title=title, outname=outname)
    title = 'Good rules vs {}'.format(xtitle)
    outname = '{}.png'.format(title.replace('/','_').replace(' ','_'))
    plot_exp(exps, xfun, getR50, selecfuns, xlabel=xlabel, xrng=xrng,
            ylabel="# rules", yrng=None, title=title, outname=outname)
    title = 'Precision vs {}'.format(xtitle)
    outname = '{}.png'.format(title.replace('/','_').replace(' ','_'))
    plot_exp(exps, xfun, getP50, selecfuns, xlabel=xlabel, xrng=xrng,
            ylabel="% rules", yrng=[-0.1,1.1], title=title, outname=outname)

dbs = set([getDB(xp) for xp in experiments])
rmaxes = set([getRMax(xp) for xp in experiments])
nis = set([getNI(xp) for xp in experiments])
bfs = set([getBF(xp) for xp in experiments])

for db in dbs:
    dbsf = lambda xp: getDB(xp) == db
    for rmax in rmaxes:
        rmaxsf = lambda xp: getRMax(xp) == rmax
        for ni in nis:
            nisf = lambda xp: getNI(xp) == ni
            plot_against(experiments, getBF, [dbsf, rmaxsf, nisf],
                    "branching factor", None,
                    "branch for {} rmax={} items={}".format(db, rmax, ni))
        for bf in bfs:
            bisf = lambda xp: getBF(xp) == bf
            plot_against(experiments, getNI, [dbsf, rmaxsf, bisf],
                    "# items", None,
                    "items for {} rmax={} branch={}".format(db, rmax, bf))
