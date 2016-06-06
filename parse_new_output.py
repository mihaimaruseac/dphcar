import sys

class Exp:
    def __init__(self, seed, n, fr, docs, mintl, maxtl, eps, epss, alpha, mis,
            nt, k):
        self.seed=seed
        self.n=n
        self.fr=fr
        self.docs=docs
        self.mintl=mintl
        self.maxtl=maxtl
        self.eps=eps
        self.epss=epss
        self.alpha=alpha
        self.mis=mis
        self.nt=nt
        self.k=k
        self.rules={}

    def record(self, c, n):
        self.rules[c] = n

    def __str__(self):
        s = '{0.n};{0.fr};{0.docs};{0.mintl};{0.maxtl};{0.eps};'.format(self)
        s += '{0.epss};{0.alpha};{0.mis};{0.nt};{0.k};'.format(self)
        for k in sorted(self.rules, reverse=True):
            s += '{:3.2f};'.format(self.rules[k])
        return s

def die(fname, line):
    print "Invalid line |{}| in file {}".format(line[:-1], fname)
    sys.exit(-1)

for fname in sys.argv[1:]:
    with open(fname, "r") as f:
        line = f.readline()
        if not line.startswith('Called with'): die(fname, line)
        line = f.readline()
        if not line.startswith('./dph'): die(fname, line)
        parts = line.split('_')
        if not parts[1].startswith('seed'): die(fname, line)
        seed=int(parts[1][4:])
        if not parts[2].startswith('n'): die(fname, line)
        n=int(parts[2][1:])
        if not parts[3].startswith('fr'): die(fname, line)
        fr=int(parts[3][2:])/100.0
        if not parts[4].startswith('docs'): die(fname, line)
        docs=int(parts[4][4:-1])*1000
        mintl=int(parts[5])
        parts=parts[-1].split()
        maxtl=int(parts[0])
        eps=float(parts[1])
        epss=float(parts[2])
        alpha=int(parts[3])
        mis=int(parts[4])
        nt=int(parts[5])
        k=int(parts[6])
        exp = Exp(seed, n, fr, docs, mintl, maxtl, eps, epss, alpha, mis, nt, k)
        line = f.readline()
        if not line.startswith('Reading graph'): die(fname, line)
        line = f.readline()
        if not line.startswith('Reading docs'): die(fname, line)
        line = f.readline()
        if not line.startswith('data-struct'): die(fname, line)
        line = f.readline()
        if not line.startswith('Running'): die(fname, line)
        line = f.readline()
        if not line.startswith('Step 1'): die(fname, line)
        line = f.readline()
        if not line.startswith('Step 2'): die(fname, line)
        line = f.readline()
        if not line.startswith('Rules saved'): die(fname, line)
        line = f.readline()
        if not line.startswith('Total time'): die(fname, line)
        line = f.readline() #ignore
        line = f.readline()
        if not line.startswith('Final histogram'): die(fname, line)
        for i in xrange(10):
            line = f.readline()
            if not line.startswith('\t'): die(fname, line)
            parts = line.split()
            exp.record(parts[0], float(parts[3])/k)
        print exp
