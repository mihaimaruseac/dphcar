import sys

fname = sys.argv[1]
with open(fname, 'r') as f:
    for line in f:
        if line.startswith('Precision'):
            db, lmax, k, eps, seed = fname.split('/')[-1].split('_')
            print '{},{},{},{},{},{}'.format(db, lmax, k, eps, seed,
                    ','.join(line.split()[1:]))
