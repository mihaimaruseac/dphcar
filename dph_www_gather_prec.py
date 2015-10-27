import sys

fname = sys.argv[1]
with open(fname, 'r') as f:
    for line in f:
        if line.startswith('Final private'):
            print '{},{}'.format(','.join(fname.split('/')[-1].split('_')),
                    ','.join([next(f).split()[-2] for x in range(10)]))
            sys.exit(0)
