import sys

dataset_dir = sys.argv[1]
prec_fname = sys.argv[2]
rule_fname = sys.argv[3]
cuts = ['10', '5', '1', 'no_cut']
if not all(map(lambda x: x.endswith('_avg.csv'), [prec_fname, rule_fname])):
    raise Exception("Use average files")
with open(prec_fname, 'r') as f, open(rule_fname, 'r') as g:
    next(f)
    next(g)
    for l1, l2 in zip(f, g):
        l1 = l1[:-1]
        l2 = l2[:-1]
        parts1 = l1.split(',')
        parts2 = l2.split(',')
        db1, lmax1, k1, eps1, seed1 = parts1[:5]
        db2, lmax2, k2, eps2, seed2 = parts2[:5]
        if any([db1 != db2,
            int(lmax1) != int(lmax2), int(k1) != int(k2),
            float(eps1) != float(eps2),
            seed1 != seed2]):
            raise Exception('Line mismatch\n\t{}\n\t{}'.format(l1, l2))
        db = db1 # = db2
        db_part_name = '{}/{}'.format(dataset_dir, db)
        for end in cuts:
            with open('{}_{}'.format(db_part_name, end), 'r') as h:
                for line in h:
                    if line.startswith('Final non-private'):
                        header = [db1, lmax1, k1, eps1, seed1, end]
                        precision = map(float, parts1[5:])
                        got = map(float, parts2[5:])
                        real = [int(next(h).split()[-2]) for x in range(10)]
                        recall = map(lambda x, y: x/y if y else 1, got, real)
                        f1 = map(lambda x, y: 2 * x * y / (x + y) if x + y else 0, precision, recall)
                        print '{},{},{},{}'.format(
                                ','.join(header),
                                ','.join(map(lambda x: '{:5.2f}'.format(x), precision)),
                                ','.join(map(lambda x: '{:5.2f}'.format(x), recall)),
                                ','.join(map(lambda x: '{:5.2f}'.format(x), f1)))
                        break # stop reading more from this file
