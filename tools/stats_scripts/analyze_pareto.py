import itertools
import sys
import matplotlib.pyplot as plt

def color(db):
    if db == 'retail': return 'r'
    if db == 'mushroom': return 'g'
    if db == 'pumsb_star': return 'b'
    return 'k'

def marker(md):
    if md == 'n_s_s': return 'o'
    if md == 'nisis': return 'v'
    if md == 'nasas': return '^'
    if md == 'n_s_d': return '>'
    if md == 'nisid': return '<'
    if md == 'nasad': return '1'
    if md == 'n_e_e': return '2'
    if md == 'nieie': return '3'
    if md == 'naeae': return '4'
    if md == 'n_e_d': return 's'
    if md == 'nieid': return 'p'
    if md == 'naead': return '*'
    if md ==   'n_d': return 'h'
    if md ==   'nid': return 'H'
    if md ==   'nad': return 'D'
    return '_'

mx=my=None
with open(sys.argv[1], 'r') as f:
    f.readline() # remove header
    for line in f:
        db, md, _, _, _, rm, ni, bf, s, _, _, _, r, r50, _, _, _\
                = line[:-1].split(',')
        plt.plot(int(r), int(r50), c=color(db.strip()), mfc=color(db.strip()),
                marker=marker(md), hold=True)
        mx = max(mx, int(r))
        my = max(my, int(r50))

xm = max(mx, my)*1.1
for c in [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
    plt.plot([0, xm], [0, xm * c], '-')
plt.xlabel("# rules")
plt.ylabel("# good rules")
plt.show()
