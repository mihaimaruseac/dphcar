import itertools
import sys
import matplotlib.pyplot as plt

def color(db):
    if db == 'retail': return 'r'
    if db == 'mushroom': return 'g'
    if db == 'pumsb_star': return 'b'
    return 'k'

def marker(md):
    if md == '3': return 'o'
    if md == '5': return 's'
    return '_'

mx=my=None
with open(sys.argv[1], 'r') as f:
    f.readline() # remove header
    for line in f:
        db, _, _, _, rm, ni, bf, s, _, _, _, r, r50, _\
                = line[:-1].split(',')
        plt.plot(int(r), int(r50), c=color(db.strip()), mfc=color(db.strip()),
                marker=marker(rm), hold=True)
        mx = max(mx, int(r))
        my = max(my, int(r50))

xm = max(mx, my)*1.1
for c in [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
    plt.plot([0, xm], [0, xm * c], '-')
plt.xlabel("# rules")
plt.ylabel("# good rules")
plt.show()
