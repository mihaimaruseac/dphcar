import itertools
import bigfloat
import math

dataset = 'kosarak'

#method = 'mt'     # distance to (m, theta)
method = 'xy'     # linear biproduct of two triangles
# DO NOT USE THESE
#method = 'tt'     # distance to (theta, theta)
#method = 'aab'    # sup(A) * sup(AB)
#method = 'amab'   # sup(A) * (m - sup(AB))
#method = 'aabtab' # sup(A) * sup(AB) * (theta - sup(AB))
#method = 'exp'    # sup(A) * exp(- |theta - sup(AB)|)

if dataset == 'kosarak':
    itemsets = {
     (6,): 601374, (3,): 450031, (11,): 364065, (1,): 197522, (218,): 88598,
     (3, 6): 265180, (6, 11): 324013, (1, 6): 132113, (6, 218): 77675,
     (3, 11): 161286, (1, 3): 84660, (3, 218): 38673, (1, 11): 91882,
     (11, 218): 61656, (1, 218): 32127,
     (3, 6, 11): 143682, (1, 3, 6): 57802, (3, 6, 218): 33978, (1, 6, 11): 86092,
     (6, 11, 218): 60630, (1, 6, 218): 30182, (1, 3, 11): 40268,
     (3, 11, 218): 26973, (1, 11, 218): 25790, (1, 3, 218): 13928,
     (1, 3, 6, 11): 37743, (3, 6, 11, 218): 26501, (1, 6, 11, 218): 25608,
     (1, 3, 11, 218): 11176, (1, 3, 6, 218): 13060,
     (1, 3, 6, 11, 218): 11091}
    #params = {'m': 601374, 'theta': 81232, 'A': 1000000, 'B': 1000}
    params = {'m': 990002, 'theta': 87743, 'A': 1000000, 'B': 1000}
elif dataset == 'test':
    itemsets = {
     (1,): 15, (2,): 8, (3,): 14, (4,): 6, (5,): 5, (6,): 11, (7,): 6, (8,): 4, (9,): 4, (10,): 5,
     (1, 2): 6, (1, 3): 10, (1, 4): 4, (1, 5): 3, (1, 6): 7, (1, 7): 3, (1, 8): 2, (1, 9): 3, (1, 10): 4,
     (2, 3): 7, (2, 4): 4, (2, 5): 3, (2, 6): 6, (2, 7): 1, (2, 8): 1, (2, 9): 1, (2, 10): 2,
     (3, 4): 4, (3, 5): 3, (3, 6): 8, (3, 7): 4, (3, 8): 2, (3, 9): 2, (3, 10): 1,
     (4, 5): 3, (4, 6): 5, (4, 7): 2, (4, 8): 0, (4, 9): 1, (4, 10): 0,
     (5, 6): 4, (5, 7): 2, (5, 8): 0, (5, 9): 0, (5, 10): 1,
     (6, 7): 2, (6, 8): 2, (6, 9): 1, (6, 10): 2,
     (7, 8): 2, (7, 9): 1, (7, 10): 1,
     (8, 9): 0, (8, 10): 0,
     (9, 10): 3,
     (1, 2, 3): 5, (1, 2, 4): 3, (1, 2, 5): 2, (1, 2, 6): 4, (1, 2, 7): 0, (1, 2, 8): 0, (1, 2, 9): 0, (1, 2, 10): 2,
     (1, 3, 4): 4, (1, 3, 5): 2, (1, 3, 6): 4, (1, 3, 7): 3, (1, 3, 8): 1, (1, 3, 9): 1, (1, 3, 10): 2,
     (1, 4, 5): 1, (1, 4, 6): 3, (1, 4, 7): 0, (1, 4, 8): 0, (1, 4, 9): 1, (1, 4, 10): 0,
     (1, 5, 6): 2, (1, 5, 7): 0
            }
    params = {'m': 15, 'theta': 12, 'A': 1000000, 'B': 1000}
else:
    raise Exception("Unknown dataset {}".format(dataset))

def valid_subsets(s):
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in
            range(1, len(s)))

def q(sup_a, sup_ab, method):
    m = params['m']
    theta = params['theta']
    A = params['A']
    B = params['B']
    if method == 'aab':
        return sup_a * sup_ab / (2 * m + 1)
    elif method == 'amab':
        return sup_a * (m - sup_ab) / (2 * m - 1)
    elif method == 'aabtab':
        return sup_a * sup_ab * (theta - sup_ab) / (2 * m * theta + theta - 1)
    elif method == 'tt' or method == 'mt':
        da = sup_a - (m if method[0] == 'm' else theta)
        dab = sup_ab - theta
        p = da * da + dab * dab
        return (A - math.sqrt(p)) / math.sqrt(2)
    elif method == 'exp':
        e = math.exp(1)
        return sup_a  * math.exp(-abs(theta - sup_ab) / B) / (m * (e - 1) + e)
    elif method == 'xy':
        def g(x, m, M, v):
            return v * max(0, 1.0 - abs(x - v + 0.0) / (M - m))
        return 0.5 * (1./theta - 1./m) * g(sup_a, theta, m, m) * g(sup_ab, theta, m, theta)
    else:
        raise Exception("Unknown method {}".format(method))

def build_rules():
    d = {}
    for i in itemsets:
        if len(i) > 1:
            for a in valid_subsets(i):
                b = list(set(i) - set(a))
                b.sort()
                b = tuple(b)
                sup_a = itemsets[a]
                sup_ab = itemsets[i]
                d[(sup_a, sup_ab)] = q(sup_a, sup_ab, method)
    return d

def print_rules(d):
    last_x = None
    for (x, y) in sorted(d.keys()):
        if x != last_x:
            if last_x: print
            last_x = x
        print x, y, d[(x, y)]

def Pr(q, Zinv = 1):
    epsilon = 0.1
    return Zinv * bigfloat.exp(epsilon * q / 2)

def print_pr_table(d):
    l = list(d.viewitems())
    l.sort(key = lambda x: x[1], reverse=True)
    z = 0.0
    for ((a, ab), q) in l:
        z += Pr(q)
    print z
    z = 1/z
    print "sup(A) sup(AB) q(A->B) Pr(A->B)"
    for ((a, ab), q) in l:
        print a, ab, q, Pr(q, z)

d = build_rules()
#print_rules(d)
print_pr_table(d)
