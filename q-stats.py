import itertools
import math

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

# DO NOT USE THESE
#method = 'aab'    # sup(A) * sup(AB)
#method = 'amab'   # sup(A) * (m - sup(AB))
#method = 'aabtab' # sup(A) * sup(AB) * (theta - sup(AB))
#method = 'exp'    # sup(A) * exp(- |theta - sup(AB)|)

#method = 'tt'     # distance to (theta, theta)
#method = 'mt'     # distance to (m, theta)
method = 'xy'     # linear biproduct of two triangles

def valid_subsets(s):
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in
            range(1, len(s)))

def q(sup_a, sup_ab, method):
    m = 601374
    theta = 81232
    A = 1000000
    B = 1000
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
        raise Exception("Uknown method {}".format(method))

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

d = build_rules()
print_rules(d)
