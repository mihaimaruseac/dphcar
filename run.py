datasets   = ['mushroom', 'retail', 'pumsb_star']#, 'kosarak']
thresholds = [        81,      881,          490]#,      9900]
seeds = [149129418924, 2548298598, 3453959835, 4958958, 58929841412,
        621094029095, 750395903, 8593592905, 90358305305, 100035693539]
k = 100

#cmd_template = "./dph data/datasets/{0}.dat r {2} {1} 0.5 0.1 {3} {6} 3 {4} {5}"
cmd_template = "./dph data/datasets/{0}.dat w {2} {1} 0.5 0.1 {3} {6} 7 {4} {5}"
redirect = " > shelves_width_7/{0}_seed_{5}_sh_{1}_bins_{2}"

for dataset, t in zip(datasets, thresholds):
    for sh in range(1,11):
        for bins in range(1,11):
            if sh * bins > k: continue
            for seed in seeds:
                print (cmd_template + redirect).format(
                        dataset, bins,sh, t, k, seed, 2 * t)
            break
