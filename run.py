datasets   = ['mushroom', 'pumsb_star', 'retail', 'kosarak']
thresholds = [        81,          490,      881,      9900]
k = 100

cmd_template = "./dph data/datasets/{0}.dat r {2} {1} 0.5 0.1 {3} 1000 5 {4}"
redirect = " > shelves/{0}_sh_{1}_bins_{2}"

for dataset, t in zip(datasets, thresholds):
    for sh in range(1,21):
        for bins in range(1,21):
            if sh * bins > k: continue
            print (cmd_template + redirect).format(dataset, bins,sh, t, k)
