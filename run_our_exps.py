lens = [3, 4, 5]
ks = [10, 25, 50, 100]
epsilons = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
seeds = [123, 2345, 32525, 4235356, 5462356, 63562367, 73523, 835, 9372835, 10251235366]
datasets = [
    #"kosarak",
    #"msnbc 0.1 0 0.1 0.9",
    "mushroom 0.1 0.5 0.1 0.1",
    #"pumsb 0.1 0.5 0.1 0.1",
    #"retail"
    ]
dataset_path = 'graphs/datasets'
out_dir = 'out_www/dph2'
cmd = './dph'

def output(dataset, rlen, k, eps, seed):
    return '{}/{}_{}_{}_{}_{}'.format(out_dir, dataset,
            rlen, k, eps, seed)

def run(dataset, rlen, k, eps, seed):
    dataset, eps_share, c0, sm, cm = dataset.split()
    print '{} {}/{} r {} {} {} {} {} {} {} {} > {}'.format(cmd,
            dataset_path, dataset, eps, eps_share, k, rlen,
            c0, sm, cm, seed,
            output(dataset, rlen, k, eps, seed))

def main():
    for dataset in datasets:
        for l in lens:
            for k in ks:
                for eps in epsilons:
                    for seed in seeds:
                        run(dataset, l, k, eps, seed)

if __name__ == '__main__':
    main()
