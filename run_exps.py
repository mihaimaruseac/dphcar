lens = [5] #3, 5, 7
ks = [10]
epsilons = [0.1, 0.25, 0.5, 0.75, 1.0]
seeds = [123, 2345, 32525, 4235356, 5462356, 63562367, 73523, 835, 9372835, 10251235366]
datasets = [
    #"kosarak",
    "msnbc 17 14795",
    "mushroom 119 23",
    "pumsb 255 22",
    #"retail"
    ]
dataset_path = 'graphs/datasets'

out_dir = 'out_www/ngram'

ngram_cmd = 'python ngram_options.py'

def output(dataset, rlen, k, eps, seed):
    return '{}/{}_{}_{}_{}_{}'.format(out_dir, dataset.split()[0],
            rlen, k, eps, seed)

def run(dataset, rlen, k, eps, seed):
    print '{} {}/{} {} {} {} {} > {}'.format(ngram_cmd,
            dataset_path, dataset,
            rlen, k, eps, seed,
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
