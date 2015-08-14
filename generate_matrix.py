import sys

data = {}

for fname in sys.argv[1:]:
    parts = fname.split('_')
    bins = int(parts[-1])
    shelves = int(parts[-3])

    with open(fname, "r") as f:
        recordMode = False
        while True:
            line = f.readline()
            if not line: break # eof
            line = line[:-1]
            if not line: continue # empty line
            if recordMode:
                parts = line.split()
                c = float(parts[0])
                r = int(parts[-2])
                data_shelves = data.get(c, {})
                data_bins = data_shelves.get(shelves, {})
                data_bins[bins] = data_bins.get(bins, 0) + r
                data_shelves[shelves] = data_bins
                data[c] = data_shelves
            if line == "Final histogram:":
                recordMode = True

for c in [0.9, 0.8, 0.7, 0.6, 0.5]:
    print c
    print ', ', ', '.join(map(str, range(1, 21)))
    data_shelves = data.get(c, {})
    for shelf in range(1, 21):
        print shelf,
        data_bins = data_shelves.get(shelf, {})
        for bin in range(1, 21):
            print ', ', data_bins.get(bin, '-'),
        print
    print
    print
