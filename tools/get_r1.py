import sys

s = 0
with open(sys.argv[1], "r") as f:
    line = f.readline()
    while line:
        line = line[:-1].split()[:-1]
        s = s + pow(2, len(line)) - 2
        line = f.readline()
print s
