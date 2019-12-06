import sys

if __name__ == '__main__':
    assert len(sys.argv) > 2, 'Enter command line argument: filename, resultFname'
    filename = sys.argv[1]
    resname = sys.argv[2]

with open(filename, 'r') as f:
    content = f.readlines()

content = [x.strip().split() for x in content]
content = [x for x in content if not x[0][0] == '#']
content = [(int(x[0]), int(x[1])) for x in content if len(x) == 2]
#print(content[:10])

mapping = dict()  # Map from text nodeIDs to consecutive numbers starting at 0
counter = 0

with open(resname, 'w+') as f:
    for u, v in content:
        if not u in mapping:
            mapping[u] = counter
            counter += 1
        if not v in mapping:
            mapping[v] = counter
            counter += 1
        f.write(str(mapping[u]) + ' ' + str(mapping[v]) + '\n')
