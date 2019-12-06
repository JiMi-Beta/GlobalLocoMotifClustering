import networkx as nx
from networkx.algorithms import community
from collections import defaultdict
import sys

def edgeConductance(G, partitions):
    sum_conductance = 0

    for partition in partitions:
        T = set(G) - partition
        num_cut_edges = nx.cut_size(G, partition, T, weight=None)
        vol_S = nx.volume(G, partition, weight=None)
        vol_T = nx.volume(G, T, weight=None)

        if vol_S != 0 and vol_T != 0:
            conductance = num_cut_edges / min(vol_S,vol_T)
        else:
            continue

        sum_conductance += conductance

    avg_conductance = sum_conductance / len(partitions)

    return avg_conductance

def modularityAndEdgeCond(directed, edgeListFile, clusterFile, resultFile, toSort):
    if directed == True:
        G = nx.read_edgelist(edgeListFile, create_using=nx.DiGraph())
    else:
        G = nx.read_edgelist(edgeListFile, create_using=nx.Graph())

    partitionsInDict = defaultdict(list)

    # each line reqresents one cluster
    for line in open(clusterFile):

        line = line.strip().split()
        #Safety measure
        line = [x for x in line if x != '']
        partitionsInDict[len(line)].append(line)


    partitions = list()
    covered= set()
    if toSort == True:
        iterator = sorted(partitionsInDict.items(), reverse=True)
        print('Sorting')
    else:
        iterator = partitionsInDict.items()
        print('Not sorting')
    # using dictionary is much faster than sorting list of lists
    for i, j in sorted(partitionsInDict.items(), reverse=True):

        # one size may refer to multiple clusters
        for partition in j:

            oneSet = set(partition)
            overlapped = covered.intersection(oneSet)

            # overlapping exsits
            if bool(overlapped):
                oneSet -= overlapped

            covered.update(oneSet)
            partitions.append(oneSet)



    modularity = community.modularity(G, partitions)
    edge_cond = edgeConductance(G, partitions)

    print("modularity: %.3f"%modularity)
    print("average edge conductance: %.3f"%edge_cond)

    with open(resultFile,'w+') as f:
        f.write('Mod: %.5f \n'%modularity)
        f.write('EC: %.5f \n'%edge_cond)

if __name__ == '__main__':
    #Arguments: Directed/Undirected, edgelist, clusterfile, resultsfile
    assert len(sys.argv) > 5, 'Insert four command line arguments.\n 1:Directed/Undirected (U/D) \n 2: edgeListFile\n 3: outputClustersFile\n 4: resultsFile'
    directed = True if sys.argv[1] == 'D' else False
    edgeListFile = sys.argv[2]
    clusterFile = sys.argv[3]
    resultFile = sys.argv[4]
    toSort = True if sys.argv[5] == 'Y' else False

    if directed == True:
        print('Directed')
    else:
        print('Undirected')

    modularityAndEdgeCond(directed, edgeListFile, clusterFile, resultFile, toSort)

