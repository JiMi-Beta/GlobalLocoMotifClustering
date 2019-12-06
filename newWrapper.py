#1. ./convert -i graph.txt -o graph.bin
#2. ./community graph.bin -l -1 > graph.tree
#3. ./community graph.bin -l -1 -q 0.0001 > graph.tree
#4. ./hierarchy graph.tree (FIRST LINE OUTPUTS LEVELS OF COMMUNITIES/DEPTH 1-BASED)
#5. ./hierarchy graph.tree -l 2 > graph_node2comm_level2 (O-BASED INDEX OF DEPTH)

import os
import sys
from timeit import timeit
import subprocess
import numpy as np
from sklearn.metrics import jaccard_similarity_score
import pandas as pd
import getopt


'''
Compute Jaccard sim between two sets
'''
def jsim(s1,s2):
    #print(s1, s2)
    inters = len(s1.intersection(s2))
    union = len(s1) + len(s2) - inters
    return inters/union
'''
Compute precision, recall and F1
'''
def evaluate(truth, pred):
    TP = len(truth.intersection(pred))
    #Use identity len(pred) = TP + FP
    FP = len(pred) - TP
    #Use identity len(truth) = FN + TP
    FN = len(truth) - TP

    
    precision, recall = TP/(TP + FP), TP/(TP + FN)
    if precision == 0 and recall == 0:
        f1 = 0
    else:
        f1 = 2 * (precision * recall) / (precision + recall)
    #print(precision, recall, f1)
    return precision, recall, f1

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:t:n:r:", ["gfile=", "tfile=", "rfile=","name=", "runtime="])
    except getopt.GetoptError:
        print('script.py -g graphPath -t ground truth communities')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-g"):
            graph_filename = arg
        elif opt in ("-t"):
            clusterout_filename =  arg
        elif opt in ("-n"):
            intermed_name = arg
        elif opt in ("-r"):
            runtime_file = arg
    assert len(sys.argv) > 1, 'Please provide one command line argument,\
                                -g graphFile -t clusteroutfile -n name -r runtimefile'
    print('COUT FILENAME:',clusterout_filename, 'GRAPH FILENAME:', graph_filename)
    convert_cmd = './convert -i %s -o graph%s.bin'%(graph_filename,intermed_name)
    communities_cmd = './community graph%s.bin -l -1 > graph%s.tree'%(intermed_name,intermed_name)
    conv_time = timeit(stmt='os.system("%s")'%convert_cmd,setup='import os', number=1)
    #print('here')
    comm_time = timeit(stmt='os.system("%s")'%communities_cmd, setup='import os', number=1)
    print('Total time: %.3f'%(conv_time + comm_time))
    #sys.exit()
    with open(runtime_file,'w+') as f:
        f.write('Runtime: %.5f \n'%(conv_time + comm_time))

    comm_info = subprocess.check_output('./hierarchy graph%s.tree'%intermed_name,shell=True).decode('utf-8').split('\n')[0]
    #filename = 'hierarchy.txt'
    #with open('hierarchy.txt','r') as f:
    #    content = f.readline()

    depth = int(comm_info.split()[-1])

    depth_index = depth - 1

    write_to_fname = "graph_node2comm_level%s_%s"%(str(depth_index),intermed_name)
    os.system("./hierarchy graph%s.tree -l %d > %s"%(intermed_name, depth_index, write_to_fname))
    #print("./hierarchy graph.tree -l 2 > %s"%write_to_fname)

    #Read output file
    with open(write_to_fname,'r') as f:
        content = f.read()

    lines = np.matrix([y for y in [x.split() for x in content.split('\n')] if len(y)>0])
    #print(lines[:,1])

    unique_communities = np.unique(lines[:,1].tolist())

    df = pd.DataFrame(lines)

    #Output communities
    communities = dict()
    with open(clusterout_filename,'w+') as f:
        for comm in unique_communities:
            #print('Community %s'%comm)
            indices = df[1] == comm
            communities[comm] = list(set(df[0][indices]))
            to_write = ''
            for node in communities[comm]:
                to_write += str(node) + ' '
            to_write += '\n'
            f.write(to_write)

    #Read ground truth communities
    '''
    with open(gt_filename,'r') as f:
        gt_content = f.read()

    gt_communities = [set(y) for y in [x.split() for x in gt_content.split('\n')] if len(y) > 0]


    macroPrec, macroRec, macroF1 = 0.0, 0.0, 0.0
    #Evaluate each output cluster against ground truth community
    for outclust_index in communities.keys():
        outclust = communities[outclust_index]

        maxSim = -1
        curr_gt = ''

        #Find GT community with maximum overlap
        for gt_index, gt_com in enumerate(gt_communities):
            sim = jsim(outclust, gt_com)
            if sim > maxSim:
                maxSim = sim
                curr_gt = gt_com

        #Evaluate precision, recall, f1
        precision, recall, f1 = evaluate(curr_gt, outclust)
        macroPrec += precision
        macroRec += recall
        macroF1 += f1

    macroPrec /= len(communities.keys())
    macroRec /= len(communities.keys())
    macroF1 /= len(communities.keys())

    print('Macro-averaged precision %.3f'%macroPrec)
    print('Macro-averaged recall %.3f'%macroRec)
    print('Macro-averaged f1 score %.3f'%macroF1)

    with open(result_filename,'w+') as f:
        f.write("Precision %.5f\n"%macroPrec)
        f.write("Recall %.5f\n"%macroRec)
        f.write("F1 %.5f\n"%macroF1)
        f.write("Runtime %.5f\n"%(conv_time + comm_time))
    '''

