from __future__ import print_function
import networkx as nx
from SolveConductance import solveConductance
from SolveILP_Max import solveNetworkProblem
#import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import operator
import heapq
import time
from APPR_fast import approximatePPR_size_fast

def FDRnet(attGraph, bound, alpha, size, time_limit, relative_gap):
    ## seed selection
    scores = nx.get_node_attributes(attGraph, 'scores')
    values = [(x,scores[x]) for x in attGraph.nodes()]
    sorted_scores = sorted(values, key=operator.itemgetter(1))
    # select seed with local FDRs less than bound
    seed_list = [x for (x,y) in sorted_scores if y <= bound]
    ## Start searching subnetworks
    #print("Number of seeds: ", len(seed_list))
    total_start = time.time()
    result = list()
    for seed in seed_list:
        #if seed in seed_used:
        print("Seed ", seed, str(seed_list.index(seed)),"/",str(len(seed_list)))
        seed_start = time.time()
        result_seed, status = denseSubgraph(attGraph, bound, alpha, size, seed, time_limit, relative_gap)
        seed_end = time.time()
        result.append((seed, seed_end-seed_start, status, result_seed))
    total_end = time.time()
    print("Total time: ", total_end - total_start)
    return result


def denseSubgraph(attGraph, bound, alpha, size, seed, time_limit, relative_gap):
    ## TO DO: check parameters
    #
    #print("start to approximate PPR...")
    #ppr_start = time.time()
    ppr, converge, res = approximatePPR_size_fast(graph=attGraph, alpha=alpha, seed=seed, size=size)
    #ppr_end = time.time()
    #print("finished PPR computation...", ppr_end-ppr_start)
    normalized_ppr = dict((x, y/float(attGraph.degree[x])) for (x, y) in ppr.items())
    nx.set_node_attributes(attGraph, normalized_ppr, 'ppr')
    nodes = heapq.nlargest(int(size), normalized_ppr, key=normalized_ppr.get)
    G = attGraph.subgraph(nodes)
    #presolve: large problem may not get solution with size > 1 within 60s
    #print("pre-solve the problem...")
    #pre_start = time.time()
    subnetwork, status = solveNetworkProblem(G, bound, seed)
    #pre_end = time.time()
    #print("pre_solve time: ",pre_end-pre_start)
    if len(subnetwork) == 1 and status == 'MIP_optimal':
        #print("pre-solve find single gene as maximum solution...")
        return subnetwork, "".join(("presolve_",status))
    else:
        #print("pre-solve can not decide the result...")
        subnetwork, status = solveConductance(G, bound, seed, attGraph, time_limit, relative_gap)
        return subnetwork, status


#    pos = nx.spring_layout(G)
#    nx.draw_networkx(G,pos=pos)
#    nx.draw_networkx_nodes(G,pos=pos,nodelist=subnetwork,node_color='b')
#    tmp_G = G.subgraph(subnetwork)
#    nx.draw(tmp_G)
#    plt.show()
#    tmp_G = G.subgraph(subnetwork)
#    print("Result is connected:", nx.is_connected(tmp_G))



def load_network_from_file(index_file, edge_file, score_file):
    # Load gene-index map
    with open(index_file) as infile:
        arrs = [l.rstrip().split() for l in infile]
        indexToGene = dict((int(arr[0]), arr[1]) for arr in arrs)

    G = nx.Graph()
    G.add_nodes_from(indexToGene.values())  # in case any nodes have degree zero

    # Load graph
    with open(edge_file) as infile:
        edges = [map(int, l.rstrip().split()[:2]) for l in infile]
    G.add_edges_from([(indexToGene[u], indexToGene[v]) for u, v in edges])

    # Load p-values or scores
    with open(score_file) as infile:
        arrs = [l.rstrip().split() for l in infile]
        geneToScores = dict((arr[0], float(arr[1])) for arr in arrs)
    # TO DO: compute local fdr from p-value
    # merge scores and set score of no-score gene as 1
    for gene in indexToGene.values():
        if gene not in geneToScores.keys():
            geneToScores[gene] = 1.0
    nx.set_node_attributes(G, geneToScores, 'scores')
    return G


# Remove self-loops, multi-edges, and restrict to the largest component
def largest_component(G):
    selfLoops = [(u, v) for u, v in G.edges() if u == v]
    G.remove_edges_from(selfLoops)
    return G.subgraph(sorted(nx.connected_components(G), key=lambda cc: len(cc), reverse=True)[0])


def run(para):
    networkName = para[0]; cancerName=para[1]; bound = float(para[2]); alpha = float(para[3]); size = int(para[4])
    time_limit = para[5]; relative_gap = para[6]
    index_file = "".join(("./network/",networkName, "_index_gene"))
    edge_file = "".join(("./network/",networkName, "_edge_list"))
    score_file = "".join(("./data/",cancerName, ".txt"))
    G_raw = load_network_from_file(index_file=index_file, edge_file=edge_file, score_file=score_file)
    G = largest_component(G_raw)
#    print(G.node['ESR2']['scores'])
#    print(G.node['LRRC59']['scores'])
#    quit()
    result = FDRnet(G, bound=bound, alpha=alpha, size=size,time_limit = time_limit, relative_gap = relative_gap)
    result_file = "".join((networkName, cancerName, "_bound", str(bound),"_alpha",str(alpha), "_size", str(size), "_time",str(time_limit),"_gap",str(relative_gap),"_conductance.txt"))
    this_folder = os.path.dirname(os.path.abspath(__file__))
    my_file = os.path.join(this_folder,"result",result_file)
    with open(my_file, 'w') as outfile:
        for (seed,time,status, r) in result:
            outfile.write("".join((seed, "\t", str(time),"\t",status,"\t"," ".join(r), "\n")))


if __name__ == "__main__":
    run(sys.argv[1:])

