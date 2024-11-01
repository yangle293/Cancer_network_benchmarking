from __future__ import print_function
import itertools
import sys
import numpy as np
from FDRnet_conductance import *

def run_final(n):
    network = ["biogrid","irefindex18","reactome21","string"]
    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    data_version = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]
    cancer = ["{}_lfdr_beta_{}_{}".format(d,j,i) for i in range(1,11) for j in beta for d in data_version]
    bound = [0.1]#[0.1,0.15,0.2,0.25]
    alpha = [0.85]
    #size = list(range(500,1010,10))
    time_limit = [100.0]
    relative_gap = [0.01]
    size = [400]
    allCombine = list(itertools.product(network, cancer,bound,alpha,size,time_limit,relative_gap))
    print(allCombine[int(n)])
    run(allCombine[int(n)])


if __name__ == "__main__":
    run_final(sys.argv[1])
