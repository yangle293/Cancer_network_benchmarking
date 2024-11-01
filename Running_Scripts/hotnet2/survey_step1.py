from __future__ import print_function
import subprocess
import sys

def run(n):
    n = int(n)
    network = ['biogrid','irefindex18','reactome21','STRING_900']
    beta = [0.48,0.47,0.38,0.33]
    num_network_permute = 100;
    num_heat_permute = 1000;
    n_cores = 1;
    sim_matrix = ["".join(('./data/networks/',net,'/',net,'_ppr_',str(b),'.h5')) for (net,b) in zip(network,beta) for i in range(1100)]
    permute_network_path = ["".join(('./data/networks/',net,'/permuted/',net,'_ppr_',str(b),'_##NUM##.h5')) for (net,b) in zip(network,beta) for i in range(1100)]
    data_version_list = ["corum_"+str(i) for i in range(5)] + ["reactome_"+str(i) for i in range(5)]
    beta_a = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    permu = [1,2,3,4,5,6,7,8,9,10]
    score = ["".join((data_version,"_q_beta_",str(b),"_",str(p))) for b in beta_a for p in permu for data_version in data_version_list] * len(network) 
    heat_file = ["".join(('./data/heats/',s,'.json')) for s in score]
    print(sim_matrix[n],permute_network_path[n],heat_file[n])
    subprocess.call(["python","../HotNet2.py","-nf",sim_matrix[n],"-pnp",permute_network_path[n],"-hf",heat_file[n],"-np",str(num_network_permute),"-hp",str(num_heat_permute),"-o","./results/","-c",str(n_cores)])

if __name__ == "__main__":
    run(sys.argv[1])
