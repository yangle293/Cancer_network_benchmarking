## python 3
import subprocess
import sys
import os
def run(n):
    n = int(n)
    data_path = "./data";
    network_path = "./network"
    result_folder = "./result"
    network_list = ["biogrid","irefindex18","reactome21","string"]
    network = [net for net in network_list for i in range(1100)]
    data_version_list = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]
    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    permu = [1,2,3,4,5,6,7,8,9,10]
    score = ["".join((data_version,"_list_beta_",str(b),"_",str(p))) for b in beta for p in permu for data_version in data_version_list] * len(network_list)

    s = score[n]
    net = network[n]
    subprocess.call(["python","netcore/netcore.py","-e",os.path.join(network_path,net+"_edge_list.tsv"),"-s",os.path.join(data_path,s+".txt"),"-pd",os.path.join(network_path,net+"_edge_permutations/"),"-o",os.path.join(result_folder,net+s+"_result.txt")])

if __name__ == "__main__":
    run(sys.argv[1])
