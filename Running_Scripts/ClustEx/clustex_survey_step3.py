import subprocess
import sys
import os

def run(n):
    n = int(n)
    data_path = "./data";
    result_folder = "./result"
    network_folder = "./network"
    network_list = ["biogrid","irefindex18","reactome21","string"]
    data_version_list = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]
    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    permu = [1,2,3,4,5,6,7,8,9,10]
    score = ["".join((data_version,"_clustex_",net,"_beta_",str(b),"_",str(p))) for b in beta for p in permu for data_version in data_version_list for net in network_list]
    networks = network_list*1100
    print(networks[n],score[n])
    #./clustex2 --gene_list irefindex9_high_%%%.txt --network irefindex9_network.txt --diffusion_kernel irefindex9_diffusion_kernel_0.01.txt --gene_importance irefindex9_high%%%_gene_importance_0.5.txt -C -n 0.03 -s 50 -j 2irefindex9_high%%%
    subprocess.call(["./clustex2","--gene_list",os.path.join(data_path,score[n]+".txt"),"--network",os.path.join(network_folder,networks[n]+"_edge_list_clustex.txt"),"--diffusion_kernel",networks[n]+"_diffusion_kernel_0.01.txt","--gene_importance",score[n]+"_gene_importance_0.5.txt","-C","-n",str(0.03),"-s",str(50),"-j",score[n]+str(2)])

if __name__ == "__main__":
    run(sys.argv[1])
