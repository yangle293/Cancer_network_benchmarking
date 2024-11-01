import subprocess
import sys
import os
# step 1: calculate diffusion matrix once
def run():
    data_path = "./data";
    result_folder = "./result"
    network_folder = "./network"
    network_list = ["biogrid","irefindex18","reactome21","string"]
    data_version_list = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]
    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    permu = [1,2,3,4,5,6,7,8,9,10]
    #score = ["".join((data_version,"_clustex_beta_",str(b),"_",str(p))) for b in beta for p in permu for data_version in data_version_list] * len(network_list)
    #./clustex2 --gene_list biogrid_%%%.txt --network biogrid_network.txt --diffusion_kernel biogrid_diffusion_kernel_0.01.txt -G -P 3 -j biogrid%%%
    #print(score[0])
    #print(os.path.join(data_path,score[0]))
    for network in network_list:
        score = "".join("corum_0_clustex_"+network+"_beta_0.01_1.txt")
        subprocess.call(["./clustex2","--gene_list",os.path.join(data_path,score),"--network",os.path.join(network_folder,network+"_edge_list_clustex.txt"),"-D","-G","-P",str(3),"-j",network])

if __name__ == "__main__":
    run()
