## python 3
import subprocess
import sys
import os

def run(n):
    n = int(n)
    data_path = "./data"
    network_folder = "./network"
    result_folder = "./result"
    network_list = ["biogrid","irefindex18","reactome21","string"]
    networks = [net for net in network_list for i in range(1100)]
    data_version_list = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]
    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    permu = [1,2,3,4,5,6,7,8,9,10]
    score = ["".join((data_version,"_list_beta_",str(b),"_",str(p))) for b in beta for p in permu for data_version in data_version_list] * len(network_list)
    network = networks[n]
    s = score[n]
    #domino --active_genes_files </path/to/dataset1,/path/to/dataset2...> --network_file </path/to/network.sif> --slices_file <slices_file.txt> --output_folder </path/to/output_folder> [-sth <slices_threshold> -mth <putative_modules_threshold>]
    subprocess.call(["domino","--active_genes_files",os.path.join(data_path,score[n]+".txt"),"--network_file",os.path.join(network_folder,network+"_edge_list.sif"),"--slices_file",os.path.join(network_folder,network+"_slices.txt"),"--output_fold",result_folder+"/"+network])

if __name__ == "__main__":
    run(sys.argv[1])
