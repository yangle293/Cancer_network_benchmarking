## python 3
import subprocess
import sys
import os

def run(n):
    n = int(n)
    data_path = "./data";
    result_folder = "./result"
    network_folder = "./network"
    network_list = ["biogrid","irefindex18","reactome21","string"]
    networks = network_list * 1100
    data_version_list = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]
    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    permu = [1,2,3,4,5,6,7,8,9,10]
    score = ["".join((data_version,"_diamond_",net,"_beta_",str(b),"_",str(p))) for b in beta for p in permu for data_version in data_version_list for net in network_list]
    #python3 DIAMOnD.py  network_file seed_file  n  alpha(optional)  outfile_name(optional)
    number_list = [100,200,300,400,500,600]
    network = networks[n]
    for number in number_list:
        subprocess.call(["python","DIAMOnD.py",os.path.join(network_folder,network+"_edge_list"),os.path.join(data_path,score[n]+".txt"),str(number),str(1),os.path.join(result_folder,score[n]+"_"+str(number)+"_module.txt")])

if __name__ == "__main__":
    run(sys.argv[1])
