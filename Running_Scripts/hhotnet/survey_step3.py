import subprocess
import sys
import os

def run(n):
    data_path = "./Data"
    intermediate_path ="./Permutations"
    results_path = "./Result"
    network_list = ["biogrid","irefindex18","reactome21","STRING_900"]
    n_permutations = 100

    network = network_list[int(n)]
    sim_matrix = network+"_similarity_matrix.h5"
    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    data_version = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]
    score = ["{}_q_beta_{}_{}".format(d,j,i) for i in range(1,11) for j in beta for d in data_version]

    for s in score:
        file_path = os.path.join("./Permutations",network+"_"+s)
        observed_edge_list = os.path.join(file_path,"hierarchy_edge_list_0_fast.tsv")
        observed_index_gene = os.path.join(file_path,"hierarchy_index_gene_0_fast.tsv")
        permutation_edge_list = [os.path.join(file_path,"hierarchy_edge_list_"+str(i)+"_fast.tsv") for i in range(1,n_permutations+1)]
        permutation_index_gene = [os.path.join(file_path,"hierarchy_index_gene_"+str(i)+"_fast.tsv") for i in range(1,n_permutations+1)]
        results = os.path.join(results_path,"clusters_"+network+"_"+s+"_fast.tsv")
        pdfs = os.path.join(results_path,"sizes_"+network+"_"+s+"_fast.pdf")
        subprocess.call(["python","process_hierarchies.py","-oelf",observed_edge_list,"-oigf",observed_index_gene,"-pelf",*permutation_edge_list,"-pigf",*permutation_index_gene,"-lsb","10","-usb","100","-cf",results,"-pl",network+" "+s,"-pf",pdfs])

if __name__ == "__main__":
    run(sys.argv[1])
