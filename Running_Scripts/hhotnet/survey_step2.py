## python 3
import subprocess
import sys

def run(n):
    n = int(n)

    #sim_matrix = "biogrid_similarity_matrix.h5"
    network_list = ["biogrid","irefindex18","reactome21","STRING_900"]
    data_version_list = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]
    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    permu = [1,2,3,4,5,6,7,8,9,10]
    score = ["".join((data_version,"_q_beta_",str(b),"_",str(p))) for b in beta for p in permu for data_version in data_version_list] * len(network_list)
    networks = [net for net in network_list for i in range(1100)]

    ##101*4400
    network = networks[n]
    s = score[n]
    file_path = "".join(("./Permutations/",network,"_",s,"/"))
    sim_matrix = network + "_similarity_matrix.h5"
    whole_input_scores = ["".join((file_path,"scores_",str(i),".tsv")) for i in range(101)]
    whole_output_edge = ["".join((file_path,"hierarchy_edge_list_",str(i),"_fast.tsv")) for i in range(101)]
    whole_output_index = ["".join((file_path,"hierarchy_index_gene_",str(i),"_fast.tsv")) for i in range(101)]
    whole_file = list(zip(whole_input_scores,whole_output_edge,whole_output_index))

    #chunks = [i for i in range(0,len(whole_file),101)]
    #this_chunk = whole_file[chunks[n]:chunks[n]+101]
    print("From",whole_file[0],"to",whole_file[-1])
    for cc in whole_file:
        subprocess.call(["python","construct_hierarchy.py","-smf",sim_matrix,"-igf",network+"_index_gene.tsv","-gsf",cc[0],"-helf",cc[1],"-higf",cc[2],"-st","0.1"])

if __name__ == "__main__":
    run(sys.argv[1])
