import os
import subprocess
import shutil


def run():
    data_path = "./Data"
    intermediate_path ="./Permutations"
    results_path = "./Result"
    network_list = ["STRING_900"]
    n_permutations = 100
    #sim_matrix = network+"_similarity_matrix.h5"

    beta = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11]
    data_version = ["corum_"+str(k) for k in range(5)] + ["reactome_"+str(k) for k in range(5)]#["sampled_complex_v1","sampled_complex_v2","sampled_path_v1","sampled_path_v2"]
    score = ["{}_q_beta_{}_{}".format(d,j,i) for i in range(1,11) for j in beta for d in data_version]
    print(score)
    ## rename the score file to .tsv files
    print("renaming score files")
    for network in network_list:
        if os.path.exists(os.path.join(data_path,network+"_edge_list")):
            os.rename(os.path.join(data_path,network+"_edge_list"),os.path.join(data_path,network+"_edge_list.tsv"))
        if os.path.exists(os.path.join(data_path,network+"_index_gene")):
            os.rename(os.path.join(data_path,network+"_index_gene"),os.path.join(data_path,network+"_index_gene.tsv"))
    for s in score:
        if os.path.exists(os.path.join(data_path,s+".txt")):
            os.rename(os.path.join(data_path,s+".txt"),os.path.join(data_path,s+".tsv"))

    ## make network directory
    print("making directories")
    # if not os.path.exists(os.path.join(intermediate_path,network)):
    #     os.mkdir(os.path.join(intermediate_path,network))

    ## make score directory
    for network in network_list:
        for s in score:
            if not os.path.exists(os.path.join(intermediate_path,network+'_'+s)):
                os.mkdir(os.path.join(intermediate_path,network+'_'+s))

    ## construct similarity sim matrix
    for network in network_list:
        if not os.path.exists(os.path.join("./",network+"_similarity_matrix.h5")):
            print("constructing similarity matrix")
            subprocess.call(["python","construct_similarity_matrix.py","-i",os.path.join("./",network+"_edge_list.tsv"),"-o",os.path.join("./",network+"_similarity_matrix.h5"),"-bof",os.path.join("./",network+"beta.txt")])

    ## Permute the data
    print("Permuting the data")
    for network in network_list:
        for s in score:
            # move s.tsv to intermediate_path and name it with scores_0.tsv
            shutil.copyfile(os.path.join(data_path,s+".tsv"), os.path.join(intermediate_path,network+'_'+s,"scores_0.tsv"))
            # run find_permutation_bins
            subprocess.call(["python","find_permutation_bins.py","-gsf",os.path.join(intermediate_path,network+'_'+s,"scores_0.tsv"),"-igf",os.path.join("./",network+'_index_gene.tsv'),"-elf",os.path.join("./",network+'_edge_list.tsv'),"-ms",str(1000),"-o",os.path.join(intermediate_path,network+'_'+s,"score_bins.tsv")])
            # for each permutation, run permute_scores
            for n in range(1,n_permutations+1):
                subprocess.call(["python","permute_scores.py","-i",os.path.join(intermediate_path,network+'_'+s,"scores_0.tsv"),"-bf",os.path.join(intermediate_path,network+'_'+s,"score_bins.tsv"),"-s",str(n),"-o",os.path.join(intermediate_path,network+'_'+s,"scores_"+str(n)+".tsv")])




if __name__ == "__main__":
    run()
