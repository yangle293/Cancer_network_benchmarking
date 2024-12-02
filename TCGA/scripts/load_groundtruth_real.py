import os
import glob
from tqdm.auto import tqdm
import pickle
import re
import networkx as nx
import csv
from collections import defaultdict
import pandas as pd
from functools import reduce
from load_results import save_results

def load_graph_structure(edge_list_file):
    with open(edge_list_file, "r") as f:
        lines = f.readlines()[1:] # Skip the first row
    graph = nx.parse_edgelist(lines, delimiter='\t', nodetype=str, data=(('weight', float),))
    return graph

def load_groundtruth(base_path):
    data_versions = ["BLCA_MERGE_v4","LUAD_MERGE_v4","LUSC_MERGE_v4","COADREAD_MERGE_v4","PRAD_MERGE_v4","HNSC_MERGE_v4",
                 "UCEC_MERGE_v4","KIPAN_MERGE_v4","BRCA_MERGE_v4"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))
    score_types = ["p","q","z","fdr"]
    network_files = {
        'biogrid': '/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/biogrid_edge_list.tsv',
        'irefindex18': '/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/irefindex18_edge_list.tsv',
        'reactome21': '/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/reactome21_edge_list.tsv',
        'string': '/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/string_edge_list.tsv',
    }


    networks = {}
    targets = {}
    scores = {}

    for key in network_files:
        n = load_graph_structure(network_files[key])
        networks[key] = n


    for data_version in data_versions:
        df_list = []
        for score_type in score_types:
            file_name = f"{data_version}_{score_type}.txt"
            file_path = os.path.join(base_path,"real",score_type,file_name)
            df = pd.read_csv(file_path,sep='\t',header=None,names=['Gene',score_type])
            df_list.append(df)

        score_df = reduce(lambda x, y: pd.merge(x, y, on = 'Gene'), df_list)
        key = (data_version)
        scores[key] = score_df

    ground_truth = {}
    ground_truth['networks'] = networks
    ground_truth['scores'] = scores
    return ground_truth



if __name__ == "__main__":
    ground_truth = load_groundtruth('/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/score')
    save_results(ground_truth,"ground_truth_real.pkl")
