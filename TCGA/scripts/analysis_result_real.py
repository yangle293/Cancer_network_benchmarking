import os
import glob
from tqdm.auto import tqdm
import pickle
import re
import networkx as nx
import csv
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import pandas as pd

# helper functions
def save_results(results, file_name):
    with open(file_name, "wb") as f:
        pickle.dump(results, f)

# def flatten_results(results):
#     flattened = []
#     for method, method_results in results.items():
#         for (data_version, network_name), subnetworks in method_results.items():
#             if isinstance(subnetworks, dict):
#                 for params, subnetwork_list in subnetworks.items():
#                     flattened.append((method, data_version, network_name,params, subnetwork_list))
#             else:
#                 flattened.append((method, data_version, network_name, None, subnetworks))
#     return flattened

# def calculate_f_score(genes1,genes2):
#     f_score = 0
#     set1 = set(genes1)
#     set2 = set(genes2)
#     # Calculate true positives (TP)
#     true_positives = len(set1.intersection(set2))
#     # Calculate false positives (FP) and false negatives (FN)
#     false_positives = len(set2.difference(set1))
#     false_negatives = len(set1.difference(set2))
#     # Calculate precision and recall
#     if true_positives + false_positives == 0:
#         precision = 1
#     else:
#         precision = true_positives / (true_positives + false_positives)

#     if true_positives + false_negatives == 0:
#         recall = 1
#     else:
#         recall = true_positives / (true_positives + false_negatives)

#     # Calculate F-score
#     if precision + recall == 0:
#         f_score = 0
#     else:
#         f_score = 2 * (precision * recall) / (precision + recall)

#     return f_score, precision, recall

# def evaluate_performance(subnetworks, cosmic):
#     fscore = 0
#     fsub_symmetric = 0
#     total_genes_pred = set().union(*subnetworks)
#     return calculate_f_score(total_genes_pred,cosmic)

# def calculate_performance(flattened_results,cosmic):
#     performance_results = []
#     for method, data_version, network_name,params, subnetworks in flattened_results:
#         subs = [s for s in subnetworks if len(s) > 1]
#         try:
#             fscore, precision,recall = evaluate_performance(subs, cosmic)
#             performance_results.append((method, data_version, network_name,params, subs, fscore, precision,recall))
#         except Exception as e:
#             print(e)
#             print(method, data_version, network_name,params, subnetworks)
#             quit()
#     # Organize the performance results into a pandas DataFrame:
#     df_raw = pd.DataFrame(performance_results, columns=["method", "data_version", "network_name", "params","subnetworks" ,"fscore","precision","recall"])
#     # reduce the params by only select max fscore
#     idx = df_raw.groupby(['method', 'data_version', 'network_name'])['fscore'].idxmax()
#     df_tmp = df_raw.loc[idx].drop(columns="params").copy()

#     performance = list(df_tmp.itertuples(index=False, name=None))
#     return performance

def calculate_lfdr(results,ground_truths):
    results_with_lfdr = []
    for method, method_results in results.items():
        for (data_version, network_name), subs in method_results.items():
            score = ground_truths['scores'][data_version]
            lfdr = []
            for sub in subs:
                tmp_lfdr = [score.loc[score['Gene'] == s, 'fdr'].iloc[0] if (score['Gene'].eq(s)).any() else 1. for s in sub]
                lfdr.append(np.mean(tmp_lfdr))
            results_with_lfdr.append((method, data_version, network_name, subs,lfdr))
    df = pd.DataFrame(results_with_lfdr, columns=["method", "data_version", "network_name", "subnetworks" ,"lfdr"])
    return df

def main():
    # load results and ground_truth
    print("loading all results and ground truth")
    with open("all_results_real.pkl", "rb") as f:
        results = pickle.load(f)
    with open("ground_truth_real.pkl", "rb") as f:
        ground_truths = pickle.load(f)

    # if os.path.isfile('lfdr_real.pkl'):
    #     print("loading lfdr...")
    #     with open('lfdr_real.pkl', 'rb') as f:
    #         df_lfdr = pickle.load(f)
    # else:
    print("Calculating lfdr...")
    df_lfdr = calculate_lfdr(results,ground_truths)
    save_results(df_lfdr,'lfdr_real.pkl')
    # # load cosmic
    # cosmic_all = pd.read_csv("../code/score/real/Census_allMon Oct 21 15_54_16 2024.csv",sep=',')
    # cosmic = list(cosmic_all['Gene Symbol'])

    # # flatten the results
    # print("flattening results...")
    # if os.path.isfile('flatten_real.pkl'):
    #     print("loading flattening results...")
    #     with open('flatten_real.pkl', 'rb') as f:
    #         flattened_results = pickle.load(f)
    # else:
    #     print("calculating flattening results...")
    #     flattened_results = flatten_results(results)
    #     save_results(flattened_results,'flatten_real.pkl')


    # # 1. evaluate general performance
    # if os.path.isfile('performance_real.pkl'):
    #     print("loading f and fsub...")
    #     with open('performance_real.pkl', 'rb') as f:
    #         performance = pickle.load(f)
    # else:
    #     print("Calculating f and fsub...")
    #     performance = calculate_performance(flattened_results,cosmic)
    #     save_results(performance,'performance_real.pkl')

    #df_raw = pd.DataFrame(performance, columns=["method", "data_version", "network_name", "subnetworks" ,"fscore","precision","recall"])



    # performance_with_lfdr = list(df_lfdr.itertuples(index=False, name=None))


    # if os.path.isfile('rlfdr_cancer.pkl'):
    #     print("loading rlfdr...")
    #     with open('rlfdr_cancer.pkl', 'rb') as f:
    #         df_rlfdr = pickle.load(f)
    # else:
    #     print("Calculating rlfdr...")
    #     df_rlfdr = calculate_real_lfdr(performance_with_lfdr,ground_truths)
    #     save_results(df_rlfdr,'rlfdr_cancer.pkl')
    # # 2. The relationship between node-level measure and detection rate
    # if os.path.isfile('node_detect_cancer.pkl'):
    #     print("loading node detection...")
    #     # with open('node_detect.pkl', 'rb') as f:
    #     #     node_detect_dict = pickle.load(f)
    # else:
    #     print("Calculating node detection...")
    #     node_detect_dict = node_detection(df_lfdr,performance,ground_truths)
    #     save_results(node_detect_dict,'node_detect_cancer.pkl')

    # # 3. The relationship between network-level measure and detection rate
    # if os.path.isfile('network_detect_cancer.pkl'):
    #     print("loading network detection...")
    #     # with open('network_detect.pkl', 'rb') as f:
    #     #     network_detect_dict = pickle.load(f)
    # else:
    #     print("Calculating network detection...")
    #     network_detect_dict = subnetwork_detection(df_lfdr,performance,ground_truths)
    #     save_results(network_detect_dict,'network_detect_cancer.pkl')

    # # # 4. Analyze the overlap between the subnetworks detected by different methods
    # # if os.path.isfile('overlap.pkl'):
    # #     print("loading overlap...")
    # #     # with open('overlap.pkl', 'rb') as f:
    # #     #     overlap_result = pickle.load(f)
    # # else:
    # #     print("Calculating overlap...")
    # #     overlap_result = network_overlap(performance,ground_truths)
    # #     save_results(overlap_result,'overlap.pkl')


if __name__ == "__main__":
    main()
