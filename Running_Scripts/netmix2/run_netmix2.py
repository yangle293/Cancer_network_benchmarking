import argparse
import csv
import pandas as pd
import numpy as np
from scipy.stats import norm
import networkx as nx
import sys
import os

from locfdr import locfdr

from src.netmix2 import *
from src.common import *

# Parse arguments
def get_parser():
    description = ''
    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-el', '--edge_list', required=False, type=str, help='edge list file')
    parser.add_argument('-gs', '--gene_score', required=False, type=str, help='gene-to-score file')
    parser.add_argument('-d', '--delta', required=False, type=float, help='delta')
    parser.add_argument('-num_edges', '--num_edges', required=False, type=int, default=175000, help='Number of edges in G_delta for computing delta')
    parser.add_argument('-density', '--density', required=False, type=float, default=0.05, help='Target edge density (quadratic)')
    parser.add_argument('-time_limit', '--time_limit', required=False, type=int, default=12, help='Time limit (hours) for netmix function')
    parser.add_argument('-v', '--verbose', required=False, type=int, default=0, help='Print the progress of running NetMix2 (0, 1, or 2)')
    parser.add_argument('-o', '--output', required=False, type=str, help='directory for netmix results')
    return parser


def run(args):
    ###########################################################
    # load network
    node_list, A_network = load_network(args.edge_list, args.verbose)
    ###########################################################
    # load p-values for each gene
    pvals_list = load_pvalues(args.gene_score, node_list, args.verbose)
    ###########################################################
    # restrict to genes in the network
    (pvals_list, node_list, A_network) = restrict_to_genes_in_network(pvals_list, node_list, A_network, args.verbose)
    ###########################################################
    # compute zscores
    zscores = compute_zscores(pvals_list)
    zscores = post_process_zscores(zscores)
    ###########################################################
    # fit local fdr
    nulltype_ind = 1
    nulltype_name = "mlest"
    r_locfdr=locfdr(zscores, nulltype=1, plot=0)
    ###########################################################
    # correct the nans in the locfdr
    nonnull_count = correct_nans_from_locfdr(r_locfdr, zscores, nulltype_name)
    if args.verbose > 0:
        print(nonnull_count)
    ##################################################
    # Compute PPR matrix
    PPR_sim_mat, PPR_mat_rowsums = compute_ppr_kernel(A_network)
    ###########################################################
    # compute G_delta
    # threshold similarity matrix
    sim_mat = PPR_sim_mat
    
    sim_mat_nodiag=sim_mat-np.diag(np.diag(sim_mat))
    sim_mat_nodiag_sorted=np.sort(sim_mat_nodiag[np.triu_indices(sim_mat_nodiag.shape[0],1)])[::-1]
    
    if args.delta:
        delta=args.delta
    else:
        delta=sim_mat_nodiag_sorted[args.num_edges]

    sim_mat_delta=1*(sim_mat_nodiag > delta)
    degs_PPR_gwas = np.sum(sim_mat_delta,0)
    num_edges_delta = int(np.sum(sim_mat_delta)/2)
    
    if args.verbose > 0:
        print("delta: {}".format(delta))
        print('number of edges in G: {}'.format(np.sum(A_network)/2))
        print('number of edges in G_delta: {}'.format(np.sum(sim_mat_delta)/2))
    ###########################################################
    # parameters for netmix
    s=nonnull_count
    alpha=s/len(node_list)
    target_edge_density = args.density
    rho=target_edge_density*(s-1)
    
    if args.verbose > 0:
        print('size of altered subnetwork: {}'.format(s))
        print('target edge density: {}'.format(target_edge_density))
    
    ###########################################################
    # RUN NETMIX
    if args.verbose > 0:
        print('running netmix')
        print("time_limit", args.time_limit)

    est_subnetwork = netmix_edgedense(sim_mat_delta, rho, zscores, alpha, edge_dense_linear=True, time_limit=3600*args.time_limit)
    ###########################################################
    # print solution
    solution_size = len(est_subnetwork)
    est_subnetwork_genes = node_list[est_subnetwork] if solution_size>0 else []

    solution_network = sim_mat_delta[np.ix_(est_subnetwork, est_subnetwork)] if solution_size>0 else None
    num_edges_in_solution = sum(sum(solution_network))/2 if solution_size>0 else 0
    solution_network_density = num_edges_in_solution/(solution_size*(solution_size-1)/2) if solution_size>0 else 0

    print("Number of vertices in altered subnetwork: {}".format(solution_size))
    print("Density of altered subnetwork: {}".format(solution_network_density))
    ###########################################################
    # write solution
    if args.output:
        write_list_to_file(args.output+'_netmix_subnetwork.tsv', est_subnetwork_genes)
        write_list_to_file(args.output+'_node_list.tsv', node_list)
    
if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
