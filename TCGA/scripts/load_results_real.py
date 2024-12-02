import os
import glob
from tqdm.auto import tqdm
import pickle
import re
import networkx as nx
import csv
from collections import defaultdict

def save_results(results, file_name):
    with open(file_name, "wb") as f:
        pickle.dump(results, f)


def load_graph_structure(edge_list_file):
    with open(edge_list_file, "r") as f:
        lines = f.readlines()[1:] # Skip the first row
    graph = nx.parse_edgelist(lines, delimiter='\t', nodetype=str, data=(('weight', float),))
    return graph

def calculate_seed_connections(results, graph):
    seed_connections = {}
    for seed, subnetwork_list in results.items():
        connection_count = 0
        for other_seed in results.keys():
            if seed != other_seed and graph.has_edge(seed, other_seed):
                connection_count += 1
        seed_connections[seed] = connection_count
    return seed_connections

def filter_subnetworks(results, seed_connections):
    sorted_seeds = sorted(seed_connections, key=seed_connections.get, reverse=True)
    filtered_results = []
    seen_seeds = set()

    for seed in sorted_seeds:
        subnetwork_list = results[seed]
        for subnetwork in subnetwork_list:
            if seed not in seen_seeds and len(subnetwork) > 1:
                filtered_results.append(subnetwork)
                seen_seeds.update(subnetwork)

    return filtered_results



def load_gene_name_mapping(index_file, name_file):
    gene_name_mapping = {}

    with open(index_file, "r") as f_index, open(name_file, "r") as f_name:
        for index_line, name_line in zip(f_index, f_name):
            index_gene1, index_gene2 = index_line.strip().split()
            name_gene1, name_gene2 = name_line.strip().split()

            gene_name_mapping[index_gene1] = name_gene1
            gene_name_mapping[index_gene2] = name_gene2

    return gene_name_mapping

def load_gene_name_mapping_clustex():
    network_files = {
    "biogrid": {
        "index": "/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/biogrid_edge_list_clustex.txt",
        "name": "/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/biogrid_edge_list.tsv",
    },
    "irefindex18": {
        "index": "/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/irefindex18_edge_list_clustex.txt",
        "name": "/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/irefindex18_edge_list.tsv",
    },
    "reactome21": {
        "index": "/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/reactome21_edge_list_clustex.txt",
        "name": "/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/reactome21_edge_list.tsv",
    },
    "string": {
        "index": "/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/string_edge_list_clustex.txt",
        "name": "/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/string_edge_list.tsv",
    },
    }

    # Load gene name mappings for each network
    network_gene_name_mappings = {
        network_name: load_gene_name_mapping(files["index"], files["name"])
        for network_name, files in network_files.items()
        }
    return network_gene_name_mappings




def load_robust_results(base_path,data_versions):
    results = {}
    #data_versions = ["corum_cancer"] + ["reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))

    progress_bar = tqdm(total=len(data_versions)*len(network_names), desc="Processing robust results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            #for beta in betas:
                #for permu in permus:
            progress_bar.update(1)
            file_name = f"{network_name}{data_version}_list_result.txt"
            file_path = os.path.join(base_path, file_name)

            if not os.path.exists(file_path):
                continue

            subnetwork = load_robust_subnetwork(file_path)
            key = (data_version, network_name)
            results[key] = subnetwork

    progress_bar.close()
    return results

def load_robust_subnetwork(file_path):
    subnetwork = set()

    with open(file_path, "r") as f:
        reader = csv.reader(f, delimiter=" ")
        for row in reader:
            subnetwork.add(row[0])
            subnetwork.add(row[1])

    return [subnetwork]


def load_regmod_results(base_path,data_versions):
    results = {}
    #data_versions = ["corum_cancer"] + ["reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))
    # b_values = [1,2,3,4,5]
    # c_values = [1,5,10]
    # f_values = [2,4,6]
    b = 3;c=5;f=4;
    progress_bar = tqdm(total=len(data_versions) * len(network_names), desc="Processing regmod results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            progress_bar.update(1)
            subnetworks = []
            # for b in b_values:
            #     for c in c_values:
            #         for f in f_values:

            file_name = f"{network_name}{data_version}_z_b{b}_c{c}_f{f}.txt"
            file_path = os.path.join(base_path, file_name)
            if os.path.exists(file_path):
                # subs = []
                with open(file_path, "r") as fp:
                    reader = csv.reader(fp, delimiter=" ")
                    for row in reader:
                        s = list(set(row))
                        subnetworks.append(s)
                    # subkey = (b,c,f)
                    # subnetworks[subkey] = subs

            key = (data_version, network_name)
            results[key] = subnetworks
    progress_bar.close()
    return results

def load_netmix2_results(base_path,data_versions):
    results = {}
    #data_versions = ["corum_cancer"] + ["reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))

    progress_bar = tqdm(total=len(data_versions) * len(network_names), desc="Processing NetMix2 results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            progress_bar.update(1)
            file_path = os.path.join(base_path, f"{network_name}{data_version}_p_netmix_subnetwork.tsv")

            if os.path.exists(file_path):
                subnetwork = load_netmix2_subnetwork(file_path)
                key = (data_version, network_name)
                results[key] = [subnetwork]

    progress_bar.close()
    return results

def load_netmix2_subnetwork(file_path):
    subnetwork = set()

    with open(file_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            subnetwork.add(row[0])

    return subnetwork


def load_netcore_results(base_path,data_versions):
    results = {}
    #data_versions = ["corum_cancer"] + ["reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))

    progress_bar = tqdm(total=len(data_versions) * len(network_names), desc="Processing Netcore results", position=0)


    for data_version in data_versions:
        for network_name in network_names:
            # for beta in betas:
            #     for permu in permus:
            progress_bar.update(1)
            file_path = os.path.join(base_path, f"{network_name}{data_version}_list_result.txtcore_norm_subnetworks.txt")

            if os.path.exists(file_path):
                subnetworks = load_netcore_subnetworks(file_path)
                key = (data_version, network_name)
                results[key] = subnetworks

    progress_bar.close()
    return results

def load_netcore_subnetworks(file_path):
    subnetworks = []

    with open(file_path, "r") as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            subnetwork = row[0].split("\t")
            subnetwork = set(subnetwork)
            subnetworks.append(subnetwork)

    return subnetworks


def load_must_results(base_path,data_versions):
    results = {}
    #data_versions = ["corum_cancer"] + ["reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))

    progress_bar = tqdm(total=len(data_versions) * len(network_names), desc="Processing MUST results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            # for beta in betas:
            #     for permu in permus:
            progress_bar.update(1)
            file_name = f"{network_name}{data_version}_list_result.txt"
            file_path = os.path.join(base_path, file_name)

            if not os.path.exists(file_path):
                continue

            subnetwork = load_must_subnetwork(file_path)
            key = (data_version, network_name)
            results[key] = subnetwork

    progress_bar.close()
    return results

def load_must_subnetwork(file_path):
    subnetwork = set()

    with open(file_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # Skip header
        for row in reader:
            subnetwork.add(row[0])
            subnetwork.add(row[1])

    return [subnetwork]



def load_hotnet2_results(base_path,data_versions):
    results = {}
    #data_versions = ["corum_cancer"] + [f"reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))

    progress_bar = tqdm(total=len(data_versions) * len(network_names), desc="Processing hotnet2 results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            progress_bar.update(1)
            dir_path = os.path.join(base_path, f"{network_name}-{data_version}_q")
            if network_name == "string":
                dir_path = dir_path.replace("string", "string_900")

            deltas = glob.glob(os.path.join(dir_path, "delta_*"))
            delta_values = [float(os.path.basename(delta).split("_")[1]) for delta in deltas]
            min_delta = float("inf")
            chosen_delta_subnetworks = None
            max_size = 0

            for delta, delta_path in zip(delta_values, deltas):
                subnetworks = load_hotnet_components(os.path.join(delta_path, "components.txt"))
                sig_subnetworks, min_sig_size = load_hotnet_significance(os.path.join(delta_path, "significance.txt"))
                filtered_subnetworks = [subnet for subnet in subnetworks if len(subnet) >= min_sig_size]

                if (len(filtered_subnetworks) > max_size) or (len(filtered_subnetworks) == max_size and delta < min_delta):
                    min_delta = delta
                    max_size = len(filtered_subnetworks)
                    chosen_delta_subnetworks = filtered_subnetworks

            key = (data_version, network_name)
            if chosen_delta_subnetworks:
                results[key] = chosen_delta_subnetworks
            else:
                continue

    progress_bar.close()
    return results

def load_hotnet_components(file_path):
    subnetworks = []

    with open(file_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            # Remove the newline character from the last gene and append the subnetwork
            row[-1] = row[-1].rstrip()
            subnetworks.append(set(row))

    return subnetworks



def load_hotnet_significance(file_path):
    min_sig_size = float("inf")
    sig_subnetworks = 0

    with open(file_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # Skip header
        for row in reader:
            size, expected, actual, p_value = int(row[0]), float(row[1]), int(row[2]), float(row[3])
            if p_value < 0.05:
                min_sig_size = min(min_sig_size, size)
                sig_subnetworks += 1

    return sig_subnetworks, min_sig_size



def load_hhotnet_results(base_path,data_versions):
    results = {}

    #data_versions = ["corum_cancer"] + ["reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))

    progress_bar = tqdm(total=len(data_versions) * len(network_names), desc="Processing hhotnet results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            # for beta in betas:
            #     for permu in permus:
            progress_bar.update(1)
            file_name = f"clusters_{network_name}_{data_version}_q_fast.tsv"
            if network_name == "string":
                file_name = file_name.replace("string", "STRING_900")

            file_path = os.path.join(base_path, file_name)

            if not os.path.exists(file_path):
                continue

            with open(file_path, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                subnetworks = []
                for row in reader:
                    if len(row) > 1:
                        subnetwork_genes = set(row)
                        subnetworks.append(subnetwork_genes)

            key = (data_version, network_name)
            results[key] = subnetworks

    progress_bar.close()
    return results


def load_fdrnet_results(base_path,data_versions):
    results = {}

    #data_versions = ["corum_cancer"] + [f"reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))
    network_files = {
        'biogrid': '/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/biogrid_edge_list.tsv',
        'irefindex18': '/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/irefindex18_edge_list.tsv',
        'reactome21': '/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/reactome21_edge_list.tsv',
        'string': '/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/networks/string_edge_list.tsv',
    }
    progress_bar = tqdm(total=len(data_versions) * len(network_names),desc="Processing fdrnet results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            graph = load_graph_structure(network_files[network_name])

            # for beta in betas:
            #     for permu in permus:
            progress_bar.update(1)
            file_name = f"{network_name}{data_version}_fdr_bound0.1_alpha0.85_size400_time100.0_gap0.01_conductance.txt"
            file_path = os.path.join(base_path, file_name)

            if not os.path.exists(file_path):
                continue

            with open(file_path, "r") as f:
                reader = csv.reader(f, delimiter="\t")
                local_results = {}
                for row in reader:
                    seed, time, status, subnetwork = row
                    if seed not in local_results:
                        local_results[seed] = []
                    subnetwork_genes = set(subnetwork.split(" "))
                    local_results[seed].append(subnetwork_genes)

                seed_connections = calculate_seed_connections(local_results, graph)
                filtered_results = filter_subnetworks(local_results, seed_connections)

            key = (data_version, network_name)
            results[key] = filtered_results

    progress_bar.close()
    return results


def load_domino_results(results_path,data_versions):
    results = {}
    #data_versions = [f"corum_cancer"] + [f"reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))

    progress_bar = tqdm(total=len(data_versions) * len(network_names),desc="Processing domino results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            # for beta in betas:
    #     for permu in permus:
            progress_bar.update(1)

            file_path = os.path.join(results_path, network_name, f"{data_version}_list/modules.out")

            if not os.path.exists(file_path):
                continue

            key = (data_version, network_name)

            if key not in results:
                results[key] = []

            with open(file_path, "r") as f:
                content = f.readlines()

                for line in content:
                    gene_names = set(line.strip().replace('[', '').replace(']', '').split(', '))
                    results[key].append(gene_names)

    progress_bar.close()
    return results


def load_diamond_results(results_path,data_versions):
    results = {}
    #data_versions = ["corum_cancer"] + ["reactome_cancer"]
    network_names = ["biogrid", "irefindex18", "reactome21", "string"]
    #betas = [round(i * 0.01, 2) for i in range(1, 12)]
    #permus = list(range(1, 11))
    size = 100#list(range(100, 701, 100))
    network_gene_name_mappings = load_gene_name_mapping_clustex()
    progress_bar = tqdm(total=len(data_versions) * len(network_names),desc="Processing diamond results", position=0)

    for data_version in data_versions:
        for network_name in network_names:
            progress_bar.update(1)

            files = glob.glob(os.path.join(results_path, f"{data_version}_diamond_{network_name}_{size}_module.txt"))
            seed_file = f"/Users/leyang/Dropbox/Projects/Subnetwork_survey/code/score/real/diamond/{data_version}_diamond_{network_name}.txt"

            for file in files:
                key = (data_version, network_name)
                # subnetworks = []
                # if key not in results:
                #     results[key] = []

                gene_name_mapping = network_gene_name_mappings[network_name]

                with open(file, "r") as f:
                    content = f.readlines()[1:]  # Skip the header line

                    gene_names = {gene_name_mapping.get(line.strip().split("\t")[1], line.strip().split("\t")[1]) for line in content}

                with open(seed_file, "r") as f:
                    content = f.readlines()
                    seed_names = {gene_name_mapping.get(line.strip().split("\t")[0], line.strip().split("\t")[0]) for line in content}

                r = gene_names.union(seed_names)
                results[key] = [r]

                # if size not in results[key]:
                #     results[key][size] = []

                            # Wrap the gene_names set with a list
                # results[key][size].append(r)

    progress_bar.close()
    return results

def load_bionet_results(results_path,data_versions):
    network_names = ['biogrid', 'irefindex18', 'reactome21', 'string']
    #data_versions = ['corum_cancer'] + ['reactome_cancer']
    #betas = [round(x, 2) for x in (0.01 * i for i in range(1, 12))]
    #permutations = list(range(1, 11))

    results = {}
    total_combinations = len(network_names) * len(data_versions)

    progress_bar = tqdm(total=total_combinations, desc="Processing bionet results", position=0)

    for network_name in network_names:
        for data_version in data_versions:
            # for beta in betas:
            #     for permu in permutations:
                    # Update progress bar
            progress_bar.update(1)

            # Search for files matching the current pattern
            file_pattern = f"{results_path}/{network_name}{data_version}_p_*_result.txt"
            files = glob.glob(file_pattern)

            # Load and process each file matching the pattern
            for file in files:
                # Extract k from the file name
                k = int(file.split("_")[-2])

                # Create a dictionary key for the current combination
                key = (data_version, network_name)

                # Initialize the subnetworks list if it doesn't exist
                if key not in results:
                    results[key] = []

                # Load the file and store the results in a dictionary
                with open(file, "r") as f:
                    content = f.readlines()
                    subnetwork = set()
                    for line in content[1:-1]:  # Skip the header and the last line
                        label, score = line.strip().split('\t')
                        if score != 'NaN':
                            subnetwork.add(label)

                    results[key].append(subnetwork)

    # Close progress bar
    progress_bar.close()
    return results

def load_clustex_results(results_path,data_versions):
    #data_versions = ['corum_cancer'] + ['reactome_cancer']
    network_names = ['biogrid', 'irefindex18', 'reactome21', 'string']
    #betas = [round(x, 2) for x in (0.01 * i for i in range(1, 12))]
    #permutations = list(range(1, 11))

    network_gene_name_mappings = load_gene_name_mapping_clustex()

    results = {}
    total_combinations = len(network_names) * len(data_versions)
    progress_bar = tqdm(total=total_combinations, desc="Processing clustex results",position=0)

    for data_version in data_versions:
        for network_name in network_names:
            # for beta in betas:
            #     for permu in permutations:
                    # Update progress bar
            progress_bar.update(1)

                    # Search for files matching the current pattern
            file_pattern = f"{results_path}/{data_version}_clustex_{network_name}2_clusters_50_0.03.txt"
            files = glob.glob(file_pattern)

            # Load and process each file matching the pattern
            for file in files:
                # Create a dictionary key for the current combination
                key = (data_version, network_name)

                # Initialize the subnetworks list if it doesn't exist
                if key not in results:
                    results[key] = []

                # Get the correct gene name mapping for the current network
                gene_name_mapping = network_gene_name_mappings[network_name]

                # Load the file and store the results in a dictionary
                with open(file, "r") as f:
                    content = f.read()
                    modules = re.findall(r'module \d+\tseed_gene_fraction [\d.]+\tnumber of genes: \d+\tgenes:\t(.+)', content)

                    for module in modules:
                        gene_indices = set(module.split('\t'))
                        gene_names = {gene_name_mapping.get(index, index) for index in gene_indices}
                        results[key].append(gene_names)
    # Close progress bar
    progress_bar.close()
    return results

if __name__ == "__main__":

    method = ['bionet','clustex','diamond','domino','fdrnet','hhotnet','hotnet2','must','netcore','netmix2','regmod','robust']
    func_dict = {'bionet':load_bionet_results,'clustex':load_clustex_results,'diamond':load_diamond_results,'domino':load_domino_results,
    'fdrnet':load_fdrnet_results,'hhotnet':load_hhotnet_results,'hotnet2':load_hotnet2_results,'must':load_must_results,
    'netcore':load_netcore_results,'netmix2':load_netmix2_results,'regmod':load_regmod_results,'robust':load_robust_results}
    results = {}
    data_versions = ["BLCA_MERGE_v4","LUAD_MERGE_v4","LUSC_MERGE_v4","COADREAD_MERGE_v4","PRAD_MERGE_v4","HNSC_MERGE_v4",
                 "UCEC_MERGE_v4","KIPAN_MERGE_v4","BRCA_MERGE_v4"]
    for m in method:
        path = f"/Users/leyang/Dropbox/Projects/Subnetwork_survey/all_results/{m}/result_real"
        r = func_dict[m](path,data_versions)
        results[m] = r
    save_results(results,"all_results_real.pkl")
