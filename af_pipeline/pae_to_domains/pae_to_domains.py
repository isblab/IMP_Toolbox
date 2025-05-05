# Author: Tristan Croll
# Github: https://github.com/tristanic/pae_to_domains

from random import random


def parse_pae_file(pae_file):
    #! modified to handle both AF2 and AF3 PAE formats

    import numpy, os
    ext = os.path.splitext(pae_file)[1]

    if ext == '.json':
        import json
        with open(pae_file, 'r') as f:
            data = json.load(f)
    elif ext == '.pkl':
        import pickle
        with open(pae_file, 'rb') as f:
            data = pickle.load(f)
    else:
        raise Exception('Invalid file format. Must be either JSON or pickle.')

    if isinstance(data, list):
        data = data[0]

    if 'pae' in data:
        # AF3
        matrix = numpy.array(data['pae'], dtype=numpy.float64)
        token_res_ids = numpy.array(data['token_res_ids'], dtype=numpy.float64)

        diffs = numpy.ediff1d(token_res_ids) # If the diffs is 0, there are consecutive values
        boundaries = numpy.concatenate(([0], (numpy.where(diffs != 0)[0] + 1), [len(token_res_ids)])) 

        # Collect ranges where repetition length is >= 2
        for start, end in zip(boundaries[:-1], boundaries[1:]):
            if end - start >= 2:
                mean_value = numpy.round(numpy.mean(matrix[start:end, start:end]), 2)
                matrix[start, start] = mean_value  # Set the matrix value to mean value
                indices_to_remove = (range(start + 1, end))
                matrix = numpy.delete(numpy.delete(matrix, indices_to_remove, axis=0), indices_to_remove, axis=1)

    elif "predicted_aligned_error" in data:
        # AF2
        matrix = numpy.array(data["predicted_aligned_error"], dtype=numpy.float64)

    #! Legacy PAE format
    # with open(pae_json_file, 'rt') as f:
    #     data = json.load(f)[0]

    # if 'residue1' in data and 'distance' in data:
    #     # Legacy PAE format, keep for backwards compatibility.
    #     r1, d = data['residue1'], data['distance']
    #     size = max(r1)
    #     matrix = numpy.empty((size, size), dtype=numpy.float64)
    #     matrix.ravel()[:] = d

    else:
        raise ValueError('Invalid PAE JSON format.')

    return matrix

def domains_from_pae_matrix_label_propagation(pae_matrix, pae_power=1, pae_cutoff=5, random_seed=1):
    try:
        import networkx
    except ImportError:
        print('ERROR: This method requires NetworkX to be installed. Please install it using "pip install networkx" '
            'in a Python >=3.6 environment and try again.')
        import sys
        sys.exit()
    import numpy
    weights = 1/pae_matrix**pae_power
    
    g = networkx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    print(pae_cutoff)
    edges = numpy.argwhere(pae_matrix < pae_cutoff)
    sel_weights = weights[edges.T[0], edges.T[1]]
    wedges = [(i,j,w) for (i,j),w in zip(edges,sel_weights)]
    g.add_weighted_edges_from(wedges)

    from networkx.algorithms import community
    
    clusters = list(community.fast_label_propagation_communities(g, weight='weight' ,seed=random_seed))
    clusters = [list(c) for c in clusters]
    # print(clusters)
    return clusters

def domains_from_pae_matrix_networkx(pae_matrix, pae_power=1, pae_cutoff=5, graph_resolution=1):
    '''
    Takes a predicted aligned error (PAE) matrix representing the predicted error in distances between each 
    pair of residues in a model, and uses a graph-based community clustering algorithm to partition the model
    into approximately rigid groups.

    Arguments:

        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should be set to some non-zero
          value to avoid divide-by-zero warnings
        * pae_power (optional, default=1): each edge in the graph will be weighted proportional to (1/pae**pae_power)
        * pae_cutoff (optional, default=5): graph edges will only be created for residue pairs with pae<pae_cutoff
        * graph_resolution (optional, default=1): regulates how aggressively the clustering algorithm is. Smaller values
          lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

    Returns: a series of lists, where each list contains the indices of residues belonging to one cluster.
    '''
    try:
        import networkx as nx
    except ImportError:
        print('ERROR: This method requires NetworkX (>=2.6.2) to be installed. Please install it using "pip install networkx" '
            'in a Python >=3.7 environment and try again.')
        import sys
        sys.exit()
    import numpy
    weights = 1/pae_matrix**pae_power

    g = nx.Graph()
    size = weights.shape[0]
    g.add_nodes_from(range(size))
    edges = numpy.argwhere(pae_matrix < pae_cutoff)
    sel_weights = weights[edges.T[0], edges.T[1]]
    wedges = [(i,j,w) for (i,j),w in zip(edges,sel_weights)]
    g.add_weighted_edges_from(wedges)

    from networkx.algorithms import community

    clusters = community.greedy_modularity_communities(g, weight='weight', resolution=graph_resolution)
    return clusters

def domains_from_pae_matrix_igraph(pae_matrix, pae_power=1, pae_cutoff=5, graph_resolution=1):
    '''
    Takes a predicted aligned error (PAE) matrix representing the predicted error in distances between each 
    pair of residues in a model, and uses a graph-based community clustering algorithm to partition the model
    into approximately rigid groups.

    Arguments:

        * pae_matrix: a (n_residues x n_residues) numpy array. Diagonal elements should be set to some non-zero
          value to avoid divide-by-zero warnings
        * pae_power (optional, default=1): each edge in the graph will be weighted proportional to (1/pae**pae_power)
        * pae_cutoff (optional, default=5): graph edges will only be created for residue pairs with pae<pae_cutoff
        * graph_resolution (optional, default=1): regulates how aggressively the clustering algorithm is. Smaller values
          lead to larger clusters. Value should be larger than zero, and values larger than 5 are unlikely to be useful.

    Returns: a series of lists, where each list contains the indices of residues belonging to one cluster.
    '''
    try:
        import igraph
    except ImportError:
        print('ERROR: This method requires python-igraph to be installed. Please install it using "pip install python-igraph" '
            'in a Python >=3.6 environment and try again.')
        import sys
        sys.exit()
    import numpy
    weights = 1/pae_matrix**pae_power

    g = igraph.Graph()
    size = weights.shape[0]
    g.add_vertices(range(size))
    edges = numpy.argwhere(pae_matrix < pae_cutoff)
    sel_weights = weights[edges.T[0], edges.T[1]]
    g.add_edges(edges)
    g.es['weight']=sel_weights

    vc = g.community_leiden(weights='weight', resolution_parameter=graph_resolution/100, n_iterations=-1)
    membership = numpy.array(vc.membership)
    from collections import defaultdict
    clusters = defaultdict(list)
    for i, c in enumerate(membership):
        clusters[c].append(i)
    clusters = list(sorted(clusters.values(), key=lambda l:(len(l)), reverse=True))
    return clusters



_defaults = {
    'output_file':  'clusters.csv',
    'pae_power':    1.0,
    'pae_cutoff':   5.0,
    'resolution':   1.0,
    'library':      'igraph'
}

if __name__ == '__main__':
    from time import time
    start_time = time()
    import argparse
    parser = argparse.ArgumentParser(description='Extract pseudo-rigid domains from an AlphaFold PAE matrix.')
    parser.add_argument('pae_file', type=str, help="Name of the PAE JSON file.")
    parser.add_argument('--output_file', type=str, default=_defaults['output_file'], help=f'Name of output file (comma-delimited text format. Default: {_defaults["output_file"]}')
    parser.add_argument('--pae_power', type=float, default=_defaults['pae_power'], help=f'Graph edges will be weighted as 1/pae**pae_power. Default: {_defaults["pae_power"]}')
    parser.add_argument('--pae_cutoff', type=float, default=_defaults['pae_cutoff'], help=f'Graph edges will only be created for residue pairs with pae<pae_cutoff. Default: {_defaults["pae_cutoff"]}')
    parser.add_argument('--resolution', type=float, default=_defaults['resolution'], help=f'Higher values lead to stricter (i.e. smaller) clusters. Default: {_defaults["resolution"]}')
    parser.add_argument('--library', type=str, default=_defaults['library'], help=f'Graph library to use. "igraph" is about 40 times faster; "networkx" is pure Python. Default: {_defaults["library"]}')
    args = parser.parse_args()
    pae = parse_pae_file(args.pae_file)
    lib = args.library
    if lib=='igraph':
        f = domains_from_pae_matrix_igraph
    else:
        f = domains_from_pae_matrix_networkx
    clusters = f(pae, pae_power=args.pae_power, pae_cutoff=args.pae_cutoff, graph_resolution=args.resolution)
    max_len = max([len(c) for c in clusters])
    clusters = [list(c) + ['']*(max_len-len(c)) for c in clusters]
    output_file = args.output_file
    with open(output_file, 'wt') as outfile:
        for c in clusters:
            outfile.write(','.join([str(e) for e in c])+'\n')
    end_time = time()
    print(f'Wrote {len(clusters)} clusters to {output_file}. Biggest cluster contains {max_len} residues. Run time was {end_time-start_time:.2f} seconds.')
    
