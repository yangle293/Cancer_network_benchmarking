from __future__ import print_function
import networkx as nx
import numpy as np
import math

# change epsilon adaptively to have a local graph large than 'size'
def approximatePPR_size_fast(graph, alpha, seed, size):
    isinstance(graph, nx.Graph)
    maxiter = 1e6
    epsilon = 1e-5
    not_converged = False
    d = np.array([1.0*val for (node, val) in graph.degree()])
    r = np.array([0.0]*len(d))
    p = np.array([0.0]*len(d))
    seed_index = list(nx.nodes(graph)).index(seed)
    r[seed_index] = 1
    resvec = r/d
    W = nx.adjacency_matrix(graph)
    
    #I = [x for (x, y) in resvec.items() if y >= epsilon]
    I = np.where(resvec > epsilon)[0]
    
    #[x for x in p.values() if x > 0.0]
    while len(np.where(p > 0.0)[0]) < size:
        
        if len(np.where(p > 0.0)[0]) != 0:
            step = size/float(len(np.where(p > 0.0)[0]))
        else:
            step = 2
        epsilon = epsilon/step
        iter_n = 0
        I = np.where(resvec > epsilon)[0]
        
        while (I.size != 0) and (iter_n < maxiter):
            k = len(I)
            iter_n = iter_n + k
            for c in I:
                n_iter = math.ceil(math.log(r[c] / (epsilon * d[c])) / math.log(2./(1-alpha))+np.finfo(float).eps)
                p, r = push_n(pp=p, res=r, uu=c, alpha_ppr=alpha, W = W,d=d, n=n_iter)
                resvec = r/d
                I = np.where(resvec > epsilon)[0]
                ind = np.argsort(resvec[I])
                I = I[ind]

    p_vector = dict(zip(nx.nodes(graph),p))
#    if iter_n > maxiter:
#        not_converged = True
    return p_vector, not_converged, r


def push_n(pp, res, uu, alpha_ppr, W,d, n):
    mult = (1 - ((1. - alpha_ppr)/2)**n) / (1 - (1 - alpha_ppr)/2)
    pp[uu] = pp[uu] + alpha_ppr * res[uu] * mult
    res = res + (1-alpha_ppr)*res[uu]/(2*d[uu])*mult*np.squeeze(W[:,uu].toarray())
    res[uu] = (1 - alpha_ppr)**n * res[uu] / (2**n)
    return pp, res


if __name__ == "__main__":
    G = nx.Graph()
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(3, 4)
    G.add_edge(4, 5)
    G.add_edge(5, 6)
    G.add_edge(6, 7)
    G.add_edge(6, 7)
    G.add_edge(7, 8)
    G.add_edge(8, 9)
    seed = 3
    alpha = 0.85
    size = 5
    p, flag, r = approximatePPR_size_fast(graph=G, seed=seed, alpha=alpha, size=size)
    print(p)
