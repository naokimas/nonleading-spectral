#    Calculate the leading eigenvector, minimizer of error e_1, minimizer of error e_2, and associated properties for the given network
#    Used for Figures 1 and 2.
#
#    Syntax
#
#    python3 stat-error-improve-spectral.py net_type nV kave alpha
#        net_type = 0:empirical (netsci-lcc), 1:PL1, 2: ER, 3: WS, 4: BA, 5: Holme-Kim
#
#    If net_type == 0, the user should manually modify the input file name
#    
#    The output, which is various statistics of the leading eigenvector, minimizer of error e_1, and minimizer of error e_2, is saved in "result.txt". 
#    Each row in result.txt corresponds to one network sample.
#    These files, appropriately renamed, are input to plot-fig-1-and-2.py
#    
#    If net_type == 0, additional output is saved in "result_2.txt".

import importlib # import a file whose name contains hyphen
import math # floor
import sys # to use input parameters
import numpy as np
import math # sqrt
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import networkx as nx

#def create_adj_mat(E):
#    """

#    Create the adjacency matrix from an edge list given by the input file 

#    Input
#        E: edge list

#    Output
#        A: adjacency matrix

#    Note
#        Node index in E must start from 1 (not from 0).

#    """

#    nV = np.int_(np.max(E[:, 0:2])+0.1**8) # E[:, 0:2] extracts the first two columns of E
#    nE = E.shape[0] # number of edges
#    A = np.zeros((nV, nV)) # initialization

#    for i in range(nE):
#        A[E[i,0].astype(int) - 1, E[i,1].astype(int) - 1] = A[E[i,1].astype(int) - 1, E[i,0].astype(int) -1] = E[i,2] # undirected network assumed
#    return A


def modified_spectral_method(A):
    """

    Calculate the performance of the one-dimensional reduction (i.e., spectral method) for each eigenvector with eigenvalue > 10^{-6}

    Input
        A: adjacency matrix

    """

    nV_tmp = A.shape[0] # number of nodes
    lam, u = np.linalg.eigh(A.T) # eigenvalues in ascending order
    # i-th column of u is the i-th normalized eigenvector
    # normalization is such that the square sum of the elements of the eigenvector = 1

    print(u[:,nV_tmp-1])
    print(lam[nV_tmp-1])
    np.savetxt('Adj_tmp.txt', A)

    kin = np.sum(A, axis = 1) # indegree. A_{ij} means j --> i. Sum of each row

    u = u[:, lam > 0.1**6] # retain eigenvectors associated with eigenvalue > 10^{-6} only
    lam = lam[lam > 0.1**6] # the corresponding eigenvalues
    num_pos_eig = lam.shape[0] # number of positive eigenvalues, precisely, eigenvalues > 10^{-6}
#    print(u.shape, lam.shape, num_pos_eig) 

    # renormalize the eigenvectors so that \sum_{i=1}^N a_i = 1
    tmp_sum = np.sum(u, axis=0) # tmp_sum[ell]: sum of all elements of the ell-th column (ell=0, 1, ..., nV_tmp - 1)
    for ell in range(num_pos_eig):
        u[:, ell] = u[:, ell] / tmp_sum[ell]
    # u[i, ell]: a_i of the ell-th eigenvector (a_1, a_2, \ldots, a_N)^{\top}

    c = u**2 # c[i, ell] = (a_i)^2 of the ell-th eigenvector

    error1 = np.zeros(num_pos_eig) # e_1
    error2 = np.zeros(num_pos_eig) # e_2
    #error4 = np.zeros(num_pos_eig) # e_4

    beta = np.zeros(num_pos_eig) # optimal beta for each eigenvector
#    num_0eigs = 0 # number of zero eigenvalues
    for ell in range(num_pos_eig): # ell-th eigenvalue/vector
#        if np.absolute(lam[ell]) < 0.1**6: # eigenvalue = 0
#            num_0eigs = num_0eigs + 1
#            beta[ell] = 10**8 # dummy
#            error1[ell] = 10**8 # dummy
#            error2[ell] = 10**8 # dummy
    #        error4[ell] = 10**8 # dummy
#        else: # eigenvalue \neq 0
        b = c[:, ell] / np.sum(c[:, ell]) # b_i in Laurence et al. Phys. Rev. X (2019), p.15, right column
        beta[ell] = np.dot(b, kin) / np.dot(u[:, ell], kin) # beta^* (i.e., optimal beta) in Eq. (A5), p.15, in Laurence et al. Phys. Rev. X (2019) 
        numerator = 0.0
        for i in range(nV_tmp):
            for j in range(i):
                numerator = numerator + c[i, ell] * c[j, ell] * (kin[i] - kin[j])**2
        error1[ell] = numerator / np.sum(c[:, ell]) # e_1
        error2[ell] = np.dot((kin - lam[ell])**2, c[:, ell]) # e_2
        #    error4[ell] = np.sum(A[ell, :]**2) + lam[ell]**2 # e_4
# errors calculated

    tmp = error1.argsort() # tmp[i]-th entry (index starting from 0) of error1 is the i-th smallest value of error1 (i value starting from 0)
    min_error1 = error1[tmp[0]]
    eig_val_optimal1 = lam[tmp[0]] # eigenvalue minimizinag e_1
    eig_vec_optimal1 = u[:, tmp[0]] # eigenvector minimizing e_1
    beta_optimal1 = beta[tmp[0]] # beta^* corresponding to the minimizer of e_1
    ranks_error1 = tmp.argsort() # rank of the leading eigenvector in terms of e_1

    tmp = error2.argsort() # similar for e_2
    min_error2 = error2[tmp[0]]
    eig_val_optimal2 = lam[tmp[0]]
    eig_vec_optimal2 = u[:, tmp[0]]
    ranks_error2 = tmp.argsort()

    return beta[num_pos_eig - 1], lam[num_pos_eig - 1], u[:, num_pos_eig - 1], min_error1, error1[num_pos_eig - 1], ranks_error1[num_pos_eig - 1], beta_optimal1, eig_val_optimal1, eig_vec_optimal1, min_error2, error2[num_pos_eig - 1], ranks_error2[num_pos_eig - 1], eig_val_optimal2, eig_vec_optimal2
#end of modified_spectral_method

# Main starts

param = sys.argv
argc = len(param)
np.set_printoptions(precision=3, suppress=True) # number of digits after the decimal point

net_type = int(param[1])
if net_type == 0: # largest connected component of the coauthorship network. It can be any network, and a user should manually modify the input file name below.
    samples = 1 # a single network
if net_type >= 1: # a model network
    samples = 200 # number of networks to be generated
    nV = int(param[2]) # number of nodes
    kave = float(param[3]) # average degree
    if net_type == 1: # PL1
        alpha = float(param[4]) # power-law exponent (denoted by gamma in the paper) = alpha + 1
        kappa = 1.0 / (kave - 1) / (alpha - 1)

# initialize vectors
eig_val_leading = np.zeros(samples)
eig_val_optimal1 = np.zeros(samples)
eig_val_optimal2 = np.zeros(samples)
beta_leading = np.zeros(samples)
beta_optimal1 = np.zeros(samples)
ratio_error1 = np.zeros(samples)
ratio_error2 = np.zeros(samples)

if net_type >= 1: # model network
    nV_actual_min = nV # smallest number of nodes in the largest connected component across network samples
    
for ind in range(samples):

    if net_type == 0: # empirical network
        Graph = nx.read_weighted_edgelist('netscience-lcc.mat', nodetype=int, comments='%',)
        # Format of the network input data is as follows.
        # The first row has "% nV nE"
        # where nV is the number of nodes, and nE is the number of edges.
        # The remaining rows have edges, one edge per row.
        # Each row has three entries (columns). The first two columns represent the two nodes connected by the edge. They are integers \in \{1, ..., nV\}. Note that the node index should start from 1, not from 0. The third column represents the edge weight.
        nV = Graph.number_of_nodes()
        nV_actual_min = nV
    elif net_type == 1: # configuration model with a power-law degree distribution generated by PL1
        if alpha > 0: # configuration model with a discretized PL1 power-law degree distribution. Power-law exponent gamma = alpha+1
            degree_sum = -1 # generate the degree sequence
            while (degree_sum % 2 == 1): # repeat the following until the degree sum becomes an even number so that we can run the configuration model
                u = 1 - np.random.uniform(0, 1, nV) # random.uniform produces u \in [0, 1)
                k_float = 1 + (u**(-1/alpha) - 1) / kappa # determine degree (in float) using inverse sampling
                degree_sequence = np.rint(k_float).astype(np.int32) # rounding to the nearest integer
                degree_sum = np.sum(degree_sequence)
        # print(degree_sequence)
        # print(np.average(degree_sequence))
        #degree_sequence = [d for n, d in G.degree()]  # degree sequence
        # print(f"Degree sequence {degree_sequence}")
        # https://networkx.org/documentation/stable/auto_examples/graph/plot_degree_sequence.html
        # degree sequence determined
        pseudoGraph=nx.configuration_model(degree_sequence)
        Graph = nx.Graph(pseudoGraph) # remove multiple edges
        Graph.remove_edges_from(nx.selfloop_edges(Graph)) # remove self-loops
    elif net_type == 2: # Erdos-Renyi
        Graph = nx.fast_gnp_random_graph(nV, kave/(nV-1))
    elif net_type == 3: # Watts-Strogatz
        pseudoGraph = nx.watts_strogatz_graph(nV, int(kave), 0.1)
        # the documentation 
        # https://networkx.org/documentation/stable/reference/generated/networkx.generators.random_graphs.watts_strogatz_graph.html
        # does not say whether the generated network has multiple edges. So, just in case, I remove it in the following command.
        Graph = nx.Graph(pseudoGraph) # remove multiple edges
    elif net_type == 4: # Barabasi-Albert (not used in the paper)
        Graph = nx.barabasi_albert_graph(nV, int(kave/2)) # the initial graph is the star with m+1 nodes
        # m = kave/2
    elif net_type == 5: # Holme-Kim
        Graph = nx.powerlaw_cluster_graph(nV, int(kave/2), 0.5)

    nV_inputnet = Graph.number_of_nodes()
    nE_inputnet = Graph.number_of_edges()
    print(nV_inputnet, 'nodes,', nE_inputnet, 'edges')
    largest_cc = max(nx.connected_components(Graph), key=len) # largest connected component
    # https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.connected_components.html?highlight=largest%20connected
    G = Graph.subgraph(largest_cc).copy()

    nV_actual = G.number_of_nodes() # number of nodes in the largest connected component

    # https://stackoverflow.com/questions/62427114/ordered-nodes-in-adjacency-matrix-when-using-nx-to-numpy-array
    A = nx.to_numpy_array(G, nodelist=sorted(G.nodes())) # transform to the adjacency matrix

    print(ind, ':', int(nV_actual), 'nodes,', int(G.number_of_edges()), 'edges')

    if nV_actual < nV_actual_min:
        nV_actual_min = nV_actual

    if net_type == 0: # network given by the input file
        A = A.transpose()

    beta_leading[ind], eig_val_leading[ind], eig_vec_leading, min_error1, v_error1, v_rank1, beta_optimal1[ind], eig_val_optimal1[ind], eig_vec_optimal1, min_error2, v_error2, v_rank2, eig_val_optimal2[ind], eig_vec_optimal2 = modified_spectral_method(A)

    ratio_error1[ind] = min_error1 / v_error1
    ratio_error2[ind] = min_error2 / v_error2

# All samples done

# ratio_error1.sort()
# ratio_error2.sort()

print(nV_actual_min)

# print(np.sum(u, axis=0))

# np.savetxt('beta.txt', beta_list, fmt='%10.6f')
        
x = np.vstack((ratio_error1, ratio_error2, beta_leading, beta_optimal1, eig_val_leading, eig_val_optimal1, eig_val_optimal2)).T
np.savetxt('result.txt', x, fmt=['%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f'])

if net_type==0:
    y = np.vstack((eig_vec_leading, eig_vec_optimal2)).T
    np.savetxt('result_2.txt', y, fmt=['%10.6f', '%10.6f'])
