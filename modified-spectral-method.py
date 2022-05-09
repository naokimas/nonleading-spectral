#    Generate networks and determine the leading eigenvector, minimizer of e_1, and minimizer of e_2 for Figs. 3-7
#
#    Syntax
#
#    python modified-spectral-method.py nV kave alpha
#        If alpha < 0, ER.
#        If nV < 0, a single network given below. In this case, the node index starts from 1 in the input file.
#
#    Outout:
#
#       E_i.txt: edge list of the i-th network instance
#       a_i.txt: three eigenvectors for the i-th network instance
#           The first column is the leading eigenvector.
#           The second column is the minimizer of e_1.
#           The third column is the minimizer of e_2.
#       out_metadata_i.txt: summary result for the i-th network instance
#       result.txt: summary numbers etc.    

import importlib # import a file whose name contains hyphen
import math # floor
import sys # to use input parameters
import numpy as np
import math # sqrt
# from scipy.sparse import csgraph # Laplacian
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import networkx as nx

def modified_spectral_method(A):
    """

    Calculate the performance of the one-dimensional reduction for each eigenvector with eigenvalue > 10^{-6}

    Input
        A: adjacency matrix

    """

    nV_tmp = A.shape[0]
    lam, u = np.linalg.eigh(A.T) # eigenvalues in ascending order
    # i-th column of u is the i-th normalized eigenvector
    # normalization is such that the square sum of the elements of the eigenvector = 1

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

    tmp = error1.argsort() # tmp[i]-th entry (index starting from 0) of error1 is the i-th smallest value of errro1 (i value starting from 0)
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

nV = int(param[1])
if nV > 0:
    kave = float(param[2])
    alpha = float(param[3]) # power-law exponent = alpha + 1
    if alpha > 0: # power-law distribution
        kappa = 1.0 / (kave - 1) / (alpha - 1)
    samples = 200 # number of networks to be generated
else: # nV < 0. Use a single network given by infilename.mat
    kave = -1 # dummy
    alpha = -1 # dummy
    samples = 1 # a single network

# initialize vectors
nV_actual = np.zeros(samples)
nE_actual = np.zeros(samples)
beta_leading = np.zeros(samples)
beta_optimal1 = np.zeros(samples)
eig_val_leading = np.zeros(samples)
eig_val_optimal1 = np.zeros(samples)
eig_val_optimal2 = np.zeros(samples)
min_error1 = np.zeros(samples)
v_error1 = np.zeros(samples)
v_rank1 = np.zeros(samples)
min_error2 = np.zeros(samples)
v_error2 = np.zeros(samples)
v_rank2 = np.zeros(samples)

for ind in range(samples):

    if nV > 0: # either scale-free networks or the ER random graph
        if alpha > 0: # configuration model with a discretized PL1 power-law degree distribution. Power-law exponent = alpha+1
            degree_sum = -1
            while (degree_sum % 2 == 1): # repeat the following until the degree sum becomes an even number so that we can run the configuration model
                u = 1 - np.random.uniform(0, 1, nV) # random.uniform produces u \in [0, 1)
                k_float = 1 + (u**(-1/alpha) - 1) / kappa # determine degree (in float) using inverse sampling
                degree_sequence = np.rint(k_float).astype(np.int32) # rounding to the nearest integer
                degree_sum = np.sum(degree_sequence)

            # degree sequence determined
            pseudoGraph=nx.configuration_model(degree_sequence)
            Graph = nx.Graph(pseudoGraph) # remove multiple edges

        else: # alpha < 0. Erdos-Renyi random graph
            Graph = nx.fast_gnp_random_graph(nV, kave/(nV-1))

#    m = 4 # new edges per new node in the Barabasi-Albert model
#    Graph = nx.barabasi_albert_graph(nV, m) # the initial graph is the star with m+1 nodes
# Graph = nx.powerlaw_cluster_graph(500, 3, 0.5) # Holme-Kim model
    else: # nV < 0. A single network given by the input file.
        Graph = nx.read_weighted_edgelist('sociop-agg.mat', nodetype=int, comments='%',)
# This should be replaced by your input file name
        nV_inputnet = Graph.number_of_nodes()
        nE_inputnet = Graph.number_of_edges()
        print(nV_inputnet, 'nodes,', nE_inputnet, 'edges')
        # columns 1 and 2: two nodes to form an edge
        # column 3: edge weight
        # E = np.loadtxt(param[1], skiprows=1)
        # E = np.loadtxt(param[1] + '.mat', skiprows=1)

    # Regardless of how Graph has been created,
    # remove self edges
    Graph.remove_edges_from(nx.selfloop_edges(Graph))

    largest_cc = max(nx.connected_components(Graph), key=len)
    # https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.connected_components.html?highlight=largest%20connected

    G = Graph.subgraph(largest_cc).copy() # largest connected component

    nV_actual[ind] = G.number_of_nodes() # number of nodes in the largest connected component
    nE_actual[ind] = G.number_of_edges() # number of edges in the largest connected component
    A = nx.to_numpy_array(G) # transform to the adjacency matrix

    print(ind, ':', int(nV_actual[ind]), 'nodes,', int(nE_actual[ind]), 'edges')

    if nV < 0: # network given by the input file
        A = A.transpose()

    beta_leading[ind], eig_val_leading[ind], eig_vec_leading, min_error1[ind], v_error1[ind], v_rank1[ind], beta_optimal1[ind], eig_val_optimal1[ind], eig_vec_optimal1, min_error2[ind], v_error2[ind], v_rank2[ind], eig_val_optimal2[ind], eig_vec_optimal2 = modified_spectral_method(A)

    if min_error1[ind] / v_error1[ind] < 0.7 and min_error2[ind] / v_error2[ind] < 0.7:
        a_list = np.vstack((eig_vec_leading, eig_vec_optimal1, eig_vec_optimal2)).T
        np.savetxt('a_' + str(ind) + '.txt', a_list, fmt='%10.6f') # three eigenvectors
        nx.write_edgelist(G, 'E_' + str(ind) + '.txt', delimiter=' ', data = False)
        outfile_metadata = open('out_metadata_' + str(ind) + '.txt', 'w')
        outfile_metadata.write('beta for leading eigenvector = ' + str(beta_leading[ind]) + '\n')
        outfile_metadata.write('beta optimal 1= ' + str(beta_optimal1[ind]) + '\n')
        outfile_metadata.write('error e_1 with leading eigenvector = ' + str(v_error1[ind]) + '\n')
        outfile_metadata.write('error e_1 optimized = ' + str(min_error1[ind]) + '\n')
        outfile_metadata.write('error e_2 with leading eigenvector = ' + str(v_error2[ind]) + '\n')
        outfile_metadata.write('error e_2 optimized = ' + str(min_error2[ind]) + '\n')
        outfile_metadata.write('eigenvalue for leading eigenvector = ' + str(eig_val_leading[ind]) + '\n')
        outfile_metadata.write('eigenvalue for e_1 optimizer = ' + str(eig_val_optimal1[ind]) + '\n')
        outfile_metadata.write('eigenvalue for e_2 optimizer = ' + str(eig_val_optimal2[ind]))
        outfile_metadata.close()

# beta_list = np.vstack((beta_best_classical, beta_best))
# print(min_error1 / error1[nV-1], ranks_error1[nV-1], min_error2 / error2[nV-1], ranks_error2[nV-1])
# print(np.sum(u, axis=0))
# np.savetxt('beta.txt', beta_list, fmt='%10.6f')

x = np.vstack((nV_actual, nE_actual, beta_leading, eig_val_leading, min_error1, v_error1, v_rank1, beta_optimal1, eig_val_optimal1, min_error2, v_error2, v_rank2, eig_val_optimal2)).T

np.savetxt('result.txt', x, fmt=['%d', '%d', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%d', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%d', '%10.3f'])