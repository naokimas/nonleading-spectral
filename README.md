# nonleading-spectral

This is a code repository for our paper on a modified spectral method to reduce dynamical systems on networks into a one-dimensional dynamics. When you use the code, please cite the following paper:

Naoki Masuda, Prosenjit Kundu.
Dimension reduction of dynamical systems on networks with leading and non-leading eigenvectors of adjacency matrices.
Preprint: arXiv:2203.13872

The code is a mix of Python and MATLAB code. Naoki Masuda developed the Python codes, and Prosenjit Kundu developed the MATLAB codes. The codes produce the numerical results shown as the figures in the paper.

Figures 1 and 2 are generated by first running stat-modified-spectral-method.py and then running plot-fig-1-and-2.py. See comments in those files for detailed instructions.

For Figures 3-7, the networks, leading eigenvector, minimizer of $$\epsilon_1$$, minimizer of $$\epsilon_2$$, the associated eigenvalues, and $$\beta^*$$ are generated by modified-spectral-method.py

Then, the (stochastic) dynamical systems are simulated with MATLAB code as follows.

Panel (a), (b) and (e) of Figures 3, 4 and 5 can be generated using run_fig345_a_b_e.m. 

Panel (c) and (d) will be generated using run_fig345_c_d.m

Figure 5

Figure 6

Figure 7
