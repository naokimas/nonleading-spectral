# generate figures 1 and 2 in the paper
#
# The output figure is saved in fig.eps

import numpy as np
import matplotlib.pyplot as plt
from scipy import linspace

data = 4
# data == 0: produce Fig. 1(a) in the paper
# data == -1: Fig. 1(b)
# data == 1: Fig. 2(a)
# data == 2: Fig. 2(b)
# data == 3: Fig. 2(c)
# data ==4: Fig. 2(d)
# data == 7: Fig. 2(e)
# data == 8: Fig. 2(f)

# necessary to get the fonts right
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

# calculate statistics to plot from result.txt (unused)
def calc_stat(r):
    samples = r.shape[0]
    nV_actual = r[:, 0]
    min_error1 = r[:, 4]
    laurence_error1 = r[:, 5]
    laurence_rank1 = r[:, 6]
    min_error2 = r[:, 9]
    laurence_error2 = r[:, 10]
    laurence_rank2 = r[:, 11]

    nV_actual_ave = np.average(nV_actual)
    nV_actual_std = np.std(nV_actual)
    fraction_not_top1 = float(np.count_nonzero(laurence_rank1))/samples # fraction of samples for which the leading eigenvector does not minimize error 1
    fraction_not_top2 = np.count_nonzero(laurence_rank2)/samples # same for error 2

#    nonzero_relative_rank_in_error1 = laurence_rank1[laurence_rank1 > 0] / (nV - num_0eigs[laurence_rank1 > 0])
#    nonzero_relative_rank_in_error2 = laurence_rank2[laurence_rank2 > 0] / (nV - num_0eigs[laurence_rank2 > 0])

    reduction_in_error1 = min_error1[laurence_rank1 > 0] / laurence_error1[laurence_rank1 > 0]
    reduction_in_error2 = min_error2[laurence_rank2 > 0] / laurence_error2[laurence_rank2 > 0]
    reduction_in_error1_ave = np.average(reduction_in_error1)
    reduction_in_error1_std = np.std(reduction_in_error1)
    reduction_in_error2_ave = np.average(reduction_in_error2)
    reduction_in_error2_std = np.std(reduction_in_error2)

    return nV_actual_ave, nV_actual_std, fraction_not_top1, fraction_not_top2, reduction_in_error1_ave, reduction_in_error1_std, reduction_in_error2_ave, reduction_in_error2_std

# "ratio-...txt" files contain the leading eigenvector, error-minimizing eigenvector, beta^*, etc.
if data==0 or data==-1: # scale-free networks from a configuration model
    r = np.loadtxt("ratio-n1000k10gamma3.5-detailed.txt") # scale-free network with N = 1000, <k> = 10, and \tilde{\gamma} = 3.5
elif data==1 or data==2: # scale-free networks from a configuration model
    r1 = np.loadtxt("ratio-n100k4gamma2.5.txt") # N = 100, <k> = 4, \tilde{\gamma} = 2.5
    r2 = np.loadtxt("ratio-n100k10gamma2.5.txt") # N = 100, <k> = 10, \tilde{\gamma} = 2.5
    r3 = np.loadtxt("ratio-n1000k4gamma2.5.txt") # N = 1000, <k> = 4, \tilde{\gamma} = 2.5
    r4 = np.loadtxt("ratio-n1000k10gamma2.5.txt") # N = 1000, <k> = 10, \tilde{\gamma} = 2.5
elif data==3 or data==4: # scale-free networks from a configuration model
    r1 = np.loadtxt("ratio-n100k4gamma3.5.txt") # N = 100, <k> = 4, \tilde{\gamma} = 3.5
    r2 = np.loadtxt("ratio-n100k10gamma3.5.txt") # N = 100, <k> = 10, \tilde{\gamma} = 3.5
    r3 = np.loadtxt("ratio-n1000k4gamma3.5.txt") # N = 1000, <k> = 4, \tilde{\gamma} = 3.5
    r4 = np.loadtxt("ratio-n1000k10gamma3.5.txt") # N = 1000, <k> = 10, \tilde{\gamma} = 3.5
elif data==5 or data==6: # unused in the paper
    r1 = np.loadtxt("ratio-n100k4ba.txt")
    r2 = np.loadtxt("ratio-n100k10ba.txt")
    r3 = np.loadtxt("ratio-n1000k4ba.txt")
    r4 = np.loadtxt("ratio-n1000k10ba.txt")
elif data==7 or data==8: # Holme-Kim model
    r1 = np.loadtxt("ratio-n100k4hk.txt") # N = 100, <k> = 4
    r2 = np.loadtxt("ratio-n100k10hk.txt") # N = 100, <k> = 10
    r3 = np.loadtxt("ratio-n1000k4hk.txt") # N = 1000, <k> = 4
    r4 = np.loadtxt("ratio-n1000k10hk.txt") # N = 1000, <k> = 10

if data >= 1:
    Nsamples = r1.shape[0]
    x1 = np.linspace(0, 1.0, Nsamples, endpoint=False) + 1.0/Nsamples # x1 = [1/Nsamples, 2/Nsamples, ..., 1] 
    # vector of x values, ranging from 0 to 1

if data==0: # \beta^* for the leading eigenvector
    plt.plot(r[:,0], r[:,2], 'ko', linewidth=2.5, markerfacecolor='none')
elif data==-1: # \beta^* for the minimizer of e_1
    plt.semilogy(r[:,0], r[:,3], 'ko', linewidth=2.5, markerfacecolor='none')
elif data==1 or data == 3 or data == 5 or data == 7: # minimizer of e_1
    rtmp = r1[:,0]
    rtmp.sort()
    plt.plot(x1, rtmp, 'k-', label=r'$N=100, \langle k\rangle = 4$', linewidth=2.5)
    rtmp = r2[:,0]
    rtmp.sort()
    plt.plot(x1, rtmp, 'r-', label=r'$N=100, \langle k\rangle = 10$', linewidth=2.5)
    rtmp = r3[:,0]
    rtmp.sort()
    plt.plot(x1, rtmp, 'k:', label=r'$N=1000, \langle k\rangle = 4$', linewidth=2.5)
    rtmp = r4[:,0]
    rtmp.sort()
    plt.plot(x1, rtmp, 'r:', label=r'$N=1000, \langle k\rangle = 10$', linewidth=2.5)
elif data==2 or data == 4 or data == 6 or data == 8: # minimizer of e_2
    rtmp = r1[:,1]
    rtmp.sort()
    plt.plot(x1, rtmp, 'k-', label=r'$N=100, \langle k\rangle = 4$', linewidth=2.5)
    rtmp = r2[:,1]
    rtmp.sort()
    plt.plot(x1, rtmp, 'r-', label=r'$N=100, \langle k\rangle = 10$', linewidth=2.5)
    rtmp = r3[:,1]
    rtmp.sort()
    plt.plot(x1, rtmp, 'k:', label=r'$N=1000, \langle k\rangle = 4$', linewidth=2.5)
    rtmp = r4[:,1]
    rtmp.sort()
    plt.plot(x1, rtmp, 'r:', label=r'$N=1000, \langle k\rangle = 10$', linewidth=2.5)

plt.xlabel('optimal/leading', fontsize=24)
if data >= 1: # Fig. 2
    plt.ylabel('cumulative distribution', fontsize=24)
    plt.legend(loc = 'upper left', numpoints = 1, labelspacing=0.25, frameon=False)
    plt.xlim(0, 1)
    plt.xticks(np.linspace(0, 1.0, 6, endpoint=True), ('0', '0.2', '0.4', '0.6', '0.8', '1'))
    plt.ylim(0, 1.01)
    plt.yticks(np.linspace(0, 1.0, 6, endpoint=True), ('0', '0.2', '0.4', '0.6', '0.8', '1'))
else: # Fig. 1
    plt.xlim(0., 1.01)
    plt.xticks(np.linspace(0, 1.0, 6, endpoint=True), ('0', '0.2', '0.4', '0.6', '0.8', '1'))
    plt.ylabel(r'$\beta^{*}$', fontsize=24, rotation=0, labelpad=12)

if data==0:
    plt.title('(a)', fontsize=28, x=-0.08, y=1.05)
    plt.ylim(1, 6)
    plt.yticks(np.linspace(1, 6, 6, endpoint=True), ('1', '2', '3', '4', '5', '6'))
elif data==-1:
    plt.ylim(1, 600000)
    plt.title('(b)', fontsize=28, x=-0.08, y=1.05)
    plt.yticks([1, 10, 100, 1000, 10000, 100000], (r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'), fontname="Arial")
elif data==1:
    plt.title(r'(a) $\tilde{\gamma} = 2.5, \epsilon_1$', fontsize=28, x=0.1, y=1.03)
elif data==2:
    plt.title(r'(b) $\tilde{\gamma} = 2.5, \epsilon_2$', fontsize=28, x=0.1, y=1.03)
    plt.legend('', frameon=False) # remove legend
elif data==3:
    plt.title(r'(c) $\tilde{\gamma} = 3.5, \epsilon_1$', fontsize=28, x=0.1, y=1.03)
    plt.legend('', frameon=False) # remove legend
elif data==4:
    plt.title(r'(d) $\tilde{\gamma} = 3.5, \epsilon_2$', fontsize=28, x=0.1, y=1.03)
    plt.legend('', frameon=False) # remove legend
elif data==5: # unused in the paper
    plt.title(r'(e-dummy) BA, $\epsilon_1$', fontsize=28, x=0.25, y=1.03)
    plt.legend('', frameon=False) # remove legend
elif data==6: # unused in the paper
    plt.title(r'(f-dummy) BA, $\epsilon_2$', fontsize=28, x=0.25, y=1.03)
    plt.legend('', frameon=False) # remove legend
elif data==7:
    plt.title(r'(e) Holme-Kim, $\epsilon_1$', fontsize=28, x=0.175, y=1.03)
    plt.legend('', frameon=False) # remove legend
elif data==8:
    plt.title(r'(f) Holme-Kim, $\epsilon_2$', fontsize=28, x=0.175, y=1.03)
    plt.legend('', frameon=False) # remove legend

plt.subplots_adjust(bottom=0.16, left=0.15, right=0.95)
plt.tick_params(labelsize=20)

plt.savefig("fig.eps")