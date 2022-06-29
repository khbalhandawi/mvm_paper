import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

from man_defs import get_man_combined, get_man
from mvm import nearest
from mvm.DOELib import scaling
from mvm.utilities import check_folder

mpl.rcParams['axes.linewidth'] = 0.5

# get the man object for the fea problem and load it
base_folder = os.path.join('data','strut_fea','opt')
img_folder = os.path.join('images','strut_fea')
check_folder(img_folder)
man = get_man_combined()
man_name = 'strut_comb'
man.load(man_name,folder=base_folder)
# for fea problem
lb = np.array([-0.009630757, 1.152203562])
ub = np.array([0.3329746, 4.3151886])

# # get the man object and load it
# base_folder = os.path.join('data','strut','opt')
# img_folder = os.path.join('images','strut')
# check_folder(img_folder)
# man = get_man()
# man_name = 'strut_s'
# man.load(man_name,folder=base_folder)
# # for simple problem
# lb = np.array([0.00, 2.5])
# ub = np.array([0.075, 3.5])

############################################################
# Optimization sweep
lean_angles = np.linspace(0.0,50.0,50)
n_designs = len(lean_angles)
n_epochs = 1000
# #---------------------------------------------------
# # Evaluate all the different designs
# for i,lean_angle in enumerate(lean_angles):

#     man.nominal_design_vector = np.array([15.0,lean_angle])
#     man.reset()
#     man.reset_outputs()
#     print(man.design_vector)

#     # Perform Monte-Carlo simulation
#     for n in range(n_epochs):
#         sys.stdout.write("Progress: %d%%   \r" % ((n / n_epochs) * 100))
#         sys.stdout.flush()

#         man.randomize()
#         man.init_decisions()
#         man.forward()
#         man.compute_impact()
#         man.compute_absorption()

#     folder = os.path.join(base_folder,'d%i'%i)
#     check_folder(folder) # create directory if it does not exist
#     man.save(man_name,folder=folder)

#---------------------------------------------------
# Load evaluations and compute bounds
I_list = []; A_list = []
I_min = []; A_min = []; I_max = []; A_max = []
for i,lean_angle in enumerate(lean_angles):

    folder = os.path.join(base_folder,'d%i'%i)
    man.load(man_name,folder=folder)

    # drop NaNs (infeasible input specifications)
    notnan = ~np.isnan(man.absorption_matrix.values).any(axis=(0,1))
    I = man.impact_matrix.values[:,:,notnan]
    A = man.absorption_matrix.values[:,:,notnan]

    I_list += [I]
    A_list += [A]

    I_min += [I.min()]; A_min += [A.min()]
    I_max += [I.max()]; A_max += [A.max()]

# # Automatic bounds
# lb = np.array([min(I_min), min(A_min)])
# ub = np.array([max(I_max), max(A_max)])

lb_n = np.array([0,0])
ub_n = np.array([1,1])

#---------------------------------------------------
# Compute distances

# create empty figure
fig, ax = plt.subplots(figsize=(7, 8))
ax.set_xlabel('Impact on performance')
ax.set_ylabel('Change absoption capability')
ax.set_xlim([lb_n[0],ub_n[0]])
ax.set_ylim([lb_n[1],ub_n[1]])

X = np.empty((0,len(man.margin_nodes)))
Y = np.empty((0,len(man.margin_nodes)))
D = np.empty((0,len(man.design_params)))
distances_mean = np.empty(0)
distances_mean_n = np.empty(0)

distances_all = np.empty((0,n_epochs))
distances_all_n = np.empty((0,n_epochs))

for i,(I,A,lean_angle) in enumerate(zip(I_list,A_list,lean_angles)):

    # Extract x and y
    x = np.mean(I,axis=(1, 2)).reshape(-1,len(man.margin_nodes))
    y = np.mean(A,axis=(1, 2)).reshape(-1,len(man.margin_nodes))

    x_n = np.mean(scaling(I,lb[0],ub[0],1),axis=(1, 2)).reshape(-1,len(man.margin_nodes))
    y_n = np.mean(scaling(A,lb[1],ub[1],1),axis=(1, 2)).reshape(-1,len(man.margin_nodes))

    dist = 0; dist_n = 0
    for x_i,y_i,x_i_n,y_i_n in zip(x,y,x_n,y_n):

        for node in range(len(x_i)):
            s = np.array([x_i[node], y_i[node]])
            pn, d = nearest(lb, ub, s)
            s_n = np.array([x_i_n[node], y_i_n[node]])
            pn_n, d_n = nearest(lb_n, ub_n, s_n)
            dist += d
            dist_n += d_n

    # Loop over each sample
    x_sample = np.mean(I,axis=(1)).reshape(-1,len(man.margin_nodes))
    y_sample = np.mean(A,axis=(1)).reshape(-1,len(man.margin_nodes))
    x_sample_n = np.mean(scaling(I,lb[0],ub[0],1),axis=(1)).reshape(-1,len(man.margin_nodes))
    y_sample_n = np.mean(scaling(A,lb[1],ub[1],1),axis=(1)).reshape(-1,len(man.margin_nodes))

    dist_sample = np.empty(n_epochs)
    dist_sample[:] = np.nan
    dist_sample_n = np.empty(n_epochs)
    dist_sample_n[:] = np.nan
    for i_sample,(x_i,y_i,x_i_n,y_i_n) in enumerate(zip(x_sample,y_sample,x_sample_n,y_sample_n)):
        dist_sample[i_sample] = 0
        dist_sample_n[i_sample] = 0
        for node in range(len(x_i)):
            s = np.array([x_i[node], y_i[node]])
            pn, d = nearest(lb, ub, s)
            s_n = np.array([x_i_n[node], y_i_n[node]])
            pn_n, d_n = nearest(lb_n, ub_n, s_n)
            dist_sample[i_sample] += d
            dist_sample_n[i_sample] += d_n


    X = np.vstack((X,x_n))
    Y = np.vstack((Y,y_n))
    D = np.vstack((D,np.array([15.0,lean_angle])))

    distances_mean = np.append(distances_mean,dist)
    distances_mean_n = np.append(distances_mean_n,dist_n)
    distances_all = np.vstack((distances_all,dist_sample))
    distances_all_n = np.vstack((distances_all_n,dist_sample_n))

    # plot the results
    color = np.random.random((1,3)) 
    ax.scatter(x_n,y_n,c=color)

plt.show()

conf_interval = 20

lower = np.nanmin(distances_all_n,axis=1)
lower_p = np.nanpercentile(distances_all_n, axis=1, q=50-conf_interval/2)
lower_s = np.nanmean(distances_all_n,axis=1) - conf_interval*np.nanstd(distances_all_n,axis=1,ddof=1)

upper = np.nanmax(distances_all_n,axis=1)
upper_p = np.nanpercentile(distances_all_n, axis=1, q=50+conf_interval/2)
upper_s = np.nanmean(distances_all_n,axis=1) + conf_interval*np.nanstd(distances_all_n,axis=1,ddof=1)

mean_d = np.nanmean(distances_all_n,axis=1)
median_d = np.nanmedian(distances_all_n,axis=1)

optimum_mean = np.array([lean_angles[np.argmax(mean_d)], np.max(mean_d)])
optimum_median = np.array([lean_angles[np.argmax(median_d)], np.max(median_d)])
optimum_unnormalized = np.array([lean_angles[np.argmax(distances_mean)], np.max(distances_mean)])
#---------------------------------------------------
# Plot sweep
palette = ['#1C758A', '#CF5044', '#BBBBBB', '#444444']
color = palette[0]

fig, ax = plt.subplots(figsize=(7, 6))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel(r'lean angle ${\theta}$ ($^o$)')
ax.set_ylabel('Distance')

a1 = ax.plot(lean_angles, median_d, color=color, linewidth = 1.0)
# a2 = ax.fill_between(lean_angles, lower_p, upper_p, color=color, alpha=0.2)
a3 = ax.axvline(optimum_median[0], color='k', linestyle='--', linewidth=0.5)
# legend_handles = [(a1[0],a2),a3]
legend_handles = [a1[0],a3]

legned_labels = ['distance from neutral', r'$\theta^{*}$']
ax.legend(legend_handles, legned_labels, fontsize = 8)
# Save plot
fig.savefig(os.path.join(img_folder,'optimization_trace.pdf'),
            format='pdf', dpi=200, bbox_inches='tight')

plt.show()