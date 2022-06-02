import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

from mvm import nearest
from mvm.DOELib import scaling
from man_defs import get_man_combined, get_man

# get man object for the fea problem and load it
folder = os.path.join('data','strut_fea_50','C1'); lean = 0.0
# folder = os.path.join('data','strut_fea_50','C2'); lean = 20.408163265306122
# folder = os.path.join('data','strut_fea_50','C3'); lean = 28.571428571428573

man = get_man_combined()
man_name = 'strut_comb'
man.load(man_name,folder=folder)

# # get man object and load it
# # folder = os.path.join('data','strut','C1'); lean = 0.0
# # folder = os.path.join('data','strut','C2'); lean = 10.0
# folder = os.path.join('data','strut','C3'); lean = 30.0

# man = get_man_combined()
# man_name = 'strut_s'
# man.load(man_name,folder=folder)

###########################################################
# drop NaNs (infeasible input specifications)
notnan = ~np.isnan(man.absorption_matrix.values).any(axis=(0,1))
I = man.impact_matrix.values[:,:,notnan]
A = man.absorption_matrix.values[:,:,notnan]

###########################################################
# create a dataframe of impact and absorptions
columns = ['node',] + [p.key for p in man.performances] + [s.key for s in man.input_specs]

n_nodes = len(man.margin_nodes)
n_spec = len(man.input_specs)
n_perf = len(man.performances)
n_samples = A.shape[2]

all_matrix = np.empty((0,1))
for i in range(n_nodes):
    all_matrix = np.vstack((all_matrix,(i+1)*np.ones((n_samples,1))))

for i in range(n_perf):
    col = I[:,i,:].reshape(-1,1)
    all_matrix = np.hstack((all_matrix,col))

for i in range(n_spec):
    col = A[:,i,:].reshape(-1,1)
    all_matrix = np.hstack((all_matrix,col))

df_matrix = pd.DataFrame(all_matrix, columns=columns)
df_matrix = df_matrix.astype({'node':'int'})

sns.pairplot(df_matrix, hue="node")
df_matrix.to_csv(os.path.join(folder,'df_matrix.csv'))

###########################################################
# create a dataframe of average impact and absorption
columns = ['node','Impact','Absorption']

# Extract x and y
x = np.nanmean(I,axis=1).ravel().reshape(-1,1) # average along performance parameters (assumes equal weighting)
y = np.nanmean(A,axis=1).ravel().reshape(-1,1) # average along input specs (assumes equal weighting)

mean_matrix = np.empty((0,1))
for i in range(n_nodes):
    mean_matrix = np.vstack((mean_matrix,(i+1)*np.ones((n_samples,1))))

mean_matrix = np.hstack((mean_matrix,x))
mean_matrix = np.hstack((mean_matrix,y))

df_mean = pd.DataFrame(mean_matrix, columns=columns)
df_mean = df_mean.astype({'node':'int'})

sns.jointplot(data=df_mean, x="Impact", 
    y="Absorption", hue="node")

df_mean.to_csv(os.path.join(folder,'df_mean.csv'))


###########################################################
# Scatter plot of average absorption and impact
# create empty figure
# fig, ax = plt.subplots(figsize=(7, 8))
# ax.set_xlabel('Impact on performance')
# ax.set_ylabel('Change absorption capability')

# # Extract x and y
# x = np.mean(I,axis=(1, 2)).ravel()  # average along performance parameters (assumes equal weighting)
# y = np.mean(A,axis=(1, 2)).ravel()  # average along input specs (assumes equal weighting)

# # plot the results
# ax.scatter(x, y)
# plt.show()

###########################################################
# calculate the distance metric to evaluate different designs
# Calculate distance metric

# Extract x and y
x = np.mean(I,axis=(1, 2)).reshape(-1,len(man.margin_nodes))
y = np.mean(A,axis=(1, 2)).reshape(-1,len(man.margin_nodes))

###########################################################
# get bounds
lb = np.array([-0.009490439, 1.143639527])
ub = np.array([0.3172583, 4.1707973])

lb = np.array([-0.0639315,  0.       ])
ub = np.array([ 0.33307271, 13.32143689])

lb_n = np.array([0, 0])
ub_n = np.array([1, 1])
x_n = scaling(x,lb[0],ub[0],1)
y_n = scaling(y,lb[1],ub[1],1)

# lb = np.array([I.min(), A.min()])
# ub = np.array([I.max(), A.max()])

dist = 0
dist_n = 0
for x_i,y_i in zip(x,y):

    for node in range(len(x_i)):
        s = np.array([x_i[node], y_i[node]])
        pn, d = nearest(lb, ub, s)
        s = scaling(s,lb,ub,1)
        pn_n, d_n = nearest(lb_n, ub_n, s)
        dist += d
        dist_n += d_n

# create empty figure
colors = ['#FF0000', '#27BE1E', '#0000FF']
fig, ax = plt.subplots(figsize=(7, 8))
ax.set_xlabel('Impact on performance')
ax.set_ylabel('Change absorption capability')

ax.scatter(x_n,y_n, c=colors[0])

for x_i,y_i in zip(x_n,y_n):
    for node in range(len(x_i)):
        s = np.array([x_i[node], y_i[node]])
        pn, _ = nearest(lb_n, ub_n, s)

        x_d = np.array([s[0], pn[0]])
        y_d = np.array([s[1], pn[1]])

        ax.plot(x_d, y_d, marker='.', linestyle='--', color=colors[0])

# ax.axis('equal')

axes = ax.axis()
ax.plot([lb_n[0], ub_n[0]], [lb_n[1], ub_n[1]], color='k', linestyle=(5, (10, 5)))

plt.show()