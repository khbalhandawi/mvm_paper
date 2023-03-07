import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

from mvm import nearest
from mvm.DOELib import scaling
from man_defs import get_man_combined_poly, get_man_combined_circ, get_man

###########################################################
#                    ANALYTICAL PROBLEM                   #
###########################################################
# # get man object and load it
# folder = os.path.join('data','strut','C1'); lean = 0.0; height = 15.0
# # folder = os.path.join('data','strut','C2'); lean = 30.0; height = 15.0
# # folder = os.path.join('data','strut','C3'); lean = 0.0; height = 20.0
# # folder = os.path.join('data','strut','C4'); lean = 30.0; height = 20.0

# # off-the shelf parts
# widths = list(range(60,120+10,10))
# # create material dictionary
# n_materials = 15

# material_dict_user = {
#     'start'   : {
#         'sigma_y' : 92, # MPa
#         'rho' : 11.95e-06, # kg/mm3
#         'cost' : 0.1 # USD/kg
#         },
#     'end'  : {
#         'sigma_y' : 828, # MPa
#         'rho' : 4.43e-06, # kg/mm3
#         'cost' : 1.10 # USD/kg
#         },
# }
# # generate 10 materials by linearly interpolating
# df = pd.DataFrame(columns=material_dict_user['start'].keys(), index=range(n_materials), dtype=float)
# df.iloc[0] = material_dict_user['start']
# df.iloc[-1] = material_dict_user['end']
# df.interpolate(method='linear',axis=0,inplace=True)
# material_dict = df.transpose().to_dict()

# man = get_man(height=height,lean=lean,materials=material_dict,widths=widths,
#     train_surrogate=False,man_folder=folder,overwrite=False,name='strut_s',num_threads=1)

# man_name = 'strut_s'
# man.load(man_name,folder=folder)

###########################################################
#                       FEA PROBLEM                       #
###########################################################

base_folder = os.path.join('data','strut_fea_poly') # choose a concept (polygonal or circumferential)
get_man_combined = get_man_combined_poly

# base_folder = os.path.join('data','strut_fea_circ') # choose a concept (polygonal or circumferential)
# get_man_combined = get_man_combined_circ

# get man object for the fea problem and load it
folder = os.path.join(base_folder,'C1'); lean = 0.0; height = 15.0
# folder = os.path.join(base_folder,'C2'); lean = 30.0; height = 15.0
# folder = os.path.join(base_folder,'C3'); lean = 0.0; height = 17.0
# folder = os.path.join(base_folder,'C4'); lean = 30.0; height = 17.0

# off-the shelf parts
widths = list(range(60,120+10,10))
# create material dictionary
n_materials = 15

material_dict_user = {
    'start'   : {
        'sigma_y' : 92, # MPa
        'rho' : 11.95e-06, # kg/mm3
        'cost' : 0.1 # USD/kg
        },
    'end'  : {
        'sigma_y' : 828, # MPa
        'rho' : 4.43e-06, # kg/mm3
        'cost' : 1.10 # USD/kg
        },
}
# generate 10 materials by linearly interpolating
df = pd.DataFrame(columns=material_dict_user['start'].keys(), index=range(n_materials), dtype=float)
df.iloc[0] = material_dict_user['start']
df.iloc[-1] = material_dict_user['end']
df.interpolate(method='linear',axis=0,inplace=True)
material_dict = df.transpose().to_dict()

man = get_man_combined(height=height,lean=lean,materials=material_dict,widths=widths,
    train_surrogate=False,man_folder=folder,overwrite=False,name='strut_comb',num_threads=1)

man_name = 'strut_comb'
man.load(man_name,folder=folder)

###########################################################
# drop NaNs (infeasible input specifications)
notnan = ~np.isnan(man.absorption_matrix.values).any(axis=(0,1))
I = man.impact_matrix.values[:,:,notnan]
I[I==0] = np.nan # remove zero impact values from averaging
for k in range(I.shape[2]): # return rows that are full of zeros
    I[np.all(np.isnan(I[:,:,k]),axis=1),:,k] = 0
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
    col[np.isnan(col)] = 0
    all_matrix = np.hstack((all_matrix,col))

for i in range(n_spec):
    col = A[:,i,:].reshape(-1,1)
    all_matrix = np.hstack((all_matrix,col))

df_matrix = pd.DataFrame(all_matrix, columns=columns)
df_matrix = df_matrix.astype({'node':'int'})

# sns.pairplot(df_matrix, hue="node")
df_matrix.to_csv(os.path.join(folder,'df_matrix.csv'))

###########################################################
# create a dataframe of average impact and absorption

# Extract x and y
x = np.nanmean(I,axis=1).ravel().reshape(-1,1) # average along performance parameters (assumes equal weighting)
y = np.nanmean(A,axis=1).ravel().reshape(-1,1) # average along input specs (assumes equal weighting)

fig1, ax1 = plt.subplots(figsize=(7, 8))
ax1.set_xlabel('Impact on performance')
ax1.set_ylabel('Frequency')

fig2, ax2 = plt.subplots(figsize=(7, 8))
ax2.set_xlabel('Change absorption capability')
ax2.set_ylabel('Frequency')

for i in range(len(man.margin_nodes)):
    x_node = np.nanmean(I,axis=1)[i,:]
    y_node = np.nanmean(A,axis=1)[i,:]
    ax1.hist(x_node,label="%i" %(i+1))
    ax2.hist(y_node,label="%i" %(i+1))

ax1.legend()
ax2.legend()
plt.show()

mean_matrix = np.empty((0,1))
for i in range(n_nodes):
    mean_matrix = np.vstack((mean_matrix,(i+1)*np.ones((n_samples,1))))

mean_matrix = np.hstack((mean_matrix,x))
mean_matrix = np.hstack((mean_matrix,y))

columns = ['node','Impact','Absorption']
df_mean = pd.DataFrame(mean_matrix, columns=columns)
df_mean = df_mean.astype({'node':'int'})

# sns.jointplot(data=df_mean, x="Impact", 
    # y="Absorption", hue="node")

df_mean.to_csv(os.path.join(folder,'df_mean.csv'))

# Total matrix generation
columns = []
total_matrix = np.empty((n_samples,0))
excess_matrix = np.empty((1000,0))

for i in range(n_spec):
    columns += ["S%i" %(i+1)]
    col = man.input_specs[i].values[notnan].reshape(-1,1)
    total_matrix = np.hstack((total_matrix,col))

for i in range(n_nodes):
    columns += ["D%i" %(i+1)]
    col = man.decisions[i].selection_values.values[1:][notnan].reshape(-1,1)
    total_matrix = np.hstack((total_matrix,col))

for i in range(n_nodes):
    columns += ["E%i" %(i+1)]
    col = man.margin_nodes[i].excess.values[1:][notnan].reshape(-1,1)
    total_matrix = np.hstack((total_matrix,col))
    col = man.margin_nodes[i].excess.values[1:].reshape(-1,1)
    excess_matrix = np.hstack((excess_matrix,col))

for i in range(n_nodes):
    columns += ["I%i" %(i+1)]
    col = np.nanmean(I,axis=1)[i].reshape(-1,1)
    total_matrix = np.hstack((total_matrix,col))

for i in range(n_nodes):
    columns += ["A%i" %(i+1)]
    col = np.nanmean(A,axis=1)[i].reshape(-1,1)
    total_matrix = np.hstack((total_matrix,col))

df_total = pd.DataFrame(total_matrix, columns=columns)
df_total.to_csv(os.path.join(folder,'df_total.csv'))
reliability = len(excess_matrix[np.all(excess_matrix > 0, axis=1),:])/1000
print("reliability is: %.3f" %reliability)
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
x = np.nanmean(I,axis=(1, 2)).reshape(-1,len(man.margin_nodes))
y = np.nanmean(A,axis=(1, 2)).reshape(-1,len(man.margin_nodes))

###########################################################
# get bounds
lb = np.array([-0.009490439, 1.143639527])
ub = np.array([0.3172583, 4.1707973])

# lb = np.array([-0.0639315,  0.       ])
# ub = np.array([ 0.33307271, 13.32143689])

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
