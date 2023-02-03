import sys
import matplotlib as mpl
import numpy as np
import os
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot, iplot, init_notebook_mode # Import offline plot functions from plotly
from scipy import stats

from mvm import nearest
from mvm.DOELib import scaling
from mvm.utilities import check_folder

from mvm import MarginNetwork
from typing import List, Union
from man_defs import get_man_combined

mpl.rcParams['axes.linewidth'] = 0.5

if __name__ == "__main__":

    # get the man object for the fea problem and load it
    # base_folder = os.path.join('data','strut_fea','opt_minexcess_deterministic')
    # img_folder = os.path.join('images','strut_fea','opt_minexcess_deterministic')
    # n_epochs = 1
    base_folder = os.path.join('data','strut_fea')
    img_folder = os.path.join('images','strut_fea')
    n_epochs = 1000

    check_folder(img_folder)

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
    height = 15.0
    lean = 0.0

    man_name = 'strut_comb'
    man = get_man_combined(height=height,lean=lean,materials=material_dict,widths=widths,
        train_surrogate=False,man_folder=base_folder,overwrite=True,name=man_name,num_threads=1)

    # for fea problem
    lb = np.array([-0.009630757, 1.152203562])
    ub = np.array([0.3329746, 4.3151886])

    lb_n = np.array([0.0, 0.0])
    ub_n = np.array([1.0, 1.0])

    designs = [[15.0,0.0],[15.0,30.0],[17.0,0.0],[17.0,30.0]]
    n_designs = len(designs)

    #---------------------------------------------------
    # Load evaluations

    columns = ['id',] + [d.key for d in man.design_params] + [decision.key for decision in man.decisions] + \
        [e.key for e in man.margin_nodes] + [p.key for p in man.performances] + ['reliability'] + \
        ['IoP %s'%e.key for e in man.margin_nodes] + ['CAC %s'%e.key for e in man.margin_nodes] + ['IoP','CAC',]

    columns_nodal = ['id',] + [d.key for d in man.design_params] + [decision.key for decision in man.decisions] + \
        ['node', 'excess', 'threshold', 'decided_value'] + ['IoP %s'%p.key for p in man.performances] + \
        ['CAC %s'%s.key for s in man.input_specs]

    df = pd.DataFrame(columns=columns)
    df_nodal = pd.DataFrame(columns=columns_nodal)

    n_nodes = len(man.margin_nodes)
    n_spec = len(man.input_specs)
    n_perf = len(man.performances)
    n_samples = n_epochs

    multi_index = []
    for id,design in enumerate(designs):
        excess_vectors = []

        sys.stdout.write("Progress: %d%%   \r" % ((id / n_designs) * 100))
        sys.stdout.flush()
        
        man.nominal_design_vector = np.array(design)

        dist = 0.0; dist_n = 0.0
        for n,node in enumerate(man.margin_nodes):
            
            multi_index += [(id,n)] # for nodal dataframe

            folder = os.path.join(base_folder,'C%i'%(id+1))
            man.load(man_name,folder=folder,results_only=False) # load the sampled data

            # calculate means
            notnan = ~np.isnan(man.absorption_matrix.values).any(axis=(0,1))
            I = man.impact_matrix.values[n,:,notnan]
            A = man.absorption_matrix.values[n,:,notnan]

            i = np.mean(I,axis=1) # absorption
            a = np.nanmean(A,axis=1) # impact
            # scaled
            i_n = scaling(i,lb[0],ub[0],1) # absorption
            a_n = scaling(a,lb[1],ub[1],1) # impact

            point = np.array([np.mean(i), np.mean(a)]) # nan mean to ignore nans
            pn, dist_node = nearest(lb, ub, point)
            dist += dist_node
            # point_n = np.array([np.mean(i_n), np.mean(a_n)])
            # pn_n, dist_node_n = nearest(lb_n, ub_n, point_n)
            # dist_n += dist_node_n

            e = np.mean(node.excess.values/node.target.values) # excess
            tt = np.mean(node.target.values) # excess
            dv = np.mean(node.decided_value.values) # excess

            # construct nodal dict
            nodal_dict = {'id':id}
            for d_i,d in enumerate(man.design_params):
                nodal_dict[d.key] = design[d_i]
            for decision_i,decision in enumerate(man.decisions):
                nodal_dict[decision.key] = stats.mode(decision.selection_values.values,keepdims=True)[0][0] # get most frequent decision
            nodal_dict['node'] = n
            nodal_dict['excess'] = e
            nodal_dict['threshold'] = tt
            nodal_dict['decided_value'] = dv
            for p_i,p in enumerate(man.performances):
                nodal_dict['IoP %s'%p.key] = i[p_i]
            for s_i,s in enumerate(man.input_specs):
                nodal_dict['CAC %s'%s.key] = a[s_i]
            nodal_dict['value'] = dist_node
            df_nodal = df_nodal.append(nodal_dict,ignore_index=True)

            excess_normalized = node.excess.values/node.target.values
            excess_vectors += [excess_normalized[1:][notnan]]

        excess_vector = np.array(excess_vectors)
        # excess
        excess_mean = np.mean(excess_vector,axis=1)

        # performance
        perf_vector = np.empty((0,n_samples))
        for p,performance in enumerate(man.performances):
            perf_vector = np.vstack((perf_vector,performance.values[1:].reshape(1,-1))) # performances
        perf_mean = np.mean(perf_vector,axis=1)

        # reliability
        reliability = notnan.sum() / n_samples

        # mean impact and absorption
        mean_i = np.mean(man.impact_matrix.values,axis=(1,2)) # absorption
        mean_a = np.nanmean(man.absorption_matrix.values,axis=(1,2)) # impact

        # construct dict
        total_dict = {'id':id}
        for d_i,d in enumerate(man.design_params):
            total_dict[d.key] = design[d_i]
        for decision_i,decision in enumerate(man.decisions):
            total_dict[decision.key] = stats.mode(decision.selection_values.values,keepdims=True)[0][0] # get most frequent decision
        for node_i,node in enumerate(man.margin_nodes):
            total_dict[node.key] = excess_mean[node_i]
        for p_i,p in enumerate(man.performances):
            total_dict[p.key] = perf_mean[p_i]
        total_dict['reliability'] = reliability
        for node_i,node in enumerate(man.margin_nodes):
            total_dict['IoP %s'%node.key] = mean_i[node_i]
            total_dict['CAC %s'%node.key] = mean_a[node_i]
        total_dict['IoP'] = np.mean(mean_i)
        total_dict['CAC'] = np.mean(mean_a)
        total_dict['value'] = dist
        df = df.append(total_dict,ignore_index=True)

    df = df.set_index('id')
    df_nodal = df_nodal.set_index('id')

    labels = ['height','lean angle','width','material','n_struts'] + \
        [e.key for e in man.margin_nodes] + ['weight','cost'] + \
        list(df.columns[len(man.design_params)+len(man.decisions)+len(man.margin_nodes)+len(man.performances):])

    labels_nodal = ['height','lean angle','width','material','n_struts'] + \
        list(df_nodal.columns[len(man.design_params)+len(man.decisions):])

    df = df.set_axis(labels, axis=1, inplace=False)
    df_nodal = df_nodal.set_axis(labels_nodal, axis=1, inplace=False)
    df
    #---------------------------------------------------
    # Create PCP

    lookup = dict()
    for key in material_dict.keys():
        if key in list(range(0,5)):
            lookup[key] = "steel"
        elif key in list(range(5,6)):
            lookup[key] = "Inconel"
        else:
            lookup[key] = "titanium"

    lookup_number = {
        "steel" : 0,
        "Inconel" : 1,
        "titanium" : 2
    }

    # PCP for concepts
    df['material_string'] = df['material'].map(lookup)
    df['material'] = df['material_string'].map(lookup_number)
    dimensions = []
    for label,column in zip(labels,df.columns):
        if label == 'material':
            dimensions += [dict(range=[df['material'].min(),df['material'].max()],
                            tickvals = list(lookup_number.values()), ticktext = list(lookup_number.keys()),
                            label='material', values=df['material']),]
        else:
            dimensions += [{
                'range': [df[column].dropna().min(),df[column].dropna().max()],
                'label': label,
                'values': df[column]
            },]

    color_col = 'value'
    data_pd = [
        go.Parcoords(
            line = dict(color=df[color_col],
                        colorscale = 'tealrose',
                        cmin = df[color_col].dropna().min(),
                        cmax = df[color_col].dropna().max(),
                        showscale = True),
            dimensions = list(dimensions)
        )
    ]
    plot(data_pd, show_link=False, filename = os.path.join(img_folder,'pcp_concepts.html'), auto_open=True)

    #-------------------------------------------------
    # PCP for nodal data
    df_nodal['material_string'] = df_nodal['material'].map(lookup)
    df_nodal['material'] = df_nodal['material_string'].map(lookup_number)
    dimensions_nodal = []
    for label,column in zip(labels_nodal,df_nodal.columns):
        if label == 'material':
            dimensions_nodal += [dict(range=[df_nodal['material'].min(),df_nodal['material'].max()],
                            tickvals = list(lookup_number.values()), ticktext = list(lookup_number.keys()),
                            label='material', values=df_nodal['material']),]
        else:
            dimensions_nodal += [{
                'range': [df_nodal[column].dropna().min(),df_nodal[column].dropna().max()],
                'label': label,
                'values': df_nodal[column]
            },]

    color_col = 'node'
    data_pd_nodal = [
        go.Parcoords(
            line = dict(color=df_nodal[color_col],
                        colorscale = 'tealrose',
                        cmin = df_nodal[color_col].dropna().min(),
                        cmax = df_nodal[color_col].dropna().max(),
                        showscale = True),
            dimensions = list(dimensions_nodal)
        )
    ]
    plot(data_pd_nodal, show_link=False, filename = os.path.join(img_folder,'pcp_nodal.html'), auto_open=True)

    #---------------------------------------------------
    # export data

    ## doe
    df_export = df.dropna().drop(['material'],axis=1) # export without NaNs
    df_export = df_export.rename(columns={'material_string':'material'})

    # shift column 'Name' to first position
    material = df_export.pop('material')
    df_export.insert(3, 'material', material)
    
    ## nodal
    df_export_nodal = df_nodal.dropna().drop(['material'],axis=1) # export without NaNs
    df_export_nodal = df_export_nodal.rename(columns={'material_string':'material'})

    # shift column 'Name' to first position
    material = df_export_nodal.pop('material')
    df_export_nodal.insert(3, 'material', material)

    # export to csv
    df_export.to_csv(os.path.join(base_folder,'concepts_data.csv'))
    df_export_nodal.to_csv(os.path.join(base_folder,'concepts_nodal_data.csv'))