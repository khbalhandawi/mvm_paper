import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os, glob, shutil
import itertools
from copy import deepcopy
import multiprocess as mp
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot, iplot, init_notebook_mode # Import offline plot functions from plotly
from scipy import stats
import random

from mvm import nearest
from mvm.DOELib import scaling
from mvm.utilities import check_folder, parallel_sampling

from mvm import MarginNetwork
from typing import List, Union
from man_defs import get_man_combined_poly, get_man_combined_circ

mpl.rcParams['axes.linewidth'] = 0.5

# must define parallelizable functions outside main
def evaluate_design_manual(i: int, design: List[Union[int,float]], mans: List[MarginNetwork], n_epochs: int,
                    base_folder: str, man_name: str, strategies: List[str], randomize: bool, process_ids: List[int]=None):

    # Select a man to forward based on process id
    if process_ids == None:
        pid = None
        man = mans[0]
    else:
        pid = mp.current_process()._identity[0] - process_ids[0]
        man = mans[pid]

    folder = os.path.join(base_folder,'d%i'%i)
    exists = check_folder(folder) # create directory if it does not exist

    if not exists:
        man.nominal_design_vector = np.array(design[:len(man.design_params)])
        man.reset()
        man.reset_inputs('all')
        man.reset_outputs()

        # Perform Monte-Carlo simulation
        for n in range(n_epochs):
            if randomize:
                man.randomize()
            man.init_decisions()
            man.decision_vector = design[len(man.design_params):]
            man.allocate_margins(strategies)
            man.forward()
            man.compute_impact()
            man.compute_absorption(max_iter=100)

        man.save(man_name,folder=folder,results_only=True)

if __name__ == "__main__":

    resume = True
    random.seed(10)

    # case 1:
    base_dir = "strut_fea_poly"
    get_man_combined = get_man_combined_poly
    strategies = ['manual',]*3
    strategy_name = 'manual'
    sampling = 'stochastic'
    n_epochs = 50
    d_labels = ['height','lean angle','width','material','n_struts']
    perf_labels = ['weight','cost']

    # # case 2:
    # base_dir = "strut_fea_circ"
    # get_man_combined = get_man_combined_circ
    # strategies = ['manual','manual','manual','min_excess']
    # strategy_name = 'manual'
    # sampling = 'stochastic'
    # n_epochs = 100
    # d_labels = ['height','lean angle','width','material','n_struts','t_shroud']
    # perf_labels = ['weight','cost','W_shroud','c_shroud']

    # # case 3:
    # base_dir = "strut_fea_poly"
    # get_man_combined = get_man_combined_poly
    # strategies = ['manual',]*3
    # strategy_name = 'manual'
    # sampling = 'deterministic'
    # n_epochs = 1
    # d_labels = ['height','lean angle','width','material','n_struts']
    # perf_labels = ['weight','cost']

    # # case 4:
    # base_dir = "strut_fea_circ"
    # get_man_combined = get_man_combined_circ
    # strategies = ['manual','manual','manual','min_excess']
    # strategy_name = 'manual'
    # sampling = 'deterministic'
    # n_epochs = 1
    # d_labels = ['height','lean angle','width','material','n_struts','t_shroud']
    # perf_labels = ['weight','cost','W_shroud','c_shroud']

    # # case 5:
    # base_dir = "strut_fea_poly"
    # get_man_combined = get_man_combined_poly
    # strategies = ['min_excess',]*3
    # strategy_name = 'minexcess'
    # sampling = 'deterministic'
    # n_epochs = 1
    # d_labels = ['height','lean angle','width','material','n_struts']
    # perf_labels = ['weight','cost']

    # # case 6:
    # base_dir = "strut_fea_circ"
    # get_man_combined = get_man_combined_circ
    # strategies = ['min_excess',]*4
    # strategy_name = 'minexcess'
    # sampling = 'deterministic'
    # n_epochs = 1
    # d_labels = ['height','lean angle','width','material','n_struts','t_shroud']
    # perf_labels = ['weight','cost','W_shroud','c_shroud']

    # get the man object for the fea problem and load it
    base_folder = os.path.join('data',base_dir,'opt_%s_%s' %(strategy_name,sampling))
    img_folder = os.path.join('images',base_dir,'opt_%s_%s' %(strategy_name,sampling))

    check_folder(img_folder)

    # off-the shelf parts
    widths = list(range(60,120+10,10))
    material_dict = {
        'Steel'   : {
            'sigma_y' : 250, # MPa
            'rho' : 10.34e-06, # kg/mm3
            'cost' : 0.09478261, # USD/kg
            },
        'Inconel'  : {
            'sigma_y' : 460,  # MPa
            'rho' : 8.19e-06,  # kg/mm3
            'cost' : 0.46,  # USD/kg
            },
        'Titanium'  : {
            'sigma_y' : 828, # MPa
            'rho' : 4.43e-06, # kg/mm3
            'cost' : 1.10 # USD/kg
            },
    }
    height = 15.0
    lean = 0.0

    man_name = 'strut_comb'
    man = get_man_combined(height=height,lean=lean,materials=material_dict,widths=widths,
        train_surrogate=False,man_folder=base_folder,overwrite=True,name=man_name,num_threads=1)

    nominal_specs = [434, 432, 80, 80]
    man.input_specs[0].value = nominal_specs[0]
    man.input_specs[1].value = nominal_specs[1]
    man.input_specs[2].value = nominal_specs[2]
    man.input_specs[3].value = nominal_specs[3]

    man.input_specs[2].inc_user = 5
    man.input_specs[3].inc_user = 5

    # for fea problem
    lb = np.array([-0.5973142383857094, 0.0])
    ub = np.array([3.0949267424315656, 7.4235687891259134])

    lb_n = np.array([0.0, 0.0])
    ub_n = np.array([1.0, 1.0])
    ############################################################
    # Optimization sweep
    universe_d = []
    # h
    lb_h = 14.0 # man.design_params[0].universe[0]
    ub_h = 20.0 # man.design_params[0].universe[1]
    universe_d += [list(np.arange(lb_h,ub_h+2,2))]

    # theta
    lb_theta = man.design_params[1].universe[0]
    ub_theta = man.design_params[1].universe[1]
    universe_d += [list(np.arange(lb_theta,ub_theta+10,10))]
    
    # Generate full-factorial DOE
    universe_decision = [d if s == 'manual' else [None,] for d,s in zip(man.universe_decision,strategies)]
    universe = universe_d + universe_decision
    design_doe = list(itertools.product(*universe))
    # design_doe = [(14, 10, 90, 'Inconel', 10),] # try a unique design
    n_designs = len(design_doe)
    print("Analyzing %i concepts" %n_designs)
    #---------------------------------------------------
    # Evaluate all the different designs
    # Parallel computation if num_threads > 1
    num_threads = 4

    man_objs = []
    for pid in range(num_threads):
        man_objs += [deepcopy(man)]

    kwargs = {'n_epochs': n_epochs, 
              'base_folder': base_folder, 
              'man_name': man_name,
              'strategies': strategies,
              'randomize': True if sampling=='stochastic' else False}

    vargs_iterator = [[i,design,] for i,design in enumerate(design_doe)]
    # Partition into chunks
    # vargs_iterator = vargs_iterator[:2000] #1
    # vargs_iterator = vargs_iterator[2000:4000] #2
    # vargs_iterator = vargs_iterator[4000:6000] #3
    # vargs_iterator = vargs_iterator[6000::] #4

    vkwargs_iterator = [{},] * len(vargs_iterator)
    fargs = []
    fkwargs = kwargs
    fkwargs['mans'] = man_objs

    if not resume:
        for item in glob.glob(base_folder+"/d*"):
            if os.path.isdir(item):
                shutil.rmtree(item)
    results = parallel_sampling(evaluate_design_manual,vargs_iterator,vkwargs_iterator,fargs,fkwargs,num_threads=num_threads) # For manual case, uncomment
    # sys.exit(0)
    
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
    for id,design in enumerate(design_doe):
        excess_vector = np.empty((0,n_samples))

        sys.stdout.write("Progress: %d%%   \r" % ((id / n_designs) * 100))
        sys.stdout.flush()

        dist = 0.0; dist_n = 0.0
        for n,node in enumerate(man.margin_nodes):
            
            multi_index += [(id,n)] # for nodal dataframe

            folder = os.path.join(base_folder,'d%i'%id)
            man.load(man_name,folder=base_folder) # load the basic MAN
            man.load(man_name,folder=folder,results_only=True) # load the sampled data

            impact_values = man.impact_matrix.values
            impact_values[impact_values==0] = np.nan # remove zero impact values from averaging
            for k in range(impact_values.shape[2]): # return rows that are full of zeros
                impact_values[np.all(np.isnan(impact_values[:,:,k]),axis=1),:,k] = 0

            # calculate means
            I = impact_values[n,:,:]
            A = man.absorption_matrix.values[n,:,:]

            i = np.nanmean(I,axis=1) # absorption
            a = np.nanmean(A,axis=1) # impact
            # scaled
            i_n = scaling(i,lb[0],ub[0],1) # absorption
            a_n = scaling(a,lb[1],ub[1],1) # impact

            point = np.array([np.nanmean(i), np.mean(a)]) # nan mean to ignore nans (only for I!!)
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
            for decision_i,(decision,strategy) in enumerate(zip(man.decisions,strategies)):
                if strategy == 'manual':
                    nodal_dict[decision.key] = design[d_i+1+decision_i]
                elif strategy == 'min_excess':
                    nodal_dict[decision.key] = stats.mode(decision.selection_values.values,keepdims=False)[0] # get most frequent decision
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
            excess_vector = np.vstack((excess_vector,excess_normalized.reshape(-1,n_samples)))

        # excess
        excess_mean = np.mean(excess_vector,axis=1)

        # performance
        perf_vector = np.empty((0,n_samples))
        for p,performance in enumerate(man.performances):
            perf_vector = np.vstack((perf_vector,performance.values.reshape(-1,n_samples))) # performances
        perf_mean = np.mean(perf_vector,axis=1)

        # reliability
        excess_feasible = excess_vector[:,(excess_vector >= 0.0).all(axis=0)]
        reliability = excess_feasible.shape[1] / excess_vector.shape[1]

        # mean impact and absorption
        mean_i = np.nanmean(impact_values,axis=(1,2)) # absorption
        mean_a = np.nanmean(man.absorption_matrix.values,axis=(1,2)) # impact

        # construct dict
        total_dict = {'id':id}
        for d_i,d in enumerate(man.design_params):
            total_dict[d.key] = design[d_i]
        for decision_i,(decision,strategy) in enumerate(zip(man.decisions,strategies)):
            if strategy == 'manual':
                total_dict[decision.key] = design[d_i+1+decision_i]
            elif strategy == 'min_excess':
                total_dict[decision.key] = stats.mode(decision.selection_values.values,keepdims=False)[0] # get most frequent decision
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
    df_bak = df
    df_nodal_bak = df_nodal

    labels = d_labels + \
        [e.key for e in man.margin_nodes] + perf_labels + \
        list(df.columns[len(man.design_params)+len(man.decisions)+len(man.margin_nodes)+len(man.performances):])

    labels_nodal = d_labels + \
        list(df_nodal.columns[len(man.design_params)+len(man.decisions):])

    df = df.set_axis(labels, axis=1, inplace=False)
    df_nodal = df_nodal.set_axis(labels_nodal, axis=1, inplace=False)

    #---------------------------------------------------
    # Create PCP

    col_name = 'material'
    new_col_name = 'dummy'
    def replace_string_col(df,col_name,new_col_name):
        # create dummy column for material 
        # https://stackoverflow.com/a/64146570
        group_vars = df[col_name].unique()
        dfg = pd.DataFrame({col_name:df[col_name].unique()})
        dfg[new_col_name] = dfg.index
        df = pd.merge(df, dfg, on = col_name, how='left')
        return df,dfg

    # PCP for DOE
    df,dfg = replace_string_col(df,col_name,new_col_name)
    dimensions = []
    for label,column in zip(labels,df.columns):
        if label == 'material':
            dimensions += [dict(range=[0,df['dummy'].max()],
                            tickvals = dfg['dummy'], ticktext = dfg['material'],
                            label='material', values=df['dummy']),]
        else:
            dimensions += [{
                'range': [df[column].dropna().min(),df[column].dropna().max()],
                'label': label,
                'values': df[column]
            },]

    color_col = 'IoP E3'
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
    plot(data_pd, show_link=False, filename = os.path.join(img_folder,'pcp_doe.html'), auto_open=True)

    #-------------------------------------------------
    # PCP for nodal data
    df_nodal,dfg_nodal = replace_string_col(df_nodal,col_name,new_col_name)
    dimensions_nodal = []
    for label,column in zip(labels_nodal,df_nodal.columns):
        if label == 'material':
            dimensions_nodal += [dict(range=[0,df_nodal['dummy'].max()],
                            tickvals = dfg_nodal['dummy'], ticktext = dfg_nodal['material'],
                            label='material', values=df_nodal['dummy']),]
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
    df_export_real = df.dropna().drop(['dummy'],axis=1) # export without NaNs
    df_export = df.drop(['material'],axis=1)
    df_export = df_export.rename(columns={'dummy':'material'})

    # shift column 'Name' to first position
    material = df_export.pop('material')
    df_export.insert(3, 'material', material)
    
    ## nodal
    df_export_nodal_real = df_nodal.dropna().drop(['dummy'],axis=1) # export without NaNs
    df_export_nodal = df_nodal.drop(['material'],axis=1)
    df_export_nodal = df_export_nodal.rename(columns={'dummy':'material'})

    # shift column 'Name' to first position
    material = df_export_nodal.pop('material')
    df_export_nodal.insert(3, 'material', material)

    # export to csv
    df_export_real.to_csv(os.path.join(base_folder,'doe_data_real.csv'))
    df_export_nodal_real.to_csv(os.path.join(base_folder,'nodal_data_real.csv'))
    df_export.to_csv(os.path.join(base_folder,'doe_data.csv'))
    df_export_nodal.to_csv(os.path.join(base_folder,'nodal_data.csv'))