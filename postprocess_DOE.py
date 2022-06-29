import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import itertools
from copy import deepcopy
import multiprocess as mp
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot, iplot, init_notebook_mode # Import offline plot functions from plotly

from mvm import nearest
from mvm.DOELib import scaling
from mvm.utilities import check_folder, parallel_sampling

from mvm import FixedParam, DesignParam, InputSpec, Behaviour, Performance, MarginNode, MarginNetwork, Decision
from mvm import GaussianFunc, UniformFunc
from typing import List

mpl.rcParams['axes.linewidth'] = 0.5

# must define parallelizable functions outside main
def get_man_DOE():
    widths = list(range(60,120+10,10))

    # define fixed parameters
    i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
    i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
    i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
    i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
    i5 = FixedParam(1.5, 'I5', description='Column effective length factor', symbol='K')
    i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')

    fixed_params = [i1, i2, i3, i4, i5, i6,]

    # define design parameters
    d1 = DesignParam(15.0, 'D1', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
    d2 = DesignParam(0.0, 'D2', universe=[0.0, 50.0], variable_type='FLOAT', description='lean angle', symbol='theta')
    design_params = [d1, d2]

    # T1,T2 distribution (Uniform)
    center = np.array([450, 425])
    Range = np.array([100, 100]) / (20 * 0.25)
    Requirement_1 = UniformFunc(center, Range, 'temp')

    # define input specifications
    s1 = InputSpec(450, 'S1', universe=[325, 550], variable_type='FLOAT', cov_index=0,
                description='nacelle temperature', distribution=Requirement_1,
                symbol='T1', inc=-1e-1, inc_type='rel')
    s2 = InputSpec(425, 'S2', universe=[325, 550], variable_type='FLOAT', cov_index=1,
                description='gas surface temperature', distribution=Requirement_1,
                symbol='T2', inc=+1e-1, inc_type='rel')

    # BX,BY distribution (Uniform)
    center = np.array([100, 100])
    Range = np.array([5, 5])
    Requirement_2 = UniformFunc(center, Range, 'force')

    # define input specifications
    s3 = InputSpec(120, 'S3', universe=[50, 200], variable_type='FLOAT', cov_index=0,
                description='bearing load x', distribution=Requirement_2,
                symbol='BX', inc=+1e-0, inc_type='rel')
    s4 = InputSpec(120, 'S4', universe=[50, 200], variable_type='FLOAT', cov_index=1,
                description='bearing load y', distribution=Requirement_2,
                symbol='BY', inc=+1e-0, inc_type='rel')

    input_specs = [s1, s2, s3, s4]
    distributions = [Requirement_1, Requirement_2]

    # define the behaviour models
    # this is the force model
    class B1(Behaviour):
        def __call__(self, T1, T2, h, theta, alpha, E, r1, r2, Ts):
            coeffs = [0.95, 1.05, 0.97]
            coeffs = 3 * [1.0, ]
            w_nominal = 60.0
            n_struts_nominal = 12

            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            force_thermal = (E * w_nominal * h * alpha) * ((T2 * coeffs[0] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(
                np.deg2rad(theta)) / length

            # force_bearing = (BY*np.cos(np.deg2rad(theta)) + BX*np.sin(np.deg2rad(theta)))*1e4/n_struts_nominal

            self.threshold = force_thermal / 1000


    # this is the buckling model
    class B2(Behaviour):
        def __call__(self, w, h, theta, E, r1, r2, K):
            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            f_buckling = ((np.pi ** 2) * E * w * (h ** 3)) / (12 * ((K * length) ** 2))
            self.decided_value = f_buckling / 1000

        def inv_call(self,decided_value, h, theta, E, r1, r2, K):
            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            f_buckling = decided_value * 1000
            w = f_buckling / (((np.pi ** 2) * E * (h ** 3)) / (12 * ((K * length) ** 2)))

            self.inverted = w


    # this is the stress model
    class B3(Behaviour):
        def __call__(self, T1, T2, h, theta, alpha, E, r1, r2, Ts):
            coeffs = [0.95, 1.05, 0.97]
            coeffs = 3 * [1.0, ]
            n_struts_nominal = 12

            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            sigma_a_thermal = (E * alpha) * ((T2 * coeffs[1] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(np.deg2rad(theta)) / length
            sigma_m_thermal = (3 / 2) * ((E * h) / (length ** 2)) * (
                    alpha * ((T2 * coeffs[2] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.sin(np.deg2rad(theta)))
            
            # F_axial = (BY*np.cos(np.deg2rad(theta)) + BX*np.sin(np.deg2rad(theta)))*1e4/n_struts_nominal
            # F_bending = (BX*np.cos(np.deg2rad(theta)) - BY*np.sin(np.deg2rad(theta)))*1e4/n_struts_nominal
            # M = F_bending*length# + (BX*1e4*r1/n_struts_nominal)

            # sigma_a_bearing = F_axial / (w*h)
            # sigma_m_bearing = M*6/(w*h**2)

            self.threshold = max([sigma_a_thermal, sigma_m_thermal])


    # this is the material model
    class B4(Behaviour):
        def __call__(self, material):

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

            chosen_mat = material_dict[material]

            self.intermediate = [chosen_mat['rho'], chosen_mat['cost']]
            self.decided_value = chosen_mat['sigma_y']


    # this is the simulation model
    class B5(Behaviour):
        def __call__(self, n_struts, w, h, theta, BX, BY, E, r1, r2, id=None):
            # Define arguments needed to build surrogate
            args = [n_struts,w,h,theta,BX,BY] # E, r1, r2 are fixed parameters
            if self.surrogate_available:
                return Behaviour.__call__(self,*args) # once surrogate is trained, the base class method will terminate here

            from fem_scirpts.nx_trs_script import run_nx_simulation
            from fem_scirpts.plot_results import postprocess_vtk, plotmesh, plotresults

            _, vtk_filename, success = run_nx_simulation(w,h,theta,n_struts,r_hub=r1,r_shroud=r2,
                youngs_modulus=E,bearing_x=BX,bearing_y=BY,pid=id)

            # base_path = os.getcwd() # Working directory
            # folder = 'Nastran_output'
            # output_path = os.path.join(base_path, 'examples', 'CAD', folder)
            # if not os.path.exists(output_path):
            #     os.makedirs(output_path)

            # plotmesh(vtk_filename,output_path)
            # plotresults(1,vtk_filename,'POINT',output_path)
            # plotresults(1,vtk_filename,'CELL',output_path)

            # displacement = postprocess_vtk(1,vtk_filename,'POINT')
            if success:
                stress = postprocess_vtk(1,vtk_filename,'CELL')
            else:
                stress = np.nan

            self.decided_value = stress # for the decision node
        

    # this is the weight model
    class B6(Behaviour):
        def __call__(self, rho, cost_coeff, w, n_struts, h, theta, r1, r2):

            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            weight = rho * w * h * length * n_struts
            cost = weight * cost_coeff
            self.performance = [weight, cost]


    b1 = B1(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B1')
    b2 = B2(n_i=0, n_p=2, n_dv=1, n_tt=0, key='B2')
    b3 = B3(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B3')
    b4 = B4(n_i=2, n_p=0, n_dv=1, n_tt=0, key='B4')
    b5 = B5(n_i=0, n_p=0, n_dv=1, n_tt=0, key='B5')
    b6 = B6(n_i=0, n_p=2, n_dv=0, n_tt=0, key='B6')

    # Define decision nodes and a model to convert to decided values
    decision_1 = Decision(universe=widths, variable_type='ENUM', key='decision_1',
                        direction='must_not_exceed', decided_value_model=b2, n_nodes=1,
                        description='the vane width selection')

    decision_2 = Decision(universe=['Steel','Inconel','Titanium'], variable_type='ENUM', key='decision_2',
                        direction='must_not_exceed', decided_value_model=b4, n_nodes=1, 
                        description='The type of material')

    decision_3 = Decision(universe=list(range(6,30+1,2)), variable_type='ENUM', key='decision_3',
                        direction='must_not_exceed', decided_value_model=None, n_nodes=1,
                        description='the number of struts')

    decisions = [decision_1, decision_2, decision_3]
    behaviours = [b1, b2, b3, b4, b5, b6,]

    # Define margin nodes
    e1 = MarginNode('E1', direction='must_not_exceed')
    e2 = MarginNode('E2', direction='must_not_exceed')
    e3 = MarginNode('E3', direction='must_not_exceed')
    margin_nodes = [e1, e2, e3,]

    # Define performances
    p1 = Performance('P1', direction='less_is_better')
    p2 = Performance('P2', direction='less_is_better')
    performances = [p1, p2,]

    # Define the MAN
    class MAN(MarginNetwork):

        def randomize(self):
            s1 = self.input_specs[0]  # T1 (stochastic)
            s2 = self.input_specs[1]  # T2 (stochastic)
            s3 = self.input_specs[2]  # BX (stochastic)
            s4 = self.input_specs[3]  # BY (stochastic)
            s1.random()
            s2.random()
            s3.random()
            s4.random()

        def forward(self,num_threads=1,recalculate_decisions=False,allocate_margin=False,strategy='min_excess',outputs=['dv','dv','dv']):

            # retrieve MAN components
            d1 = self.design_params[0]  # h
            d2 = self.design_params[1]  # theta

            s1 = self.input_specs[0]  # T1 (stochastic)
            s2 = self.input_specs[1]  # T2 (stochastic)
            s3 = self.input_specs[2]  # BX (stochastic)
            s4 = self.input_specs[3]  # BY (stochastic)

            i1 = self.fixed_params[0]  # alpha
            i2 = self.fixed_params[1]  # E
            i3 = self.fixed_params[2]  # r1
            i4 = self.fixed_params[3]  # r2
            i5 = self.fixed_params[4]  # K
            i6 = self.fixed_params[5]  # Ts

            b1 = self.behaviours[0]  # calculates axial force
            b2 = self.behaviours[1]  # calculates buckling load
            b3 = self.behaviours[2]  # calculates bending and axial stresses
            b4 = self.behaviours[3]  # convert material index to yield stress, density, and cost
            b5 = self.behaviours[4]  # calculates von mises stress from a simulation
            b6 = self.behaviours[5]  # calculates weight and cost

            decision_1 = self.decisions[0]  # select the width of the vane based on the maximum supported buckling load
            decision_2 = self.decisions[1]  # select the material based on the max stress
            decision_3 = self.decisions[2]  # select the number of struts based on center displacement and max stress

            e1 = self.margin_nodes[0]  # margin against buckling (F,F_buckling)
            e2 = self.margin_nodes[1]  # margin against axial or bending failure (max(sigma_a,sigma_m),sigma_y)
            e3 = self.margin_nodes[2]  # margin against yielding due to bearing loads(sigma,sigma_y)

            p1 = self.performances[0]  # weight
            p2 = self.performances[1]  # cost

            # Execute behaviour models
            # T1, T2, h, theta, alpha, E, r1, r2, Ts
            b1(s1.value, s2.value, d1.value, d2.value, i1.value, i2.value, i3.value, i4.value, i6.value)
            # Execute decision node for width: w, h, theta, E, r1, r2, K
            args = [
                self.design_params[0].value, # h
                self.design_params[1].value, # theta
                self.fixed_params[1].value, # E
                self.fixed_params[2].value, # r1
                self.fixed_params[3].value, # r2
                self.fixed_params[4].value, # K
            ]
            decision_1(b1.threshold, recalculate_decisions, allocate_margin, strategy, num_threads, outputs[0],*args)
            # invert decided value: decided_value, h, theta, E, r1, r2, K
            b2.inv_call(decision_1.output_value, d1.value, d2.value, i2.value, i3.value, i4.value, i5.value)

            # T1, T2, h, theta, alpha, E, r1, r2, Ts)
            b3(s1.value, s2.value, d1.value, d2.value, i1.value, i2.value, i3.value, i4.value, i6.value)
            # Execute decision node for material and translate to yield stress: material
            decision_2(b3.threshold, recalculate_decisions, allocate_margin, strategy, num_threads, outputs[1])
            # invert decided value: decided_value, h, theta, E, r1, r2, K
            b4.inv_call(decision_2.output_value)
            
            # args = [
            #     b2.inverted, # w
            #     self.design_params[0].value, # h
            #     self.design_params[1].value, # theta
            #     self.input_specs[2].value, # BX
            #     self.input_specs[3].value, # BY
            #     self.fixed_params[1].value, # E
            #     self.fixed_params[2].value, # r1
            #     self.fixed_params[3].value, # r2
            # ]
            # kwargs = {
            #     'id' : self.key
            # }
            # # Execute decision node for struts: n_struts, w, h, theta, BX, BY, E, r1, r2
            # decision_3(decision_2.output_value, override_decisions, recalculate_decisions, num_threads, outputs[2], *args, **kwargs)
            # # invert decided value: decided_value, w, h, theta, BX, BY, E, r1, r2
            # b5.inv_call(decision_3.output_value, b2.inverted, d1.value, d2.value, s3.value, s4.value, i2.value, i3.value, i4.value)

            # calculate n_struts sigma_y, w, h, theta, BX, BY, E, r1, r2
            b5.inv_call(decision_2.output_value, b2.inverted, d1.value, d2.value, s3.value, s4.value, i2.value, i3.value, i4.value)
            # Execute decision node for struts: n_struts, w, h, theta, BX, BY, E, r1, r2
            decision_3(b5.inverted, recalculate_decisions, allocate_margin, strategy, num_threads, outputs[2])

            # Compute excesses
            e1(b1.threshold, decision_1.decided_value)
            e2(b3.threshold, decision_2.decided_value)
            e3(b5.inverted, decision_3.decided_value)

            # Compute performances
            # rho, cost_coeff, w, n_struts, h, theta, r1, r2
            b6(b4.intermediate[0], b4.intermediate[1], b2.inverted, decision_3.output_value, d1.value, d2.value, i3.value, i4.value)
            p1(b6.performance[0])
            p2(b6.performance[1])


    man = MAN(design_params, input_specs, fixed_params,
              behaviours, decisions, margin_nodes, performances, 'MAN_1')

    return man

def evaluate_design(i: int,design, mans: List[MarginNetwork], n_epochs: int,
                    base_folder: str, man_name: str, process_ids: List[int]=None):

    # Select a man to forward based on process id
    if process_ids == None:
        pid = None
        man = mans[0]
    else:
        pid = mp.current_process()._identity[0] - process_ids[0]
        man = mans[pid]

    man.nominal_design_vector = np.array(design[:len(man.design_params)])
    man.reset()
    man.reset_inputs('all')
    man.reset_outputs()

    # Perform Monte-Carlo simulation
    for n in range(n_epochs):
        man.randomize()
        man.init_decisions()
        man.decision_vector = design[len(man.design_params):]
        man.allocate_margins('manual')
        man.forward()
        man.compute_impact()
        man.compute_absorption()

    folder = os.path.join(base_folder,'d%i'%i)
    check_folder(folder) # create directory if it does not exist
    man.save(man_name,folder=folder,results_only=True)

if __name__ == "__main__":

    # get the man object for the fea problem and load it
    base_folder = os.path.join('data','strut_fea','opt')
    img_folder = os.path.join('images','strut_fea','opt')
    check_folder(img_folder)
    man = get_man_DOE()
    man_name = 'strut_comb'
    man.load(man_name,folder=base_folder)

    man.input_specs[2].inc_user = 20
    man.input_specs[3].inc_user = 20

    # for fea problem
    lb = np.array([-0.5973142383857094, 0.0])
    ub = np.array([3.0949267424315656, 7.4235687891259134])

    # # get the man object and load it
    # base_folder = os.path.join('data','strut','opt')
    # img_folder = os.path.join('images','strut','opt')
    # check_folder(img_folder)
    # man = get_man()
    # man_name = 'strut_s'
    # man.load(man_name,folder=base_folder)
    # # for simple problem
    # lb = np.array([0.00, 2.5])
    # ub = np.array([0.075, 3.5])

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

    universe = universe_d + man.universe_decision
    design_doe = list(itertools.product(*universe))
    n_designs = len(design_doe)
    n_epochs = 50

    # 1026, 1027, 1028, 1029
    # doe_reduced = [design_doe[d] for d in [1973,1972]] # debug 1553,1973,1552,1972
    #---------------------------------------------------
    # Evaluate all the different designs
    # Parallel computation if num_threads > 1
    num_threads = 8

    man_objs = []
    for pid in range(num_threads):
        man_objs += [deepcopy(man)]

    kwargs = {'n_epochs': n_epochs, 
              'base_folder': base_folder, 
              'man_name': man_name}

    vargs_iterator = [[i,design,] for i,design in enumerate(design_doe)]
    vkwargs_iterator = [{},] * n_designs
    fargs = []
    fkwargs = kwargs
    fkwargs['mans'] = man_objs

    # d1553 d1973 d1552 d1972

    # results = parallel_sampling(evaluate_design,vargs_iterator,vkwargs_iterator,fargs,fkwargs,num_threads=num_threads)
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

            # calculate means
            I = man.impact_matrix.values[n,:,:]
            A = man.absorption_matrix.values[n,:,:]

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
                nodal_dict[decision.key] = design[d_i+1+decision_i]
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
        mean_i = np.mean(man.impact_matrix.values,axis=(1,2)) # absorption
        mean_a = np.nanmean(man.absorption_matrix.values,axis=(1,2)) # impact

        # construct dict
        total_dict = {'id':id}
        for d_i,d in enumerate(man.design_params):
            total_dict[d.key] = design[d_i]
        for decision_i,decision in enumerate(man.decisions):
            total_dict[decision.key] = design[d_i+1+decision_i]
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
    plot(data_pd, show_link=False, filename = os.path.join(img_folder,'pcp_doe'), auto_open=True)

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
    plot(data_pd_nodal, show_link=False, filename = os.path.join(img_folder,'pcp_nodal'), auto_open=True)


    #---------------------------------------------------
    # export data

    ## doe
    df_export = df.drop(['material'],axis=1)
    df_export = df_export.rename(columns={'dummy':'material'})

    # shift column 'Name' to first position
    material = df_export.pop('material')
    df_export.insert(3, 'material', material)
    
    ## nodal
    df_export_nodal = df_nodal.drop(['material'],axis=1)
    df_export_nodal = df_export_nodal.rename(columns={'dummy':'material'})

    # shift column 'Name' to first position
    material = df_export_nodal.pop('material')
    df_export_nodal.insert(3, 'material', material)

    # export to csv
    df_export.to_csv(os.path.join(base_folder,'doe_data.csv'))
    df_export_nodal.to_csv(os.path.join(base_folder,'nodal_data.csv'))