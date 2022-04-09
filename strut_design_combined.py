import sys

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from random import randint

from mvm import Design
from mvm import FixedParam, DesignParam, InputSpec, Behaviour, Performance, MarginNode, MarginNetwork, Decision
from mvm import GaussianFunc, UniformFunc

if __name__ == "__main__":

    num_threads = 1

    # define fixed parameters
    i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
    i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
    i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
    i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
    i5 = FixedParam(1.0, 'I5', description='Column effective length factor', symbol='K')
    i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')

    fixed_params = [i1, i2, i3, i4, i5, i6,]

    # define design parameters
    d1 = DesignParam(100.0, 'D1', universe=[70.0, 130.0], variable_type='FLOAT', description='vane length', symbol='w')
    d2 = DesignParam(15.0, 'D2', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
    d3 = DesignParam(10.0, 'D3', universe=[0.0, 30.0], variable_type='FLOAT', description='lean angle', symbol='theta')
    design_params = [d1, d2, d3]

    # T1,T2 distribution (Uniform)
    center = np.array([450, 425])
    Range = np.array([100, 100]) / (20 * 0.25)
    Requirement_1 = UniformFunc(center, Range, 'temp')

    # define input specifications
    s1 = InputSpec(450, 'S1', universe=[325, 550], variable_type='FLOAT', cov_index=0,
                description='nacelle temperature', distribution=Requirement_1,
                symbol='T1', inc=-1e-0, inc_type='rel')
    s2 = InputSpec(425, 'S2', universe=[325, 550], variable_type='FLOAT', cov_index=1,
                description='gas surface temperature', distribution=Requirement_1,
                symbol='T2', inc=+1e-0, inc_type='rel')

    # BX,BY distribution (Uniform)
    center = np.array([0, 0])
    Range = np.array([450, 450])
    Requirement_2 = UniformFunc(center, Range, 'force')

    # define input specifications
    s3 = InputSpec(200, 'S3', universe=[0, 450], variable_type='FLOAT', cov_index=0,
                description='bearing load x', distribution=Requirement_2,
                symbol='BX', inc=+1e-0, inc_type='rel')
    s4 = InputSpec(200, 'S4', universe=[0, 450], variable_type='FLOAT', cov_index=1,
                description='bearing load y', distribution=Requirement_2,
                symbol='BY', inc=+1e-0, inc_type='rel')

    input_specs = [s1, s2, s3, s4]

    # define the behaviour models

    # this is the length model
    class B1(Behaviour):
        def __call__(self, theta, r1, r2):
            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)
            self.intermediate = length


    # this is the weight model
    class B2(Behaviour):
        def __call__(self, rho, w, h, L, n_struts, cost_coeff):
            weight = rho * w * h * L * n_struts
            cost = weight * cost_coeff
            self.performance = [weight, cost]


    # this is the axial stress model
    class B3(Behaviour):
        def __call__(self, alpha, E, r1, r2, Ts, T1, T2, w, h, theta, L):
            coeffs = [0.95, 1.05, 0.97]
            coeffs = 3 * [1.0, ]

            force = (E * w * h * alpha) * ((T2 * coeffs[0] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(
                np.deg2rad(theta)) / L
            sigma_a = (E * alpha) * ((T2 * coeffs[1] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(np.deg2rad(theta)) / L
            self.threshold = [force / 1000, sigma_a]


    # this is the bending stress model
    class B4(Behaviour):
        def __call__(self, alpha, E, r1, r2, Ts, T1, T2, h, theta, L):
            coeffs = [0.95, 1.05, 0.97]
            coeffs = 3 * [1.0, ]

            sigma_m = (3 / 2) * ((E * h) / (L ** 2)) * (
                    alpha * ((T2 * coeffs[2] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.sin(np.deg2rad(theta)))
            self.threshold = sigma_m


    # this is the buckling model
    class B5(Behaviour):
        def __call__(self, E, K, w, h, L):
            f_buckling = ((np.pi ** 2) * E * w * (h ** 3)) / (12 * ((K * L) ** 2))
            self.decided_value = f_buckling / 1000


    class B6(Behaviour):
        def __call__(self, material):

            material_dict = {
                'Inconel'   : {
                    'sigma_y' : 92, # MPa
                    'rho' : 11.95e-06, # kg/mm3
                    'cost' : 0.1 # USD/kg
                    },
                'Titanium'  : {
                    'sigma_y' : 828, # MPa
                    'rho' : 4.43e-06, # kg/mm3
                    'cost' : 1.10 # USD/kg
                    },
            }

            # generate 10 materials by linearly interpolating
            df = pd.DataFrame(columns=material_dict['Inconel'].keys(), index=range(11), dtype=float)
            df.iloc[0] = material_dict['Inconel']
            df.iloc[-1] = material_dict['Titanium']
            df.interpolate(method='linear',axis=0,inplace=True)
            material_dict_interp = df.transpose().to_dict()

            chosen_mat = material_dict_interp[material]

            self.intermediate = [chosen_mat['rho'], chosen_mat['cost']]
            self.decided_value = [chosen_mat['sigma_y'],chosen_mat['sigma_y'],]

            return self.decided_value


    # this is the simulation model
    class B7(Behaviour):
        def __call__(self, n_struts, w, h, theta, BX, BY, E, r1, r2, id=None):
            # Define arguments needed to build surrogate
            args = [n_struts,w,h,theta,BX,BY]
            if self.surrogate_available:
                return Behaviour.__call__(self,*args) # once surrogate is trained, the base class method will terminate here

            from fem_scirpts.nx_trs_script import run_nx_simulation
            from fem_scirpts.plot_results import postprocess_vtk, plotmesh, plotresults

            _, vtk_filename = run_nx_simulation(w,h,theta,n_struts,r_hub=r1,r_shroud=r2,
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
            stress = postprocess_vtk(1,vtk_filename,'CELL')

            self.decided_value = stress # for the decision node

            return self.decided_value


    b1 = B1('B1')
    b2 = B2('B2')
    b3 = B3('B3')
    b4 = B4('B4')
    b5 = B5('B5')
    b6 = B6('B6')
    b7 = B7('B7')

    # Define decision nodes and a model to convert to decided values
    decision_1 = Decision(universe=list(range(11)), variable_type='ENUM', key='decision_1',
                        direction=['must_not_exceed','must_not_exceed'], decided_value_model=b6, 
                        n_nodes=2, description='The type of material')

    decision_2 = Decision(universe=[6,18], variable_type='INT', key='decision_2',
                        direction='must_exceed', decided_value_model=b7, n_nodes=1,
                        description='the number of struts')

    decisions = [decision_1, decision_2,]
    behaviours = [b1, b2, b3, b4, b5, b6, b7,]

    # Define margin nodes
    e1 = MarginNode('E1', direction='must_not_exceed')
    e2 = MarginNode('E2', direction='must_not_exceed')
    e3 = MarginNode('E3', direction='must_not_exceed')
    e4 = MarginNode('E4', direction='must_exceed')
    margin_nodes = [e1, e2, e3, e4,]

    # Define performances
    p1 = Performance('P1', direction='less_is_better')
    p2 = Performance('P2', direction='less_is_better')
    performances = [p1, p2,]

    # Define the MAN
    class MAN(MarginNetwork):

        def randomize(self):
            Requirement_1.random()
            Requirement_2.random()
            s1.random()
            s2.random()
            s3.random()
            s4.random()

        def forward(self, num_threads=1,recalculate_decisions=False,override_decisions=False):

            def random_with_N_digits(n):
                range_start = 10**(n-1)
                range_end = (10**n)-1
                return randint(range_start, range_end)

            # retrieve MAN components
            d1 = self.design_params[0]  # w
            d2 = self.design_params[1]  # h
            d3 = self.design_params[2]  # theta

            s1 = self.input_specs[0]  # T1 (stochastic)
            s2 = self.input_specs[1]  # T2 (stochastic)
            s3 = self.input_specs[2]  # T3 (stochastic)
            s4 = self.input_specs[3]  # T4 (stochastic)

            i1 = self.fixed_params[0]  # alpha
            i2 = self.fixed_params[1]  # E
            i3 = self.fixed_params[2]  # r1
            i4 = self.fixed_params[3]  # r2
            i5 = self.fixed_params[4]  # K
            i6 = self.fixed_params[5]  # Ts

            b1 = self.behaviours[0]  # calculates length
            b2 = self.behaviours[1]  # calculates weight and cost
            b3 = self.behaviours[2]  # calculates axial force and stress
            b4 = self.behaviours[3]  # calculates bending stress
            b5 = self.behaviours[4]  # calculates buckling load
            b6 = self.behaviours[5]  # convert material index to yield stress, density, and cost
            b7 = self.behaviours[6]  # runs the nx simulation and returns center displacement and max stress

            decision_1 = self.decisions[0]  # select a material based on maximum bending and axial stress
            decision_2 = self.decisions[1]  # select the number of struts based on center displacement and max stress

            e1 = self.margin_nodes[0]  # margin against buckling (F,F_buckling)
            e2 = self.margin_nodes[1]  # margin against axial failure (sigma_a,sigma_y)
            e3 = self.margin_nodes[2]  # margin against bending failure (sigma_m,sigma_y)
            e4 = self.margin_nodes[3]  # margin against yielding due to bearing loads(sigma,sigma_y)

            p1 = self.performances[0]  # weight
            p2 = self.performances[1]  # cost

            # Execute behaviour models
            # theta, r1, r2
            b1(d3.value, i3.value, i4.value)
            # alpha, E, r1, r2, Ts, T1, T2, w, h, theta, L
            b3(i1.value, i2.value, i3.value, i4.value, i6.value, s1.value, s2.value, d1.value, d2.value, d3.value, b1.intermediate)
            # alpha, E, r1, r2, Ts, T1, T2, h, theta, L
            b4(i1.value, i2.value, i3.value, i4.value, i6.value, s1.value, s2.value, d2.value, d3.value, b1.intermediate)
            # E, K, w, h, L
            b5(i2.value, i5.value, d1.value, d2.value, b1.intermediate)

            # Execute decision node and translation model
            decision_1([b3.threshold[1], b4.threshold], override_decisions, recalculate_decisions, num_threads)
            # E, K, w, h, L
            b6(decision_1.selection_value)

            args = [
                self.design_params[0].value, # w
                self.design_params[1].value, # h
                self.design_params[2].value, # theta
                self.input_specs[2].value, # BX
                self.input_specs[3].value, # BY
                self.fixed_params[1].value, # E
                self.fixed_params[2].value, # r1
                self.fixed_params[3].value, # r2
            ]

            kwargs = {
                'id' : self.key
                # 'id' : random_with_N_digits(10)
            }

            decision_2(decision_1.decided_value[0], override_decisions, recalculate_decisions, num_threads, *args, **kwargs)

            # Compute excesses
            e1(b3.threshold[0], b5.decided_value)
            e2(b3.threshold[1], decision_1.decided_value[0])
            e3(b4.threshold, decision_1.decided_value[1])
            e4(decision_1.decided_value[0], decision_2.decided_value)

            # Compute performances
            # (rho, w, h, L, n_struts, cost_coeff)
            b2(b6.intermediate[0], d1.value, d2.value, b1.intermediate, decision_2.selection_value, b6.intermediate[1])
            p1(b2.performance[0])
            p2(b2.performance[1])


    man = MAN(design_params, input_specs, fixed_params,
            behaviours, decisions, margin_nodes, performances, 'MAN_1')

    # # Create surrogate model for estimating threshold performance
    # man.train_performance_surrogate(n_samples=100, sm_type='LS', sampling_freq=1, 
    #     num_threads=num_threads, bandwidth=[1e-3, ])
    # man.save('strut_comb',folder='strut_fea')

    # # Train surrogate for B7
    # variable_dict = {
    #     'n_struts' : {'type' : 'INT', 'limits' : decision_2.universe},
    #     'w' : {'type' : 'FLOAT', 'limits' : d1.universe},
    #     'h' : {'type' : 'FLOAT', 'limits' : d2.universe},
    #     'theta' : {'type' : 'FLOAT', 'limits' : d3.universe},
    #     'BX' : {'type' : 'FLOAT', 'limits' : s3.universe},
    #     'BY' : {'type' : 'FLOAT', 'limits' : s4.universe},
    #     'E' : {'type' : 'fixed', 'limits' : i2.value},
    #     'r1' : {'type' : 'fixed', 'limits' : i3.value},
    #     'r2' : {'type' : 'fixed', 'limits' : i4.value},
    # }

    # b7.train_surrogate(variable_dict,n_outputs=1,n_samples=1000,num_threads=num_threads)
    # man.save('strut_comb',folder='strut_fea')
    # sys.exit(0)

    # load the MAN
    man.load('strut_comb',folder='strut_fea')
    # man.train_performance_surrogate(ext_samples=(man.xt,man.yt), sm_type='LS')

    man.init_decisions(num_threads=num_threads)
    man.forward()
    man.view_perf(e_indices=[0, 1], p_index=0)
    man.view_perf(e_indices=[0, 1], p_index=1)

    # Perform Monte-Carlo simulation
    n_epochs = 1000
    for n in range(n_epochs):
        sys.stdout.write("Progress: %d%%   \r" % ((n / n_epochs) * 100))
        sys.stdout.flush()

        man.randomize()
        man.init_decisions(num_threads=num_threads)
        man.forward()
        man.compute_impact()
        man.compute_absorption(num_threads=num_threads)

    man.save('strut_fem_100',folder='strut_fea')

    # View distribution of excess
    man.margin_nodes[0].excess.view(xlabel='E1')
    man.margin_nodes[1].excess.view(xlabel='E2')
    man.margin_nodes[2].excess.view(xlabel='E3')
    man.margin_nodes[3].excess.view(xlabel='E4')

    # View distribution of Impact on Performance
    man.impact_matrix.view(0, 0, xlabel='E1,P1')
    man.impact_matrix.view(1, 0, xlabel='E2,P1')
    man.impact_matrix.view(2, 0, xlabel='E3,P1')
    man.impact_matrix.view(3, 0, xlabel='E4,P1')
    man.impact_matrix.view(0, 1, xlabel='E1,P2')
    man.impact_matrix.view(1, 1, xlabel='E2,P2')
    man.impact_matrix.view(2, 1, xlabel='E3,P2')
    man.impact_matrix.view(3, 1, xlabel='E4,P2')

    man.deterioration_vector.view(0, xlabel='S1')
    man.deterioration_vector.view(1, xlabel='S2')
    man.deterioration_vector.view(1, xlabel='S3')
    man.deterioration_vector.view(1, xlabel='S4')

    man.absorption_matrix.view(3, 0, xlabel='E4,S1')
    man.absorption_matrix.view(3, 1, xlabel='E4,S1')
    man.absorption_matrix.view(3, 2, xlabel='E4,S2')
    man.absorption_matrix.view(3, 3, xlabel='E4,S3')

    man.absorption_matrix.view(0, 3, xlabel='E1,S4')
    man.absorption_matrix.view(1, 3, xlabel='E1,S4')
    man.absorption_matrix.view(2, 3, xlabel='E1,S4')
    man.absorption_matrix.view(3, 3, xlabel='E1,S4')

    man.absorption_matrix.view(3, 0, xlabel='E4,S1')
    man.absorption_matrix.view(3, 1, xlabel='E4,S2')
    man.absorption_matrix.view(3, 2, xlabel='E4,S3')
    man.absorption_matrix.view(3, 3, xlabel='E4,S4')

    # display the margin value plot
    man.compute_mvp('scatter')

    sys.exit(0)

    # Effect of alternative designs
    n_designs = 100
    n_epochs = 10
    lb = np.array(man.universe_d)[:, 0]
    ub = np.array(man.universe_d)[:, 1]
    design_doe = Design(lb, ub, n_designs, 'LHS')

    # create empty figure
    fig, ax = plt.subplots(figsize=(7, 8))
    ax.set_xlabel('Impact on performance')
    ax.set_ylabel('Change absoption capability')

    for d, design in enumerate(design_doe.unscale()):
        man.nominal_design_vector = design
        man.reset()
        man.reset_outputs()

        # Perform Monte-Carlo simulation
        for n in range(n_epochs):
            sys.stdout.write("Progress: %d%%   \r" % ((d * n_epochs + n) / (n_designs * n_epochs) * 100))
            sys.stdout.flush()

            man.randomize()
            man.init_decisions(w=d1.value, h=d2.value, theta=d3.value, E=i1.value,
                r1=i2.value, r2=i3.value, BX=s1.value, BY=s2.value)
            man.forward()
            man.compute_impact()
            man.compute_absorption()

        # Extract x and y
        x = np.mean(man.impact_matrix.values,
                    axis=(1, 2)).ravel()  # average along performance parameters (assumes equal weighting)
        y = np.mean(man.absorption_matrix.values,
                    axis=(1, 2)).ravel()  # average along input specs (assumes equal weighting)

        # plot the results
        color = np.random.random((1, 3))
        ax.scatter(x, y, c=color)

    plt.show()
