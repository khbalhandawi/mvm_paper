import sys

import matplotlib.pyplot as plt
import numpy as np
import os

from mvm import Design
from mvm import FixedParam, DesignParam, InputSpec, Behaviour, Performance, MarginNode, MarginNetwork, Decision
from mvm import GaussianFunc, UniformFunc
from mvm import nearest

if __name__ == '__main__':

    num_threads = 1

    # define fixed parameters
    i1 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
    i2 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
    i3 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
    i4 = FixedParam(460.0, 'I5', description='yield strength', symbol='sigma_y')
    i5 = FixedParam(8.19e-06, 'I6', description='density', symbol='rho')
    i6 = FixedParam(0.46, 'I7', description='cost per unit mass', symbol='cost')
    i7 = FixedParam(0.4, 'I8', description='maximum allowed deflection', symbol='d_max')

    fixed_params = [i1, i2, i3, i4, i5, i6, i7,]

    # define design parameters
    d1 = DesignParam(100.0, 'D1', universe=[70.0, 130.0], variable_type='FLOAT', description='vane length', symbol='w')
    d2 = DesignParam(15.0, 'D2', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
    d3 = DesignParam(10.0, 'D3', universe=[0.0, 30.0], variable_type='FLOAT', description='lean angle', symbol='theta')
    design_params = [d1, d2, d3]

    # # BX,BY distribution (Gaussian)
    # center = np.array([450, 425])
    # Sigma = np.array([
    #     [100, 25],
    #     [75, 100],
    # ]) / (20 * 0.1)
    # Requirement = GaussianFunc(center, Sigma, 'temp')

    # BX,BY distribution (Uniform)
    center = np.array([0, 0])
    Range = np.array([450, 450])
    Requirement = UniformFunc(center, Range, 'temp')

    # define input specifications
    s1 = InputSpec(200, 'S1', universe=[-450, 450], variable_type='FLOAT', cov_index=0,
                description='bearing load x', distribution=Requirement,
                symbol='BX', inc=+1e-0, inc_type='rel')
    s2 = InputSpec(200, 'S2', universe=[-450, 450], variable_type='FLOAT', cov_index=1,
                description='bearing load y', distribution=Requirement,
                symbol='BY', inc=+1e-0, inc_type='rel')
    input_specs = [s1, s2]


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


    # this is the simulation model
    class B3(Behaviour):
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

            displacement = postprocess_vtk(1,vtk_filename,'POINT')
            stress = postprocess_vtk(1,vtk_filename,'CELL')

            self.decided_value = [stress, displacement] # for the decision node

            return [stress, displacement]


    b1 = B1('B1')
    b2 = B2('B2')
    b3 = B3('B3')

    # Define decision nodes and a model to convert to decided values
    decision_1 = Decision(universe=[6,18], variable_type='INT', key='decision_1',
                        direction=['must_exceed','must_exceed'], decided_value_model=b3, n_nodes=2,
                        description='the number of struts')

    decisions = [decision_1, ]
    behaviours = [b1, b2, b3,]
        
    # Define margin nodes
    e1 = MarginNode('E1', direction='must_exceed')
    e2 = MarginNode('E2', direction='must_exceed')
    margin_nodes = [e1, e2,]

    # Define performances
    p1 = Performance('P1', direction='less_is_better')
    p2 = Performance('P2', direction='less_is_better')
    performances = [p1, p2]


    # Define the MAN
    class MAN(MarginNetwork):

        def randomize(self):
            Requirement.random()
            s1.random()
            s2.random()

        def forward(self, num_threads=1,recalculate_decisions=False,override_decisions=False):
            # retrieve MAN components
            d1 = self.design_params[0]  # w
            d2 = self.design_params[1]  # h
            d3 = self.design_params[2]  # theta

            s1 = self.input_specs[0]  # T1 (stochastic)
            s2 = self.input_specs[1]  # T2 (stochastic)

            i1 = self.fixed_params[0]  # E
            i2 = self.fixed_params[1]  # r1
            i3 = self.fixed_params[2]  # r2
            i4 = self.fixed_params[3]  # sigma_y
            i5 = self.fixed_params[4]  # rho
            i6 = self.fixed_params[5]  # cost
            i7 = self.fixed_params[6]  # d_max

            b1 = self.behaviours[0]  # calculates length
            b2 = self.behaviours[1]  # calculates weight and cost
            b3 = self.behaviours[2]  # runs the nx simulation and returns center displacement and max stress

            decision_1 = self.decisions[0]  # select the number of struts based on center displacement and max stress

            e1 = self.margin_nodes[0]  # margin against yielding (sigma,sigma_y)
            e2 = self.margin_nodes[1]  # margin against deflection (d,d_max)

            p1 = self.performances[0]  # weight
            p2 = self.performances[1]  # cost

            # Execute behaviour models
            b1(d3.value, i2.value, i3.value)

            # Execute decision node and translation model
            args = [
                self.design_params[0].value, # w
                self.design_params[1].value, # h
                self.design_params[2].value, # theta
                self.input_specs[0].value, # BX
                self.input_specs[1].value, # BY
                self.fixed_params[0].value, # E
                self.fixed_params[1].value, # r1
                self.fixed_params[2].value, # r2
            ]

            kwargs = {
                'id' : self.key
            }

            decision_1([i4.value,i7.value], override_decisions, 
                recalculate_decisions, num_threads, *args, **kwargs)

            # Compute excesses
            e1(i4.value, decision_1.decided_value[0])
            e2(i7.value, decision_1.decided_value[1])

            # Compute performances
            # (rho, w, h, L, n_struts, cost_coeff)
            b2(i5.value, d1.value, d2.value, b1.intermediate, decision_1.selection_value, i6.value)
            p1(b2.performance[0])
            p2(b2.performance[1])


    man = MAN(design_params, input_specs, fixed_params,
            behaviours, decisions, margin_nodes, performances, 'MAN_1')

    # Create surrogate model sfor estimating threshold performance
    # man.train_performance_surrogate(n_samples=1000, sm_type='LS', sampling_freq=1, num_threads=num_threads)
    # man.save('strut_fem')
    # # sys.exit(0)

    # # Train surrogate for B3
    # variable_dict = {
    #     'n_struts' : {'type' : 'INT', 'limits' : decision_1.universe},
    #     'w' : {'type' : 'FLOAT', 'limits' : d1.universe},
    #     'h' : {'type' : 'FLOAT', 'limits' : d2.universe},
    #     'theta' : {'type' : 'FLOAT', 'limits' : d3.universe},
    #     'BX' : {'type' : 'FLOAT', 'limits' : s1.universe},
    #     'BY' : {'type' : 'FLOAT', 'limits' : s2.universe},
    #     'E' : {'type' : 'fixed', 'limits' : i1.value},
    #     'r1' : {'type' : 'fixed', 'limits' : i2.value},
    #     'r2' : {'type' : 'fixed', 'limits' : i3.value},
    # }

    # b3.train_surrogate(variable_dict,n_outputs=2,n_samples=1000,num_threads=num_threads)
    # man.save('strut_fem')
    # sys.exit(0)

    # load the MAN
    man.load('strut_fem')
    # man.train_performance_surrogate(ext_samples=(man.xt,man.yt), sm_type='LS')

    man.init_decisions(num_threads=num_threads)
    man.forward()
    man.view_perf(e_indices=[0, 1], p_index=0)
    man.view_perf(e_indices=[0, 1], p_index=1)

    # Perform Monte-Carlo simulation
    n_epochs = 100
    for n in range(n_epochs):
        sys.stdout.write("Progress: %d%%   \r" % ((n / n_epochs) * 100))
        sys.stdout.flush()

        man.randomize()
        man.init_decisions(num_threads=num_threads)
        man.forward()
        man.compute_impact()
        man.compute_absorption(num_threads=num_threads)

    man.save('strut_fem_100')

    # View distribution of excess
    man.margin_nodes[0].excess.view(xlabel='E1')
    man.margin_nodes[1].excess.view(xlabel='E2')

    # View distribution of Impact on Performance
    man.impact_matrix.view(0, 0, xlabel='E1,P1')
    man.impact_matrix.view(1, 0, xlabel='E2,P1')
    man.impact_matrix.view(0, 1, xlabel='E1,P2')
    man.impact_matrix.view(1, 1, xlabel='E2,P2')

    man.deterioration_vector.view(0, xlabel='S1')
    man.deterioration_vector.view(1, xlabel='S2')

    man.absorption_matrix.view(0, 0, xlabel='E1,S1')
    man.absorption_matrix.view(1, 0, xlabel='E2,S1')

    man.absorption_matrix.view(0, 1, xlabel='E1,S2')
    man.absorption_matrix.view(1, 1, xlabel='E2,S2')

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
