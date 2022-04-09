import sys

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from mvm import Design
from mvm import FixedParam, DesignParam, InputSpec, Behaviour, Performance, MarginNode, MarginNetwork, Decision
from mvm import GaussianFunc, UniformFunc

# define fixed parameters
i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
i5 = FixedParam(1.0, 'I5', description='Column effective length factor', symbol='K')
i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')

fixed_params = [i1, i2, i3, i4, i5, i6]

# define design parameters
d1 = DesignParam(100.0, 'D1', universe=[70.0, 130.0], variable_type='FLOAT', description='vane length', symbol='w')
d2 = DesignParam(15.0, 'D2', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
d3 = DesignParam(10.0, 'D3', universe=[0.0, 30.0], variable_type='FLOAT', description='lean angle', symbol='theta')
design_params = [d1, d2, d3]

# # T1,T2 distribution (Gaussian)
# center = np.array([450, 425])
# Sigma = np.array([
#     [100, 25],
#     [75, 100],
# ]) / (20 * 0.1)
# Requirement = GaussianFunc(center, Sigma, 'temp')

# T1,T2 distribution (Uniform)
center = np.array([450, 425])
Range = np.array([100, 100]) / (20 * 0.25)
Requirement = UniformFunc(center, Range, 'temp')

# define input specifications
s1 = InputSpec(center[0], 'S1', universe=[325, 550], variable_type='FLOAT', cov_index=0,
               description='nacelle temperature', distribution=Requirement,
               symbol='T1', inc=-1e-0, inc_type='rel')
s2 = InputSpec(center[1], 'S2', universe=[325, 550], variable_type='FLOAT', cov_index=1,
               description='gas surface temperature', distribution=Requirement,
               symbol='T2', inc=+1e-0, inc_type='rel')
input_specs = [s1, s2]


# define the behaviour models

# this is the length model
class B1(Behaviour):
    def __call__(self, theta, r1, r2):
        length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)
        self.intermediate = length


# this is the weight model
class B2(Behaviour):
    def __call__(self, rho, w, h, L, cost_coeff):
        weight = rho * w * h * L
        cost = weight * cost_coeff
        self.performance = [weight, cost]


coeffs = [0.95, 1.05, 0.97]
coeffs = 3 * [1.0, ]


# this is the axial stress model
class B3(Behaviour):
    def __call__(self, alpha, E, r1, r2, Ts, T1, T2, w, h, theta, L):
        force = (E * w * h * alpha) * ((T2 * coeffs[0] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(
            np.deg2rad(theta)) / L
        sigma_a = (E * alpha) * ((T2 * coeffs[1] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(np.deg2rad(theta)) / L
        self.threshold = [force / 1000, sigma_a]


# this is the bending stress model
class B4(Behaviour):
    def __call__(self, alpha, E, r1, r2, Ts, T1, T2, h, theta, L):
        sigma_m = (3 / 2) * ((E * h) / (L ** 2)) * (
                alpha * ((T2 * coeffs[2] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.sin(np.deg2rad(theta)))
        self.threshold = sigma_m


# this is the buckling model
class B5(Behaviour):
    def __call__(self, E, K, w, h, L):
        f_buckling = ((np.pi ** 2) * E * w * (h ** 3)) / (12 * ((K * L) ** 2))
        self.decided_value = f_buckling / 1000


b1 = B1('B1')
b2 = B2('B2')
b3 = B3('B3')
b4 = B4('B4')
b5 = B5('B5')


# class B6(Behaviour):
#     def __call__(self, material):

#         material_dict = {
#             'Inconel'   : {
#                 'sigma_y' : 460, # MPa
#                 'rho' : 8.19e-06, # kg/mm3
#                 'cost' : 0.46 # USD/kg
#                 },
#             'Titanium'  : {
#                 'sigma_y' : 828, # MPa
#                 'rho' : 4.43e-06, # kg/mm3
#                 'cost' : 1.10 # USD/kg
#                 },
#         }

#         chosen_mat = material_dict[material]

#         self.intermediate = [chosen_mat['rho'], chosen_mat['cost']]
#         self.decided_value = chosen_mat['sigma_y']

#         return self.decided_value


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


b6 = B6('B6')

# Define decision nodes and a model to convert to decided values
# decision_1 = Decision(universe=['Inconel', 'Titanium'], variable_type='ENUM', key='decision_1',
#                       direction='must_not_exceed', decided_value_model=b6, description='The type of material')
# decision_1 = Decision(universe=list(range(11)), variable_type='ENUM', key='decision_1',
#                       direction='must_not_exceed', decided_value_model=b6, 
#                       description='The type of material')
decision_1 = Decision(universe=list(range(11)), variable_type='ENUM', key='decision_1',
                      direction=['must_not_exceed','must_not_exceed'], decided_value_model=b6, 
                      n_nodes=2, description='The type of material')

decisions = [decision_1, ]
behaviours = [b1, b2, b3, b4, b5, b6]

# Define margin nodes
e1 = MarginNode('E1', direction='must_not_exceed')
e2 = MarginNode('E2', direction='must_not_exceed')
e3 = MarginNode('E3', direction='must_not_exceed')
margin_nodes = [e1, e2, e3]

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

    def forward(self, num_threads=1, recalculate_decisions=False, override_decisions=False):
        # retrieve MAN components
        d1 = self.design_params[0]  # w
        d2 = self.design_params[1]  # h
        d3 = self.design_params[2]  # theta

        s1 = self.input_specs[0]  # T1 (stochastic)
        s2 = self.input_specs[1]  # T2 (stochastic)

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

        decision_1 = self.decisions[0]  # select a material based on maximum bending and axial stress

        e1 = self.margin_nodes[0]  # margin against buckling (F,F_buckling)
        e2 = self.margin_nodes[1]  # margin against axial failure (sigma_a,sigma_y)
        e3 = self.margin_nodes[2]  # margin against bending failure (sigma_m,sigma_y)

        p1 = self.performances[0]  # weight
        p2 = self.performances[1]  # cost

        # Execute behaviour models
        b1(d3.value, i3.value, i4.value)
        b3(i1.value, i2.value, i3.value, i4.value, i6.value, s1.value, s2.value, d1.value, d2.value, d3.value,
           b1.intermediate)
        b4(i1.value, i2.value, i3.value, i4.value, i6.value, s1.value, s2.value, d2.value, d3.value, b1.intermediate)
        b5(i2.value, i5.value, d1.value, d2.value, b1.intermediate)

        # Execute decision node and translation model
        decision_1([b3.threshold[1], b4.threshold],override_decisions,recalculate_decisions,num_threads)
        b6(decision_1.selection_value)

        # Compute excesses
        e1(b3.threshold[0], b5.decided_value)
        e2(b3.threshold[1], b6.decided_value[0])
        e3(b4.threshold, b6.decided_value[1])

        # Compute performances
        b2(b6.intermediate[0], d1.value, d2.value, b1.intermediate, b6.intermediate[1])
        p1(b2.performance[0])
        p2(b2.performance[1])


man = MAN(design_params, input_specs, fixed_params,
          behaviours, decisions, margin_nodes, performances, 'MAN_1')

# Create surrogate model for estimating threshold performance
man.train_performance_surrogate(n_samples=700, sm_type='LS', bandwidth=[1e-3, ], sampling_freq=1)
man.save('strut_s')

# # load the MAN
# man.load('strut_s')

man.init_decisions()
man.forward()
man.view_perf(e_indices=[0, 1], p_index=0)
man.view_perf(e_indices=[0, 2], p_index=0)
man.view_perf(e_indices=[1, 2], p_index=0)
man.view_perf(e_indices=[0, 1], p_index=1)
man.view_perf(e_indices=[0, 2], p_index=1)
man.view_perf(e_indices=[1, 2], p_index=1)

# Perform Monte-Carlo simulation
n_epochs = 1000
for n in range(n_epochs):
    sys.stdout.write("Progress: %d%%   \r" % ((n / n_epochs) * 100))
    sys.stdout.flush()

    man.randomize()
    man.init_decisions()
    man.forward()
    man.compute_impact()
    man.compute_absorption()

# View distribution of excess
man.margin_nodes[0].excess.view(xlabel='E1')
man.margin_nodes[1].excess.view(xlabel='E2')
man.margin_nodes[2].excess.view(xlabel='E3')

# View distribution of Impact on Performance
man.impact_matrix.view(0, 0, xlabel='E1,P1')
man.impact_matrix.view(1, 0, xlabel='E2,P1')
man.impact_matrix.view(2, 0, xlabel='E3,P1')
man.impact_matrix.view(0, 1, xlabel='E1,P2')
man.impact_matrix.view(1, 1, xlabel='E2,P2')
man.impact_matrix.view(2, 1, xlabel='E3,P2')

man.deterioration_vector.view(0, xlabel='S1')
man.deterioration_vector.view(1, xlabel='S2')

man.absorption_matrix.view(0, 0, xlabel='E1,S1')
man.absorption_matrix.view(1, 0, xlabel='E2,S1')
man.absorption_matrix.view(2, 0, xlabel='E3,S1')

man.absorption_matrix.view(0, 1, xlabel='E1,S2')
man.absorption_matrix.view(1, 1, xlabel='E2,S2')
man.absorption_matrix.view(2, 1, xlabel='E3,S2')

# display the margin value plot
man.compute_mvp('scatter')

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
        man.init_decisions()
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
