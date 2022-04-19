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
i5 = FixedParam(1.5, 'I5', description='Column effective length factor', symbol='K')
i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')

fixed_params = [i1, i2, i3, i4, i5, i6]

# define design parameters
d1 = DesignParam(15.0, 'D2', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
d2 = DesignParam(10.0, 'D3', universe=[0.0, 30.0], variable_type='FLOAT', description='lean angle', symbol='theta')
design_params = [d1, d2]

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
               symbol='T1', inc=-1e-1, inc_type='rel')
s2 = InputSpec(center[1], 'S2', universe=[325, 550], variable_type='FLOAT', cov_index=1,
               description='gas surface temperature', distribution=Requirement,
               symbol='T2', inc=+1e-1, inc_type='rel')
input_specs = [s1, s2]


# define the behaviour models
# this is the force model
class B1(Behaviour):
    def __call__(self, T1, T2, h, theta, alpha, E, r1, r2, Ts):
        coeffs = [0.95, 1.05, 0.97]
        coeffs = 3 * [1.0, ]
        w_nominal = 100.0

        length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

        force = (E * w_nominal * h * alpha) * ((T2 * coeffs[0] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(
            np.deg2rad(theta)) / length
        self.threshold = force / 1000


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

        length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

        sigma_a = (E * alpha) * ((T2 * coeffs[1] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(np.deg2rad(theta)) / length
        sigma_m = (3 / 2) * ((E * h) / (length ** 2)) * (
                alpha * ((T2 * coeffs[2] * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.sin(np.deg2rad(theta)))
        self.threshold = max([sigma_a, sigma_m])


# this is the material model
class B4(Behaviour):
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
        df = pd.DataFrame(columns=material_dict['Inconel'].keys(), index=range(15), dtype=float)
        df.iloc[0] = material_dict['Inconel']
        df.iloc[-1] = material_dict['Titanium']
        df.interpolate(method='linear',axis=0,inplace=True)
        material_dict_interp = df.transpose().to_dict()

        chosen_mat = material_dict_interp[material]

        self.intermediate = [chosen_mat['rho'], chosen_mat['cost']]
        self.decided_value = chosen_mat['sigma_y']


# this is the weight model
class B5(Behaviour):
    def __call__(self, rho, cost_coeff, w, h, theta, r1, r2):
        
        length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

        weight = rho * w * h * length
        cost = weight * cost_coeff
        self.performance = [weight, cost]


b1 = B1(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B1')
b2 = B2(n_i=0, n_p=2, n_dv=1, n_tt=0, key='B2')
b3 = B3(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B3')
b4 = B4(n_i=2, n_p=0, n_dv=1, n_tt=0, key='B4')
b5 = B5(n_i=0, n_p=2, n_dv=0, n_tt=0, key='B5')


# Define decision nodes and a model to convert to decided values
decision_1 = Decision(universe=list(range(60,120+5,10)), variable_type='ENUM', key='decision_1',
                    direction='must_not_exceed', decided_value_model=b2, n_nodes=1,
                    description='the vane width selection')

decision_2 = Decision(universe=list(range(15)), variable_type='ENUM', key='decision_2',
                    direction='must_not_exceed', decided_value_model=b4, n_nodes=1, 
                    description='The type of material')

decisions = [decision_1, decision_2]
behaviours = [b1, b2, b3, b4, b5]

# Define margin nodes
e1 = MarginNode('E1', direction='must_not_exceed')
e2 = MarginNode('E2', direction='must_not_exceed')
margin_nodes = [e1, e2]

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

    def forward(self, num_threads=1, recalculate_decisions=False, override_decisions=False,outputs=['dv','dv','dv']):
        # retrieve MAN components
        d1 = self.design_params[0]  # h
        d2 = self.design_params[1]  # theta

        s1 = self.input_specs[0]  # T1 (stochastic)
        s2 = self.input_specs[1]  # T2 (stochastic)

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
        b5 = self.behaviours[4]  # calculates weight and cost

        decision_1 = self.decisions[0]  # select the width of the vane based on the maximum supported buckling load
        decision_2 = self.decisions[1]  # select the number of struts based on center displacement and max stress

        e1 = self.margin_nodes[0]  # margin against buckling (F,F_buckling)
        e2 = self.margin_nodes[1]  # margin against axial or bending failure (max(sigma_a,sigma_m),sigma_y)

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
        decision_1(b1.threshold, override_decisions, recalculate_decisions, num_threads, outputs[0],*args)
        # invert decided value: decided_value, h, theta, E, r1, r2, K
        b2.inv_call(decision_1.output_value, d1.value, d2.value, i2.value, i3.value, i4.value, i5.value)

        # T1, T2, h, theta, alpha, E, r1, r2, Ts)
        b3(s1.value, s2.value, d1.value, d2.value, i1.value, i2.value, i3.value, i4.value, i6.value)
        # Execute decision node for material and translate to yield stress: material
        decision_2(b3.threshold, override_decisions, recalculate_decisions, num_threads, outputs[1])
        # invert decided value: decided_value, h, theta, E, r1, r2, K
        b4.inv_call(decision_2.output_value)

        # Compute excesses
        e1(b1.threshold, decision_1.decided_value)
        e2(b3.threshold, decision_2.decided_value)

        # Compute performances
        # rho, cost_coeff, w, h, theta, r1, r2
        b5(b4.intermediate[0], b4.intermediate[1], b2.inverted, d1.value, d2.value, i3.value, i4.value)
        p1(b5.performance[0])
        p2(b5.performance[1])


man = MAN(design_params, input_specs, fixed_params,
          behaviours, decisions, margin_nodes, performances, 'MAN_1')

# train material surrogate
variable_dict = {
    'material' : {'type' : 'ENUM', 'limits' : decision_2.universe},
}
b4.train_surrogate(variable_dict,n_samples=50,sm_type='KRG')
b4.train_inverse(sm_type='LS')

# Create surrogate model for estimating threshold performance
man.save('strut_s',folder='strut')

man.init_decisions()
man.forward()

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

man.save('strut_s',folder='strut')

# # load the MAN
# man.load('strut_s_100',folder='strut')

# View distribution of excess
man.margin_nodes[0].excess.view(xlabel='E1')
man.margin_nodes[1].excess.view(xlabel='E2')

# View distribution of Impact on Performance
man.impact_matrix.view(0, 0, xlabel='E1,P1')
man.impact_matrix.view(1, 0, xlabel='E2,P1')

man.deterioration_vector.view(0, xlabel='S1')
man.deterioration_vector.view(1, xlabel='S2')

man.absorption_matrix.view(0, 0, xlabel='E1,S1')
man.absorption_matrix.view(1, 0, xlabel='E2,S1')

man.absorption_matrix.view(0, 1, xlabel='E1,S2')
man.absorption_matrix.view(1, 1, xlabel='E2,S2')

impact_matrix = np.nanmean(man.impact_matrix.values,axis=(2,))  # average along performance parameters (assumes equal weighting)
absorption_matrix = np.nanmean(man.absorption_matrix.values,axis=(2,))  # average along input specs (assumes equal weighting)

rows = ['E%i'%i for i in range(absorption_matrix.shape[0])]
cols = ['S%i'%i for i in range(absorption_matrix.shape[1])]
plt.imshow(absorption_matrix, cmap='hot', interpolation='nearest')
plt.xticks(range(0, absorption_matrix.shape[1]),cols)
plt.yticks(range(0, absorption_matrix.shape[0]),rows)

plt.show()

rows = ['E%i'%i for i in range(impact_matrix.shape[0])]
cols = ['S%i'%i for i in range(impact_matrix.shape[1])]
plt.imshow(impact_matrix, cmap='hot', interpolation='nearest')
plt.xticks(range(0, impact_matrix.shape[1]),cols)
plt.yticks(range(0, impact_matrix.shape[0]),rows)

plt.show()

# display the margin value plot
man.compute_mvp('scatter')
# sys.exit(0)

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

X = np.empty((0,len(man.margin_nodes)))
Y = np.empty((0,len(man.margin_nodes)))
D = np.empty((0,len(man.design_params)))

for i,design in enumerate(design_doe.unscale()):
    man.nominal_design_vector = design
    man.reset()
    man.reset_outputs()

    # Display progress bar
    sys.stdout.write("Progress: %d%%   \r" %((i/n_designs)* 100))
    sys.stdout.flush()

    # Perform MAN computations
    man.init_decisions()
    man.forward()
    man.compute_impact()
    man.compute_absorption()

    # Extract x and y
    x = np.mean(man.impact_matrix.values,axis=(1,2)).ravel() # average along performance parameters (assumes equal weighting)
    y = np.mean(man.absorption_matrix.values,axis=(1,2)).ravel() # average along input specs (assumes equal weighting)

    if not all(np.isnan(y)):
        X = np.vstack((X,x))
        Y = np.vstack((Y,y))
        D = np.vstack((D,design))

    # plot the results
    color = np.random.random((1,3))
    ax.scatter(x,y,c=color)

plt.show()