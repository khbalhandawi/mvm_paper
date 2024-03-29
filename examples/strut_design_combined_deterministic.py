import sys

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from random import randint

from mvm import Design
from mvm import FixedParam, DesignParam, InputSpec, Behaviour, Performance, MarginNode, MarginNetwork, Decision
from mvm import GaussianFunc, UniformFunc
from mvm import nearest

num_threads = 1

# define fixed parameters
i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
i5 = FixedParam(2.5, 'I5', description='Column effective length factor', symbol='K')
i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')

fixed_params = [i1, i2, i3, i4, i5, i6,]

# define design parameters
d1 = DesignParam(15.0, 'D1', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
d2 = DesignParam(10.0, 'D2', universe=[0.0, 50.0], variable_type='FLOAT', description='lean angle', symbol='theta')
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

# define the behaviour models
# this is the force model
class B1(Behaviour):
    def __call__(self, T1, T2, h, theta, alpha, E, r1, r2, Ts):
        coeffs = [0.95, 1.05, 1.05]
        # coeffs = 3 * [1.0, ]
        w_nominal = 60.0
        n_struts_nominal = 12

        length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

        force_thermal = (E * w_nominal * h * alpha) * (((T2 ** coeffs[0]) * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(
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
        coeffs = [0.95, 1.01, 1.01]
        # coeffs = 3 * [1.0, ]
        n_struts_nominal = 12

        length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

        sigma_a_thermal = (E * alpha) * (((T2 ** coeffs[1]) * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.cos(np.deg2rad(theta)) / length
        sigma_m_thermal = (3 / 2) * ((E * h) / (length ** 2)) * (
                alpha * (((T2 ** coeffs[2]) * r2) - (T1 * r1) - (Ts * (r2 - r1))) * np.sin(np.deg2rad(theta)))
        
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
            'dummy'   : {
                'sigma_y' : 92, # MPa
                'rho' : 11.95e-06, # kg/mm3
                'cost' : 0.1 # USD/kg
            },
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

        # generate 10 materials by linearly interpolating
        df = pd.DataFrame(columns=material_dict['dummy'].keys(), index=range(15), dtype=float)
        df.iloc[0] = material_dict['dummy']
        df.iloc[-1] = material_dict['Titanium']
        df.interpolate(method='linear',axis=0,inplace=True)
        material_dict_interp = df.transpose().to_dict()
        chosen_mat = material_dict_interp[material]

        # chosen_mat = material_dict[material]

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
decision_1 = Decision(universe=list(range(60,120+1,1)), variable_type='ENUM', key='decision_1',
                    direction='must_not_exceed', decided_value_model=b2, n_nodes=1,
                    description='the vane width selection')

decision_2 = Decision(universe=list(range(15)), variable_type='ENUM', key='decision_2',
                    direction='must_not_exceed', decided_value_model=b4, n_nodes=1, 
                    description='The type of material')

decision_3 = Decision(universe=list(range(6,30+1,1)), variable_type='ENUM', key='decision_3',
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

    def forward(self, num_threads=1,recalculate_decisions=False,allocate_margin=False,strategy='min_excess',outputs=['dv','dv','dv']):

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
        # decision_3(decision_2.output_value, recalculate_decisions, allocate_margin, strategy, num_threads, outputs[2], *args, **kwargs)
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

base_folder = os.path.join('data','strut_fea')

# # train material surrogate
# variable_dict = {
#     'material' : {'type' : 'ENUM', 'limits' : decision_2.universe},
# }
# b4.train_surrogate(variable_dict,n_samples=50,sm_type='KRG')
# b4.train_inverse(sm_type='LS')

# # Train surrogate for bearing load case
# variable_dict = {
#     'n_struts' : {'type' : 'INT', 'limits' : [decision_3.universe[0],decision_3.universe[-1]]},
#     'w' : {'type' : 'FLOAT', 'limits' : [decision_1.universe[0],decision_1.universe[-1]]},
#     'h' : {'type' : 'FLOAT', 'limits' : d1.universe},
#     'theta' : {'type' : 'FLOAT', 'limits' : d2.universe},
#     'BX' : {'type' : 'FLOAT', 'limits' : s3.universe},
#     'BY' : {'type' : 'FLOAT', 'limits' : s4.universe},
#     'E' : {'type' : 'fixed', 'limits' : i2.value},
#     'r1' : {'type' : 'fixed', 'limits' : i3.value},
#     'r2' : {'type' : 'fixed', 'limits' : i4.value},
# }
# b5.train_surrogate(variable_dict,n_samples=1000,num_threads=num_threads)
# man.save('strut_comb',folder=base_folder)
# b5.train_inverse('n_struts',sm_type='KRG', bandwidth=[1e-3])
# man.save('strut_comb',folder=base_folder)
# sys.exit(0)

# load just the behaviour models
# b4.load(os.path.join(base_folder,'strut_comb'))
# b5.load(os.path.join(base_folder,'strut_comb'))
# b5.train_inverse('n_struts', sm_type='KRG', bandwidth=[1e-3])
# man.save('strut_comb_d',folder=base_folder)

# load the MAN
man.load('strut_comb',folder=base_folder)

man.init_decisions(num_threads=num_threads)
man.decision_vector = [60, 10, 20]
man.allocate_margins('manual')
man.forward()

# Run a forward pass of the MAN
man.compute_impact()
man.compute_absorption()

# View value of excess
print(man.excess_vector)

# View Impact on Performance
print(man.impact_matrix.value)

# View Deterioration
print(man.deterioration_vector.value)

# View Absorption capability
print(man.absorption_matrix.value)

# display the margin value plot
d = man.compute_mvp('scatter', show_neutral=True)

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

# Effect of alternative designs
# designs = [[70,10,6+i] for i in range(25)]
designs = [[60+5*i,10,20] for i in range(13)]


# create empty figure
fig, ax = plt.subplots(figsize=(7, 8))
ax.set_xlabel('Impact on performance')
ax.set_ylabel('Change absorption capability')

X = np.empty((0, len(man.margin_nodes)))
Y = np.empty((0, len(man.margin_nodes)))
E = np.empty((0, len(man.margin_nodes)))
TT = np.empty((0, len(man.margin_nodes)))
DV = np.empty((0, len(man.margin_nodes)))
TT_A = np.empty((0, len(man.margin_nodes)))

for i, design in enumerate(designs):
    man.reset()
    man.reset_outputs()

    # Display progress bar
    sys.stdout.write("Progress: %d%%   \r" % ((i / len(designs)) * 100))
    sys.stdout.flush()

    # Perform MAN computations
    man.init_decisions()
    man.decision_vector = design
    man.allocate_margins('manual')
    man.forward()
    man.compute_impact()
    man.compute_absorption()

    # Extract x and y
    x = np.mean(man.impact_matrix.values,
                axis=(1, 2)).ravel()  # average along performance parameters (assumes equal weighting)
    y = np.mean(man.absorption_matrix.values,
                axis=(1, 2)).ravel()  # average along input specs (assumes equal weighting)

    # if not all(np.isnan(y)):
    #     X = np.vstack((X, x))
    #     Y = np.vstack((Y, y))
    X = np.vstack((X, x))
    Y = np.vstack((Y, y))
    E = np.vstack((E, man.excess_vector))
    TT = np.vstack((TT, man.tt_vector))
    DV = np.vstack((DV, man.dv_vector))

    man.forward(outputs=['tt',]*len(man.decisions)) # get absorption threshold
    TT_A = np.vstack((TT_A, man.tt_vector))

    # plot the results
    color = np.random.random((1, 3))
    ax.scatter(x, y, c=color)

plt.show()

widths = [designs[i][0] for i in range(13)]

# create empty figure for CAC
fig, ax = plt.subplots(figsize=(7, 6))
ax.set_xlabel('Decided width w')
ax.set_ylabel('Change absorption capability (E3)')
ax.plot(widths,Y[:,2])
plt.show()

# create empty figure for IoP
fig, ax = plt.subplots(figsize=(7, 6))
ax.set_xlabel('Decided width w')
ax.set_ylabel('Impact on performance (E1)')
ax.plot(widths,X[:,0])
plt.show()


# Calculate distance metric
p1 = np.array([X.min(), Y.min()])
p2 = np.array([X.max(), Y.max()])
p1 = p1 - 0.1 * abs(p2 - p1)
p2 = p2 + 0.1 * abs(p2 - p1)

distances = np.empty(0)
for i, (x, y) in enumerate(zip(X, Y)):

    dist = 0
    for node in range(len(x)):
        s = np.array([x[node], y[node]])
        pn, d = nearest(p1, p2, s)
        dist += d

    distances = np.append(distances, dist)

distances_sorted = distances[np.argsort(-distances)]
X_sorted = X[np.argsort(-distances)]
Y_sorted = Y[np.argsort(-distances)]
D_sorted = list(np.array(designs)[np.argsort(-distances)])

# display distances only for three designs
X_plot = X_sorted[[0,int(len(X)/2),-1], :]
Y_plot = Y_sorted[[0,int(len(X)/2),-1], :]

# create empty figure
colors = ['#FF0000', '#27BE1E', '#0000FF']
fig, ax = plt.subplots(figsize=(7, 8))
ax.set_xlabel('Impact on performance')
ax.set_ylabel('Change absoption capability')
ax.set_xlim(p1[0], p2[0])
ax.set_ylim(p1[1], p2[1])

p = ax.plot([0, 1], [0, 1], transform=ax.transAxes, color='k', linestyle=(5, (10, 5)))

for i, (x, y) in enumerate(zip(X_plot, Y_plot)):

    ax.scatter(x, y, c=colors[i])

    dist = 0
    for node in range(len(x)):
        s = np.array([x[node], y[node]])
        pn, d = nearest(p1, p2, s)
        dist += d

        x_d = [s[0], pn[0]]
        y_d = [s[1], pn[1]]
        ax.plot(x_d, y_d, marker='.', linestyle='--', color=colors[i])

plt.show()

# create empty figure
fig, ax = plt.subplots(figsize=(7, 8))
ax.set_xlabel('Overall distance from the neutral line')
ax.set_ylabel('Frequency')
ax.hist(distances,bins=len(distances))

plt.show()

best_i = np.argmax(distances)
best_design = designs[best_i]

print('---------------------------------')
result = 'The best design:\n'
for value,d_object in zip(best_design,man.design_params):
    result += d_object.symbol + ' = ' + '%.3f'%value + '\n'

result += '---------------------------------'
print(result)

worst_i = np.argmin(distances)
worst_design = designs[worst_i]

print('---------------------------------')
result = 'The worst design:\n'
for value,d_object in zip(worst_design,man.design_params):
    result += d_object.symbol + ' = ' + '%.3f'%value + '\n'

result += '---------------------------------'
print(result)