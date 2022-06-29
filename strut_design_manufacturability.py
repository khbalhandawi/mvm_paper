import sys

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from mvm import Design
from mvm import FixedParam, DesignParam, InputSpec, Behaviour, Performance, MarginNode, MarginNetwork, Decision
from mvm import GaussianFunc, UniformFunc
from mvm import nearest
from mvm.utilities import check_folder

# define fixed parameters
i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
i5 = FixedParam(1.5, 'I5', description='Column effective length factor', symbol='K')
i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')

fixed_params = [i1, i2, i3, i4, i5, i6]

# define design parameters
d1 = DesignParam(15.0, 'D1', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
d2 = DesignParam(10.0, 'D2', universe=[0.0, 30.0], variable_type='FLOAT', description='lean angle', symbol='theta')
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
               symbol='T1', inc=-1e-2, inc_type='rel')
s2 = InputSpec(center[1], 'S2', universe=[325, 550], variable_type='FLOAT', cov_index=1,
               description='gas surface temperature', distribution=Requirement,
               symbol='T2', inc=+1e-2, inc_type='rel')
input_specs = [s1, s2]


# define the behaviour models
# this is the force model
class B1(Behaviour):
    def __call__(self, T1, T2, h, theta, alpha, E, r1, r2, Ts):
        coeffs = [0.95, 1.05, 0.97]
        coeffs = 3 * [1.0, ]
        w_nominal = 60.0

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
            'dummy'   : {
                'sigma_y' : 92, # MPa
                'rho' : 11.95e-06, # kg/mm3
                'cost' : -0.18, # USD/kg
                'ED_mat' : 94.6,
                'C_ps' : 0.4248,
                'T_melt' : 1671, # K
                'C_pl' : 0.633
                },
            'Inconel'  : {
                'sigma_y' : 460,  # MPa
                'rho' : 8.19e-06,  # kg/mm3
                'cost' : 0.46,  # USD/kg
                'ED_mat' : 144.8,
                'C_ps' : 0.4889,
                'T_melt' : 1806, # K
                'C_pl' : 0.6915
                },
            'Titanium'  : {
                'sigma_y' : 828, # MPa
                'rho' : 4.43e-06, # kg/mm3
                'cost' : 1.10, # USD/kg
                'ED_mat' : 195,
                'C_ps' : 0.553,
                'T_melt' : 1941, # K
                'C_pl' : 0.75
                },
        }

        chosen_mat = material_dict[material]

        self.intermediate = [chosen_mat['rho'], chosen_mat['cost'], chosen_mat['ED_mat']]
        self.decided_value = chosen_mat['sigma_y']


# this is the weight model
class B5(Behaviour):
    def __call__(self, rho, cost_coeff, w, h, theta, r1, r2):
        
        length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

        weight = rho * w * h * length
        cost = weight * cost_coeff
        self.performance = [weight, cost]

# this is the AM manufacturability model
class B6(Behaviour):
    def __call__(self, rho, ED_mat, w, h, theta, r1, r2):
        
        # fixed parameters
        N_L	                = 383.7588405
        Lt	                = 0.5
        a	                = 0.000556
        C_mach	            = 40
        C_inert	            = 5
        G_cons	            = 100
        C_energy	        = 30
        P_cons	            = 10
        V_cham	            = 26600000
        P_m	                = 30
        hatch	            = 0.1
        beta	            = 0.18
        Laser_P	            = 155
        alpha	            = 5
        Laser_S	            = 300
        E_D	                = 5.093819785

        # calculations
        L                   = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)
        A_tot	            = 2*((w*h)+(L*h)+(L*w))
        A_s	                = 2*L*(h+w)
        V_tot	            = w*h*L
        W_mat	            = rho*w*h*L
        AM_Build_Time	    = N_L*(A_tot+A_s)*a
        K_u	                = V_tot/V_cham
        C_exec	            = AM_Build_Time*(C_mach+(C_inert*G_cons)+(C_energy*P_cons*K_u))
        C_mat	            = W_mat*P_m
        C_manufacturing	    = C_exec+C_mat
        K	                = E_D/ED_mat
        Manufacturability	= C_manufacturing*K

        self.performance = Manufacturability


b1 = B1(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B1')
b2 = B2(n_i=0, n_p=2, n_dv=1, n_tt=0, key='B2')
b3 = B3(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B3')
b4 = B4(n_i=3, n_p=0, n_dv=1, n_tt=0, key='B4')
b5 = B5(n_i=0, n_p=2, n_dv=0, n_tt=0, key='B5')
b6 = B6(n_i=0, n_p=1, n_dv=0, n_tt=0, key='B6')

# Define decision nodes and a model to convert to decided values
decision_1 = Decision(universe=list(range(60,120+5,10)), variable_type='ENUM', key='decision_1',
                    direction='must_not_exceed', decided_value_model=b2, n_nodes=1,
                    description='the vane width selection')

decision_2 = Decision(universe=['dummy', 'Inconel', 'Titanium'], variable_type='ENUM', key='decision_2',
                    direction='must_not_exceed', decided_value_model=b4, n_nodes=1, 
                    description='The type of material')

# Define margin nodes
e1 = MarginNode('E1', direction='must_not_exceed')
e2 = MarginNode('E2', direction='must_not_exceed')


# Define performances
p1 = Performance('P1', direction='less_is_better')
p2 = Performance('P2', direction='less_is_better')
p3 = Performance('P3', direction='less_is_better')


# Define the MAN
class MAN_AM(MarginNetwork):

    def randomize(self):
        Requirement.random()
        s1.random()
        s2.random()

    def forward(self, num_threads=1, recalculate_decisions=False, allocate_margin=False, strategy='min_excess', outputs=['dv','dv','dv']):
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
        b6 = self.behaviours[5]  # calculates manufacturability

        decision_1 = self.decisions[0]  # select the width of the vane based on the maximum supported buckling load
        decision_2 = self.decisions[1]  # select the number of struts based on center displacement and max stress

        e1 = self.margin_nodes[0]  # margin against buckling (F,F_buckling)
        e2 = self.margin_nodes[1]  # margin against axial or bending failure (max(sigma_a,sigma_m),sigma_y)

        p1 = self.performances[0]  # weight
        p2 = self.performances[1]  # cost
        p3 = self.performances[2]  # manufacturability

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
        # invert decided value: sigma_y
        b4.inv_call(decision_2.output_value)

        # Compute excesses
        e1(b1.threshold, decision_1.decided_value)
        e2(b3.threshold, decision_2.decided_value)

        # Compute performances
        # rho, cost_coeff, w, h, theta, r1, r2
        b5(b4.intermediate[0], b4.intermediate[1], b2.inverted, d1.value, d2.value, i3.value, i4.value)
        # rho, ED_mat, w, h, theta, r1, r2
        b6(b4.intermediate[0], b4.intermediate[2], b2.inverted, d1.value, d2.value, i3.value, i4.value)

        p1(b5.performance[0])
        p2(b5.performance[1])
        p3(b6.performance)


decisions = [decision_1, decision_2]
behaviours = [b1, b2, b3, b4, b5, b6]
margin_nodes = [e1, e2]
performances = [p1, p2, p3]

man_AM = MAN_AM(design_params, input_specs, fixed_params,
          behaviours, decisions, margin_nodes, performances, 'MAN_1')

# train material surrogate
variable_dict = {
    'material' : {'type' : 'ENUM', 'limits' : decision_2.universe},
}
b4.train_surrogate(variable_dict,n_samples=50,sm_type='KRG')
b4.train_inverse(sm_type='LS')

# Create surrogate model for estimating threshold performance
check_folder('strut_manufacturability')
man_AM.save('strut_AM',folder='strut_manufacturability')

man_AM.init_decisions()
man_AM.allocate_margins()
man_AM.forward()

# Run a forward pass of the MAN
man_AM.compute_impact()
man_AM.compute_absorption()

# View Impact on Performance
impact_AM = man_AM.impact_matrix.value

man_AM.reset()

################################################################
# Manufacturing by casting


# this is the material model
class B7(Behaviour):
    def __call__(self, material):

        material_dict = {
            'dummy'   : {
                'sigma_y' : 92, # MPa
                'rho' : 11.95e-06, # kg/mm3
                'cost' : -0.18, # USD/kg
                'ED_mat' : 94.6,
                'C_ps' : 0.4248,
                'T_melt' : 1671, # K
                'C_pl' : 0.633
                },
            'Inconel'  : {
                'sigma_y' : 460,  # MPa
                'rho' : 8.19e-06,  # kg/mm3
                'cost' : 0.46,  # USD/kg
                'ED_mat' : 144.8,
                'C_ps' : 0.4889,
                'T_melt' : 1806, # K
                'C_pl' : 0.6915
                },
            'Titanium'  : {
                'sigma_y' : 828, # MPa
                'rho' : 4.43e-06, # kg/mm3
                'cost' : 1.10, # USD/kg
                'ED_mat' : 195,
                'C_ps' : 0.553,
                'T_melt' : 1941, # K
                'C_pl' : 0.75
                },
        }

        chosen_mat = material_dict[material]

        self.intermediate = [chosen_mat['rho'], chosen_mat['cost'], chosen_mat['C_ps'], chosen_mat['T_melt'], chosen_mat['C_pl']]
        self.decided_value = chosen_mat['sigma_y']


# this is the casting manufacturability model
class B8(Behaviour):
    def __call__(self, rho, cost_coeff, C_ps, T_melt, C_pl, w, h, theta, r1, r2):
        
        # fixed parameters
        F_Loss	        = 1.103
        C_unit_mat	    = cost_coeff
        F_Loss_energy	= 1.9314
        T_room	        = 295.2
        T_tap	        = 2088
        C_ac	        = 80
        C_unit_energy	= 30
        C_s	            = 12

        # calculations
        L                   = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)
        W_mat	            = rho*w*h*L
        C_material	        = C_unit_mat*W_mat*F_Loss
        C_energy	        = C_unit_energy*W_mat*F_Loss_energy*(C_ps*(T_melt-T_room)+C_pl*(T_tap-T_melt))
        V_cast 	            = w*h*L
        C_tooling	        = np.exp((0.629*(V_cast/10E9))+(0.048*C_ac)+(0.023*C_s)+0.739)
        Manufacturability	= C_material+C_energy+C_tooling

        self.performance = Manufacturability


b7 = B7(n_i=5, n_p=0, n_dv=1, n_tt=0, key='B7')
b8 = B8(n_i=0, n_p=1, n_dv=0, n_tt=0, key='B8')

# Define the MAN
class MAN_casting(MarginNetwork):

    def randomize(self):
        s1 = self.input_specs[0]  # T1 (stochastic)
        s2 = self.input_specs[1]  # T2 (stochastic)
        s1.random()
        s2.random()

    def forward(self, num_threads=1, recalculate_decisions=False, allocate_margin=False, strategy='min_excess', outputs=['dv','dv','dv']):
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
        b7 = self.behaviours[3]  # convert material index to yield stress, density, and cost
        b5 = self.behaviours[4]  # calculates weight and cost
        b8 = self.behaviours[5]  # calculates manufacturability

        decision_1 = self.decisions[0]  # select the width of the vane based on the maximum supported buckling load
        decision_2 = self.decisions[1]  # select the number of struts based on center displacement and max stress

        e1 = self.margin_nodes[0]  # margin against buckling (F,F_buckling)
        e2 = self.margin_nodes[1]  # margin against axial or bending failure (max(sigma_a,sigma_m),sigma_y)

        p1 = self.performances[0]  # weight
        p2 = self.performances[1]  # cost
        p3 = self.performances[2]  # manufacturability

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
        # invert decided value: sigma_y
        b7.inv_call(decision_2.output_value)

        # Compute excesses
        e1(b1.threshold, decision_1.decided_value)
        e2(b3.threshold, decision_2.decided_value)

        # Compute performances
        # rho, cost_coeff, w, h, theta, r1, r2
        b5(b7.intermediate[0], b7.intermediate[1], b2.inverted, d1.value, d2.value, i3.value, i4.value)
        # rho, cost_coeff, C_ps, T_melt, C_pl, w, h, theta, r1, r2)
        b8(b7.intermediate[0], b7.intermediate[1], b7.intermediate[2], b7.intermediate[3], b7.intermediate[4], 
            b2.inverted, d1.value, d2.value, i3.value, i4.value)

        p1(b5.performance[0])
        p2(b5.performance[1])
        p3(b8.performance)


decisions = [decision_1, decision_2]
behaviours = [b1, b2, b3, b7, b5, b8]
margin_nodes = [e1, e2]
performances = [p1, p2, p3]

man_casting = MAN_casting(design_params, input_specs, fixed_params,
          behaviours, decisions, margin_nodes, performances, 'MAN_2')


# train material surrogate
variable_dict = {
    'material' : {'type' : 'ENUM', 'limits' : decision_2.universe},
}
b7.train_surrogate(variable_dict,n_samples=50,sm_type='KRG')
b7.train_inverse(sm_type='LS')

# Create surrogate model for estimating threshold performance
check_folder('strut_manufacturability')
man_casting.save('strut_casting',folder='strut_manufacturability')

man_casting.init_decisions()
man_casting.allocate_margins()
man_casting.forward()

# Run a forward pass of the MAN
man_casting.compute_impact()
man_casting.compute_absorption()

# View Impact on Performance
impact_casting = man_casting.impact_matrix.value

################################################################
# Plot results

import pandas as pd
import numpy as np
import altair as alt

e1_values = np.hstack((impact_AM[0,:].reshape(-1,1),impact_casting[0,:].reshape(-1,1)))
e2_values = np.hstack((impact_AM[1,:].reshape(-1,1),impact_casting[1,:].reshape(-1,1)))

e1=pd.DataFrame(e1_values,index=["weight","raw mat. cost","manuf."],columns=["AM","Casting"])
e2=pd.DataFrame(e2_values,index=["weight","raw mat. cost","manuf."],columns=["AM","Casting"])

def prep_df(df, name):
    df = df.stack().reset_index()
    df.columns = ['facet', 'Concept', 'Impact']
    df['Margin'] = name
    return df

df1 = prep_df(e1, 'E1')
df2 = prep_df(e2, 'E2')
df = pd.concat([df1, df2,])

df.to_csv(os.path.join('strut_manufacturability','impact_data.csv'))