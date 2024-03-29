
import numpy as np
import pandas as pd
import os

from mvm import Design
from mvm import FixedParam, DesignParam, InputSpec, Behaviour, Performance, MarginNode, MarginNetwork, Decision
from mvm import GaussianFunc, UniformFunc


def get_man_combined_circ(height=15.0, lean=0.0, materials=None, widths=list(range(60,120+10,10)), 
    train_surrogate=True, man_folder=None, name='strut_comb', overwrite=True, num_threads=1):

    # default material selection
    if materials == None:
        materials = {
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

    # define fixed parameters
    i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
    i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
    i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
    i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
    i5 = FixedParam(1.5, 'I5', description='Column effective length factor', symbol='K')
    i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')
    i7 = FixedParam(201.0, 'I7', description='shroud width', symbol='ws')

    fixed_params = [i1, i2, i3, i4, i5, i6, i7]

    # define design parameters
    d1 = DesignParam(height, 'D1', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
    d2 = DesignParam(lean, 'D2', universe=[0.0, 50.0], variable_type='FLOAT', description='lean angle', symbol='theta')
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
        def __call__(self, w, h, theta, E, r1, r2, K, **kwargs):
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


    class B7(Behaviour):
        def __call__(self, BX, BY):
            offset = -21960.00
            scaling = 177.17

            self.threshold = np.sqrt(BX ** 2 + BY ** 2) * scaling + offset


    # this is the material model
    class B4(Behaviour):
        def __call__(self, material, **kwargs):

            chosen_mat = materials[material]

            self.intermediate = [chosen_mat['rho'], chosen_mat['cost']]
            self.decided_value = chosen_mat['sigma_y']


    # this is the simulation model
    class B5(Behaviour):
        def __call__(self, n_struts, w, h, theta, BX, BY, E, r1, r2, **kwargs):
            # Define arguments needed to build surrogate
            args = [n_struts,w,h,theta,BX,BY] # E, r1, r2 are fixed parameters
            if self.surrogate_available:
                return Behaviour.__call__(self,*args) # once surrogate is trained, the base class method will terminate here

            from fem_scirpts.nx_trs_script import run_nx_simulation
            from fem_scirpts.plot_results import postprocess_vtk, plotmesh, plotresults

            try:
                _, vtk_filename, success = run_nx_simulation(w,h,theta,n_struts,r_hub=r1,r_shroud=r2,
                    youngs_modulus=E,bearing_x=BX,bearing_y=BY,pid=kwargs['id'])
            except:
                vtk_filename = None
                success = False

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
        

    # this is the shroud model
    class B8(Behaviour):
        def __call__(self, ts, sigma_y, n_struts, theta, r1, r2, ws, **kwargs):
            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            B_max = (sigma_y * ws * ts**2) / (6 * (length * np.sin(np.deg2rad(theta)) + (r2 - r1) * np.sin(2 * np.pi / n_struts)))

            self.decided_value = B_max

        def inv_call(self, decided_value, sigma_y, n_struts, theta, r1, r2, ws):
            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            # B = max(decided_value - offset,1e-3)
            B = decided_value
            A = ((6 * B) / (ws * sigma_y))
            ts = np.sqrt(A * (length * np.sin(np.deg2rad(theta)) + (r2 - r1) * np.sin(2 * np.pi / n_struts)))

            self.inverted = ts


    # this is the weight model
    class B6(Behaviour):
        def __call__(self, rho, cost_coeff, w, n_struts, h, theta, r1, r2):

            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            weight = rho * w * h * length * n_struts
            cost = weight * cost_coeff
            self.performance = [weight, cost]


    class B9(Behaviour):
        def __call__(self, ts, rho, cost_coeff, r2, ws):
            
            Ws = rho * np.pi * ts * (ts + 2 * r2) * ws
            Cs = Ws * cost_coeff

            self.performance = [Ws, Cs]

    b1 = B1(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B1')
    b2 = B2(n_i=0, n_p=0, n_dv=1, n_tt=0, key='B2')
    b3 = B3(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B3')
    b4 = B4(n_i=2, n_p=0, n_dv=1, n_tt=0, key='B4')
    b7 = B7(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B7')
    b5 = B5(n_i=0, n_p=0, n_dv=1, n_tt=0, key='B5')
    b6 = B6(n_i=0, n_p=2, n_dv=0, n_tt=0, key='B6')
    b8 = B8(n_i=0, n_p=0, n_dv=1, n_tt=0, key='B8')
    b9 = B9(n_i=0, n_p=2, n_dv=0, n_tt=0, key='B9')

    # Define decision nodes and a model to convert to decided values
    decision_1 = Decision(universe=widths, variable_type='ENUM', key='decision_1',
                        direction='must_not_exceed', decided_value_model=b2, n_nodes=1,
                        description='the vane width selection')

    decision_2 = Decision(universe=list(materials.keys()), variable_type='ENUM', key='decision_2',
                        direction='must_not_exceed', decided_value_model=b4, n_nodes=1, 
                        description='The type of material')

    decision_3 = Decision(universe=list(range(6,30+1,2)), variable_type='ENUM', key='decision_3',
                        direction='must_not_exceed', decided_value_model=None, n_nodes=1,
                        description='the number of struts')

    decision_4 = Decision(universe=list(np.arange(2, 8+0.5, 0.5)), variable_type='ENUM', key='decision_4',
                        direction='must_not_exceed', decided_value_model=b8, n_nodes=1,
                        description='the shroud thickness')

    decisions = [decision_1, decision_2, decision_3, decision_4]
    behaviours = [b1, b2, b3, b4, b5, b6, b7, b8, b9]

    # Define margin nodes
    e1 = MarginNode('E1', direction='must_not_exceed')
    e2 = MarginNode('E2', direction='must_not_exceed')
    e3 = MarginNode('E3', direction='must_not_exceed')
    e4 = MarginNode('E4', direction='must_not_exceed')
    margin_nodes = [e1, e2, e3, e4]

    # Define performances
    p1 = Performance('P1', direction='less_is_better')
    p2 = Performance('P2', direction='less_is_better')
    p3 = Performance('P3', direction='less_is_better')
    p4 = Performance('P4', direction='less_is_better')
    performances = [p1, p2, p3, p4]

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

        def forward(self,num_threads=1,recalculate_decisions=False,allocate_margin=False,strategy=['min_excess',]*len(margin_nodes),outputs=['dv',]*len(margin_nodes)):

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
            i7 = self.fixed_params[6]  # ws

            b1 = self.behaviours[0]  # calculates axial force
            b2 = self.behaviours[1]  # calculates buckling load
            b3 = self.behaviours[2]  # calculates bending and axial stresses
            b4 = self.behaviours[3]  # convert material index to yield stress, density, and cost
            b5 = self.behaviours[4]  # calculates von mises stress from a simulation
            b6 = self.behaviours[5]  # calculates weight and cost
            b7 = self.behaviours[6]  # calculates shroud load
            b8 = self.behaviours[7]  # calculates shroud load carrying capacity
            b9 = self.behaviours[8]  # calculates shroud weight and cost

            decision_1 = self.decisions[0]  # select the width of the vane based on the maximum supported buckling load
            decision_2 = self.decisions[1]  # select the material based on the max stress
            decision_3 = self.decisions[2]  # select the number of struts based on center displacement and max stress
            decision_4 = self.decisions[3]  # select the thickness of the shroud based on bending stress

            e1 = self.margin_nodes[0]  # margin against buckling (F,F_buckling)
            e2 = self.margin_nodes[1]  # margin against axial or bending failure (max(sigma_a,sigma_m),sigma_y)
            e3 = self.margin_nodes[2]  # margin against yielding due to bearing loads(sigma,sigma_y)
            e4 = self.margin_nodes[3]  # margin against yielding in the shroud due to bearing loads(sigma_s,sigma_y)

            p1 = self.performances[0]  # weight
            p2 = self.performances[1]  # cost
            p3 = self.performances[2]  # shroud weight
            p4 = self.performances[3]  # shroud cost

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
            decision_1(b1.threshold, recalculate_decisions, allocate_margin, strategy[0], num_threads, outputs[0], *args)
            # invert decided value: decided_value, h, theta, E, r1, r2, K
            b2.inv_call(decision_1.output_value, d1.value, d2.value, i2.value, i3.value, i4.value, i5.value)

            # T1, T2, h, theta, alpha, E, r1, r2, Ts)
            b3(s1.value, s2.value, d1.value, d2.value, i1.value, i2.value, i3.value, i4.value, i6.value)
            # Execute decision node for material and translate to yield stress: material
            decision_2(b3.threshold, recalculate_decisions, allocate_margin, strategy[1], num_threads, outputs[1])
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
            # decision_3(decision_2.output_value, recalculate_decisions, allocate_margin, strategy[2], num_threads, outputs[2], *args, **kwargs)
            # # invert decided value: decided_value, w, h, theta, BX, BY, E, r1, r2
            # b5.inv_call(decision_3.output_value, b2.inverted, d1.value, d2.value, s3.value, s4.value, i2.value, i3.value, i4.value)

            # calculate n_struts: sigma_y, w, h, theta, BX, BY, E, r1, r2
            b5.inv_call(decision_2.output_value, b2.inverted, d1.value, d2.value, s3.value, s4.value, i2.value, i3.value, i4.value)
            # Execute decision node for struts: n_struts, w, h, theta, BX, BY, E, r1, r2
            decision_3(b5.inverted, recalculate_decisions, allocate_margin, strategy[2], num_threads, outputs[2])

            # BX, BY
            b7(s3.value, s4.value)
            # Execute decision node for shroud thickness: ts, n_struts, theta, BX, BY, r1, r2, ws
            args = [
                decision_2.output_value, # sigma_y
                b5.inverted, # n_struts
                self.design_params[1].value, # theta
                self.fixed_params[2].value, # r1
                self.fixed_params[3].value, # r2
                self.fixed_params[6].value, # ws
            ]
            decision_4(b7.threshold, recalculate_decisions, allocate_margin, strategy[3], num_threads, outputs[3], *args)
            # invert decided value: decided_value, sigma_y, n_struts, theta, r1, r2, ws
            b8.inv_call(decision_4.output_value, decision_2.output_value, b5.inverted, d2.value, i3.value, i4.value, i7.value)

            # Compute excesses
            e1(b1.threshold, decision_1.decided_value)
            e2(b3.threshold, decision_2.decided_value)
            e3(b5.inverted, decision_3.decided_value)
            e4(b7.threshold, decision_4.decided_value)

            # Compute performances
            # rho, cost_coeff, w, n_struts, h, theta, r1, r2
            b6(b4.intermediate[0], b4.intermediate[1], b2.inverted, decision_3.output_value, d1.value, d2.value, i3.value, i4.value)
            # ts, rho, cost_coeff, r2, ws
            b9(b8.inverted, b4.intermediate[0], b4.intermediate[1], i4.value, i7.value)

            p1(b6.performance[0])
            p2(b6.performance[1])
            p3(b9.performance[0])
            p4(b9.performance[1])


    man = MAN(design_params, input_specs, fixed_params,
            behaviours, decisions, margin_nodes, performances, name)

    if man_folder is None:
        man_folder = '/'

    if train_surrogate:
        # train material surrogate
        variable_dict = {
            'material' : {'type' : 'ENUM', 'limits' : decision_2.universe},
        }
        b4.train_surrogate(variable_dict,n_samples=50,sm_type='KRG')
        b4.train_inverse(sm_type='LS')

        # Train surrogate for bearing load case
        variable_dict = {
            'n_struts' : {'type' : 'INT', 'limits' : [decision_3.universe[0],decision_3.universe[-1]]},
            'w' : {'type' : 'FLOAT', 'limits' : [decision_1.universe[0],decision_1.universe[-1]]},
            'h' : {'type' : 'FLOAT', 'limits' : d1.universe},
            'theta' : {'type' : 'FLOAT', 'limits' : d2.universe},
            'BX' : {'type' : 'FLOAT', 'limits' : s3.universe},
            'BY' : {'type' : 'FLOAT', 'limits' : s4.universe},
            'E' : {'type' : 'fixed', 'limits' : i2.value},
            'r1' : {'type' : 'fixed', 'limits' : i3.value},
            'r2' : {'type' : 'fixed', 'limits' : i4.value},
        }
        b5.train_surrogate(variable_dict,n_samples=1000,num_threads=num_threads)
        man.save(name,folder=man_folder)
        b5.train_inverse('n_struts',sm_type='LS', bandwidth=[1e-3])
        man.save(name,folder=man_folder)
    else:
        # load just the behaviour models
        b4.load(os.path.join(man_folder,name))
        b5.load(os.path.join(man_folder,name))
        b5.train_inverse('n_struts', sm_type='LS', bandwidth=[1e-3])
    
    if overwrite:
        man.save(name,folder=man_folder)

    return man

def get_man_combined_poly(height=15.0, lean=0.0, materials=None, widths=list(range(60,120+10,10)), 
    train_surrogate=True, man_folder=None, name='strut_comb', overwrite=True, num_threads=1):

    # default material selection
    if materials == None:
        materials = {
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

    # define fixed parameters
    i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
    i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
    i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
    i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
    i5 = FixedParam(1.5, 'I5', description='Column effective length factor', symbol='K')
    i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')

    fixed_params = [i1, i2, i3, i4, i5, i6,]

    # define design parameters
    d1 = DesignParam(height, 'D1', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
    d2 = DesignParam(lean, 'D2', universe=[0.0, 50.0], variable_type='FLOAT', description='lean angle', symbol='theta')
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

            chosen_mat = materials[material]

            self.intermediate = [chosen_mat['rho'], chosen_mat['cost']]
            self.decided_value = chosen_mat['sigma_y']


    # this is the simulation model
    class B5(Behaviour):
        def __call__(self, n_struts, w, h, theta, BX, BY, E, r1, r2, **kwargs):
            # Define arguments needed to build surrogate
            args = [n_struts,w,h,theta,BX,BY] # E, r1, r2 are fixed parameters
            if self.surrogate_available:
                return Behaviour.__call__(self,*args) # once surrogate is trained, the base class method will terminate here

            from fem_scirpts.nx_trs_script import run_nx_simulation
            from fem_scirpts.plot_results import postprocess_vtk, plotmesh, plotresults

            try:
                _, vtk_filename, success = run_nx_simulation(w,h,theta,n_struts,r_hub=r1,r_shroud=r2,
                    youngs_modulus=E,bearing_x=BX,bearing_y=BY,pid=kwargs['id'])
            except:
                vtk_filename = None
                success = False

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

    decision_2 = Decision(universe=list(materials.keys()), variable_type='ENUM', key='decision_2',
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

        def forward(self,num_threads=1,recalculate_decisions=False,allocate_margin=False,strategy=['min_excess',]*len(margin_nodes),outputs=['dv',]*len(margin_nodes)):

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
            decision_1(b1.threshold, recalculate_decisions, allocate_margin, strategy[0], num_threads, outputs[0],*args)
            # invert decided value: decided_value, h, theta, E, r1, r2, K
            b2.inv_call(decision_1.output_value, d1.value, d2.value, i2.value, i3.value, i4.value, i5.value)

            # T1, T2, h, theta, alpha, E, r1, r2, Ts)
            b3(s1.value, s2.value, d1.value, d2.value, i1.value, i2.value, i3.value, i4.value, i6.value)
            # Execute decision node for material and translate to yield stress: material
            decision_2(b3.threshold, recalculate_decisions, allocate_margin, strategy[1], num_threads, outputs[1])
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
            # decision_3(decision_2.output_value, recalculate_decisions, allocate_margin, strategy[2], num_threads, outputs[2], *args, **kwargs)
            # # invert decided value: decided_value, w, h, theta, BX, BY, E, r1, r2
            # b5.inv_call(decision_3.output_value, b2.inverted, d1.value, d2.value, s3.value, s4.value, i2.value, i3.value, i4.value)

            # calculate n_struts sigma_y, w, h, theta, BX, BY, E, r1, r2
            b5.inv_call(decision_2.output_value, b2.inverted, d1.value, d2.value, s3.value, s4.value, i2.value, i3.value, i4.value)
            # Execute decision node for struts: n_struts, w, h, theta, BX, BY, E, r1, r2
            decision_3(b5.inverted, recalculate_decisions, allocate_margin, strategy[2], num_threads, outputs[2])

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
            behaviours, decisions, margin_nodes, performances, name)

    if man_folder is None:
        man_folder = '/'

    if train_surrogate:
        # train material surrogate
        variable_dict = {
            'material' : {'type' : 'ENUM', 'limits' : decision_2.universe},
        }
        b4.train_surrogate(variable_dict,n_samples=50,sm_type='KRG')
        b4.train_inverse(sm_type='LS')

        # Train surrogate for bearing load case
        variable_dict = {
            'n_struts' : {'type' : 'INT', 'limits' : [decision_3.universe[0],decision_3.universe[-1]]},
            'w' : {'type' : 'FLOAT', 'limits' : [decision_1.universe[0],decision_1.universe[-1]]},
            'h' : {'type' : 'FLOAT', 'limits' : d1.universe},
            'theta' : {'type' : 'FLOAT', 'limits' : d2.universe},
            'BX' : {'type' : 'FLOAT', 'limits' : s3.universe},
            'BY' : {'type' : 'FLOAT', 'limits' : s4.universe},
            'E' : {'type' : 'fixed', 'limits' : i2.value},
            'r1' : {'type' : 'fixed', 'limits' : i3.value},
            'r2' : {'type' : 'fixed', 'limits' : i4.value},
        }
        b5.train_surrogate(variable_dict,n_samples=1000,num_threads=num_threads)
        man.save(name,folder=man_folder)
        b5.train_inverse('n_struts',sm_type='LS', bandwidth=[1e-3])
        man.save(name,folder=man_folder)
    else:
        # load just the behaviour models
        b4.load(os.path.join(man_folder,name))
        b5.load(os.path.join(man_folder,name))
        b5.train_inverse('n_struts', sm_type='LS', bandwidth=[1e-3])
    
    if overwrite:
        man.save(name,folder=man_folder)

    return man

def get_man(height=15.0, lean=0.0, materials=None, widths=list(range(60,120+5,10)), 
    train_surrogate=True, man_folder=None, overwrite=True, name='strut_s', num_threads=1):

    # default material selection
    if materials == None:
        materials = {
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

    # define fixed parameters
    i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
    i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
    i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
    i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
    i5 = FixedParam(1.5, 'I5', description='Column effective length factor', symbol='K')
    i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')

    fixed_params = [i1, i2, i3, i4, i5, i6]

    # define design parameters
    d1 = DesignParam(height, 'D2', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
    d2 = DesignParam(lean, 'D3', universe=[10.0, 50.0], variable_type='FLOAT', description='lean angle', symbol='theta')
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
        def __call__(self, w, h, theta, E, r1, r2, K, **kwargs):
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
        def __call__(self, material, **kwargs):

            chosen_mat = materials[material]

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
    decision_1 = Decision(universe=widths, variable_type='ENUM', key='decision_1',
                        direction='must_not_exceed', decided_value_model=b2, n_nodes=1,
                        description='the vane width selection')

    decision_2 = Decision(universe=list(materials.keys()), variable_type='ENUM', key='decision_2',
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
            s1 = self.input_specs[0]  # T1 (stochastic)
            s2 = self.input_specs[1]  # T2 (stochastic)
            s1.random()
            s2.random()

        def forward(self,num_threads=1,recalculate_decisions=False,allocate_margin=False,strategy=['min_excess',]*len(margin_nodes),outputs=['dv',]*len(margin_nodes)):
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
            decision_1(b1.threshold, recalculate_decisions, allocate_margin, strategy[0], num_threads, outputs[0],*args)
            # invert decided value: decided_value, h, theta, E, r1, r2, K
            b2.inv_call(decision_1.output_value, d1.value, d2.value, i2.value, i3.value, i4.value, i5.value)

            # T1, T2, h, theta, alpha, E, r1, r2, Ts)
            b3(s1.value, s2.value, d1.value, d2.value, i1.value, i2.value, i3.value, i4.value, i6.value)
            # Execute decision node for material and translate to yield stress: material
            decision_2(b3.threshold, recalculate_decisions, allocate_margin, strategy[1], num_threads, outputs[1])
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
            behaviours, decisions, margin_nodes, performances, name)

    if train_surrogate:
        # train material surrogate
        variable_dict = {
            'material' : {'type' : 'ENUM', 'limits' : decision_2.universe},
        }
        b4.train_surrogate(variable_dict,n_samples=50,sm_type='KRG')
        b4.train_inverse(sm_type='LS')
    else:
        # load just the behaviour models
        b4.load(os.path.join(man_folder,name))


    # Create surrogate model for estimating threshold performance
    if overwrite:
        man.save(name,folder=man_folder)

    return man

def get_man_combined_circ_old(height=15.0, lean=0.0, materials=None, widths=list(range(60,120+10,10)), 
    train_surrogate=True, man_folder=None, name='strut_comb', overwrite=True, num_threads=1):

    # default material selection
    if materials == None:
        materials = {
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

    # define fixed parameters
    i1 = FixedParam(7.17E-06, 'I1', description='Coefficient of thermal expansion', symbol='alpha')
    i2 = FixedParam(156.3E3, 'I2', description='Youngs modulus', symbol='E')
    i3 = FixedParam(346.5, 'I3', description='Radius of the hub', symbol='r1')
    i4 = FixedParam(536.5, 'I4', description='Radius of the shroud', symbol='r2')
    i5 = FixedParam(1.5, 'I5', description='Column effective length factor', symbol='K')
    i6 = FixedParam(25.0, 'I6', description='ambient temperature', symbol='Ts')
    i7 = FixedParam(201.0, 'I7', description='shroud width', symbol='ws')

    fixed_params = [i1, i2, i3, i4, i5, i6, i7]

    # define design parameters
    d1 = DesignParam(height, 'D1', universe=[5.0, 20.0], variable_type='FLOAT', description='vane height', symbol='h')
    d2 = DesignParam(lean, 'D2', universe=[0.0, 50.0], variable_type='FLOAT', description='lean angle', symbol='theta')
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
        def __call__(self, w, h, theta, E, r1, r2, K, **kwargs):
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


    class B7(Behaviour):
        def __call__(self, BX, BY):
            offset = -21960.00
            scaling = 177.17

            self.threshold = np.sqrt(BX ** 2 + BY ** 2) * scaling + offset


    # this is the material model
    class B4(Behaviour):
        def __call__(self, material, **kwargs):

            chosen_mat = materials[material]

            self.intermediate = [chosen_mat['rho'], chosen_mat['cost']]
            self.decided_value = chosen_mat['sigma_y']


    # this is the simulation model
    class B5(Behaviour):
        def __call__(self, n_struts, w, h, theta, BX, BY, E, r1, r2, **kwargs):
            # Define arguments needed to build surrogate
            args = [n_struts,w,h,theta,BX,BY] # E, r1, r2 are fixed parameters
            if self.surrogate_available:
                return Behaviour.__call__(self,*args) # once surrogate is trained, the base class method will terminate here

            from fem_scirpts.nx_trs_script import run_nx_simulation
            from fem_scirpts.plot_results import postprocess_vtk, plotmesh, plotresults

            try:
                _, vtk_filename, success = run_nx_simulation(w,h,theta,n_struts,r_hub=r1,r_shroud=r2,
                    youngs_modulus=E,bearing_x=BX,bearing_y=BY,pid=kwargs['id'])
            except:
                vtk_filename = None
                success = False

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
        

    # this is the shroud model
    class B8(Behaviour):
        def __call__(self, ts, sigma_y, n_struts, theta, r1, r2, ws, **kwargs):
            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            B_max = (sigma_y * ws * ts**2) / (6 * (length * np.sin(np.deg2rad(theta)) + (r2 - r1) * np.sin(2 * np.pi / n_struts)))

            self.decided_value = B_max

        def inv_call(self, decided_value, sigma_y, n_struts, theta, r1, r2, ws):
            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            # B = max(decided_value - offset,1e-3)
            B = decided_value
            A = ((6 * B) / (ws * sigma_y))
            ts = np.sqrt(A * (length * np.sin(np.deg2rad(theta)) + (r2 - r1) * np.sin(2 * np.pi / n_struts)))

            self.inverted = ts


    # this is the weight model
    class B6(Behaviour):
        def __call__(self, rho, cost_coeff, w, n_struts, h, theta, r1, r2):

            length = -r1 * np.cos(np.deg2rad(theta)) + np.sqrt(r2 ** 2 - (r1 * np.sin(np.deg2rad(theta))) ** 2)

            weight = rho * w * h * length * n_struts
            cost = weight * cost_coeff
            self.performance = [weight, cost]

    # this is the shroud weight model
    class B9(Behaviour):
        def __call__(self, ts, rho, r2, ws):
            
            Ws = rho * np.pi * ts * (ts + 2 * r2) * ws

            self.performance = Ws
    # class B9(Behaviour):
    #     def __call__(self, ts, rho, cost_coeff, r2, ws):
            
    #         Ws = rho * np.pi * ts * (ts + 2 * r2) * ws
    #         Cs = Ws * cost_coeff

    #         self.performance = [Ws, Cs]

    b1 = B1(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B1')
    b2 = B2(n_i=0, n_p=0, n_dv=1, n_tt=0, key='B2')
    b3 = B3(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B3')
    b4 = B4(n_i=2, n_p=0, n_dv=1, n_tt=0, key='B4')
    b7 = B7(n_i=0, n_p=0, n_dv=0, n_tt=1, key='B7')
    b5 = B5(n_i=0, n_p=0, n_dv=1, n_tt=0, key='B5')
    b6 = B6(n_i=0, n_p=2, n_dv=0, n_tt=0, key='B6')
    b8 = B8(n_i=0, n_p=0, n_dv=1, n_tt=0, key='B8')
    b9 = B9(n_i=0, n_p=1, n_dv=0, n_tt=0, key='B9')

    # Define decision nodes and a model to convert to decided values
    decision_1 = Decision(universe=widths, variable_type='ENUM', key='decision_1',
                        direction='must_not_exceed', decided_value_model=b2, n_nodes=1,
                        description='the vane width selection')

    decision_2 = Decision(universe=list(materials.keys()), variable_type='ENUM', key='decision_2',
                        direction='must_not_exceed', decided_value_model=b4, n_nodes=1, 
                        description='The type of material')

    decision_3 = Decision(universe=list(range(6,30+1,2)), variable_type='ENUM', key='decision_3',
                        direction='must_not_exceed', decided_value_model=None, n_nodes=1,
                        description='the number of struts')

    decision_4 = Decision(universe=list(np.arange(2, 8+0.5, 0.5)), variable_type='ENUM', key='decision_4',
                        direction='must_not_exceed', decided_value_model=b8, n_nodes=1,
                        description='the shroud thickness')

    decisions = [decision_1, decision_2, decision_3, decision_4]
    behaviours = [b1, b2, b3, b4, b5, b6, b7, b8, b9]

    # Define margin nodes
    e1 = MarginNode('E1', direction='must_not_exceed')
    e2 = MarginNode('E2', direction='must_not_exceed')
    e3 = MarginNode('E3', direction='must_not_exceed')
    e4 = MarginNode('E4', direction='must_not_exceed')
    margin_nodes = [e1, e2, e3, e4]

    # Define performances
    p1 = Performance('P1', direction='less_is_better')
    p2 = Performance('P2', direction='less_is_better')
    p3 = Performance('P3', direction='less_is_better')
    p4 = Performance('P4', direction='less_is_better')
    performances = [p1, p2, p3]

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

        def forward(self,num_threads=1,recalculate_decisions=False,allocate_margin=False,strategy=['min_excess',]*len(margin_nodes),outputs=['dv',]*len(margin_nodes)):

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
            i7 = self.fixed_params[6]  # ws

            b1 = self.behaviours[0]  # calculates axial force
            b2 = self.behaviours[1]  # calculates buckling load
            b3 = self.behaviours[2]  # calculates bending and axial stresses
            b4 = self.behaviours[3]  # convert material index to yield stress, density, and cost
            b5 = self.behaviours[4]  # calculates von mises stress from a simulation
            b6 = self.behaviours[5]  # calculates weight and cost
            b7 = self.behaviours[6]  # calculates shroud load
            b8 = self.behaviours[7]  # calculates shroud load carrying capacity
            b9 = self.behaviours[8]  # calculates shroud weight and cost

            decision_1 = self.decisions[0]  # select the width of the vane based on the maximum supported buckling load
            decision_2 = self.decisions[1]  # select the material based on the max stress
            decision_3 = self.decisions[2]  # select the number of struts based on center displacement and max stress
            decision_4 = self.decisions[3]  # select the thickness of the shroud based on bending stress

            e1 = self.margin_nodes[0]  # margin against buckling (F,F_buckling)
            e2 = self.margin_nodes[1]  # margin against axial or bending failure (max(sigma_a,sigma_m),sigma_y)
            e3 = self.margin_nodes[2]  # margin against yielding due to bearing loads(sigma,sigma_y)
            e4 = self.margin_nodes[3]  # margin against yielding in the shroud due to bearing loads(sigma_s,sigma_y)

            p1 = self.performances[0]  # weight
            p2 = self.performances[1]  # cost
            p3 = self.performances[2]  # shroud weight
            # p4 = self.performances[3]  # shroud cost

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
            decision_1(b1.threshold, recalculate_decisions, allocate_margin, strategy[0], num_threads, outputs[0], *args)
            # invert decided value: decided_value, h, theta, E, r1, r2, K
            b2.inv_call(decision_1.output_value, d1.value, d2.value, i2.value, i3.value, i4.value, i5.value)

            # T1, T2, h, theta, alpha, E, r1, r2, Ts)
            b3(s1.value, s2.value, d1.value, d2.value, i1.value, i2.value, i3.value, i4.value, i6.value)
            # Execute decision node for material and translate to yield stress: material
            decision_2(b3.threshold, recalculate_decisions, allocate_margin, strategy[1], num_threads, outputs[1])
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
            # decision_3(decision_2.output_value, recalculate_decisions, allocate_margin, strategy[2], num_threads, outputs[2], *args, **kwargs)
            # # invert decided value: decided_value, w, h, theta, BX, BY, E, r1, r2
            # b5.inv_call(decision_3.output_value, b2.inverted, d1.value, d2.value, s3.value, s4.value, i2.value, i3.value, i4.value)

            # calculate n_struts: sigma_y, w, h, theta, BX, BY, E, r1, r2
            b5.inv_call(decision_2.output_value, b2.inverted, d1.value, d2.value, s3.value, s4.value, i2.value, i3.value, i4.value)
            # Execute decision node for struts: n_struts, w, h, theta, BX, BY, E, r1, r2
            decision_3(b5.inverted, recalculate_decisions, allocate_margin, strategy[2], num_threads, outputs[2])

            # BX, BY
            b7(s3.value, s4.value)
            # Execute decision node for shroud thickness: ts, n_struts, theta, BX, BY, r1, r2, ws
            args = [
                decision_2.output_value, # sigma_y
                b5.inverted, # n_struts
                self.design_params[1].value, # theta
                self.fixed_params[2].value, # r1
                self.fixed_params[3].value, # r2
                self.fixed_params[6].value, # ws
            ]
            decision_4(b7.threshold, recalculate_decisions, allocate_margin, strategy[3], num_threads, outputs[3], *args)
            # invert decided value: decided_value, sigma_y, n_struts, theta, r1, r2, ws
            b8.inv_call(decision_4.output_value, decision_2.output_value, b5.inverted, d2.value, i3.value, i4.value, i7.value)

            # Compute excesses
            e1(b1.threshold, decision_1.decided_value)
            e2(b3.threshold, decision_2.decided_value)
            e3(b5.inverted, decision_3.decided_value)
            e4(b7.threshold, decision_4.decided_value)

            # Compute performances
            # rho, cost_coeff, w, n_struts, h, theta, r1, r2
            b6(b4.intermediate[0], b4.intermediate[1], b2.inverted, decision_3.output_value, d1.value, d2.value, i3.value, i4.value)
            # # ts, rho, cost_coeff, r2, ws
            # b9(b8.inverted, b4.intermediate[0], b4.intermediate[1], i4.value, i7.value)
            # ts, rho, r2, ws
            b9(b8.inverted, b4.intermediate[0], i4.value, i7.value)

            p1(b6.performance[0])
            p2(b6.performance[1])
            p3(b9.performance)
            # p3(b9.performance[0])
            # p4(b9.performance[1])


    man = MAN(design_params, input_specs, fixed_params,
            behaviours, decisions, margin_nodes, performances, name)

    if man_folder is None:
        man_folder = '/'

    if train_surrogate:
        # train material surrogate
        variable_dict = {
            'material' : {'type' : 'ENUM', 'limits' : decision_2.universe},
        }
        b4.train_surrogate(variable_dict,n_samples=50,sm_type='KRG')
        b4.train_inverse(sm_type='LS')

        # Train surrogate for bearing load case
        variable_dict = {
            'n_struts' : {'type' : 'INT', 'limits' : [decision_3.universe[0],decision_3.universe[-1]]},
            'w' : {'type' : 'FLOAT', 'limits' : [decision_1.universe[0],decision_1.universe[-1]]},
            'h' : {'type' : 'FLOAT', 'limits' : d1.universe},
            'theta' : {'type' : 'FLOAT', 'limits' : d2.universe},
            'BX' : {'type' : 'FLOAT', 'limits' : s3.universe},
            'BY' : {'type' : 'FLOAT', 'limits' : s4.universe},
            'E' : {'type' : 'fixed', 'limits' : i2.value},
            'r1' : {'type' : 'fixed', 'limits' : i3.value},
            'r2' : {'type' : 'fixed', 'limits' : i4.value},
        }
        b5.train_surrogate(variable_dict,n_samples=1000,num_threads=num_threads)
        man.save(name,folder=man_folder)
        b5.train_inverse('n_struts',sm_type='LS', bandwidth=[1e-3])
        man.save(name,folder=man_folder)
    else:
        # load just the behaviour models
        b4.load(os.path.join(man_folder,name))
        b5.load(os.path.join(man_folder,name))
        b5.train_inverse('n_struts', sm_type='LS', bandwidth=[1e-3])
    
    if overwrite:
        man.save(name,folder=man_folder)

    return man