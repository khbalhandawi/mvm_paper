import matplotlib.pyplot as plt
import numpy as np

from mvm import Design
from mvm import Distribution, GaussianFunc
from mvm import InputSpec, Behaviour, MarginNode, MarginNetwork
from mvm import TriangularFunc, FuzzySet, FuzzyRule, FuzzySystem

np.warnings.filterwarnings('error', category=np.VisibleDeprecationWarning)

# Generate universe variables
lb = np.array([250, 480, 1.0])
ub = np.array([450, 680, 6.0])
labels = ['$T_1$', '$T_2$', '$n_\mathrm{safety}$']

universe = np.linspace(lb, ub, 100)  # grid for all variables

# Define shape of triangular membership functions
shapes_lo = np.array([lb, lb, lb + (ub - lb) / 2])
shapes_md = np.array([lb + (ub - lb) / 4, lb + (ub - lb) / 2, ub - (ub - lb) / 4])
shapes_hi = np.array([lb + 1 * (ub - lb) / 2, ub, ub])

# Generate fuzzy membership functions
fuzzy_sets = []
for i in range(len(lb)):
    lo = TriangularFunc(universe[:, i], labels[i])
    lo.set_func(shapes_lo[0, i], shapes_lo[1, i], shapes_lo[2, i])

    md = TriangularFunc(universe[:, i], labels[i])
    md.set_func(shapes_md[0, i], shapes_md[1, i], shapes_md[2, i])

    hi = TriangularFunc(universe[:, i], labels[i])
    hi.set_func(shapes_hi[0, i], shapes_hi[1, i], shapes_hi[2, i])

    fuzzyset_i = FuzzySet(lo, md, hi, labels[i])

    fuzzy_sets += [fuzzyset_i]

# label each fuzzy set
temp_T1 = fuzzy_sets[0]
temp_T2 = fuzzy_sets[1]
n_safety = fuzzy_sets[2]

# Visualize these universes and membership functions
temp_T1.view()
temp_T2.view()
n_safety.view()

# Define fuzzy rules
rule1 = FuzzyRule([{'fun1': temp_T1.lo, 'fun2': temp_T2.hi, 'operator': 'AND'}, ], n_safety.lo, label='R1')
rule2 = FuzzyRule([{'fun1': temp_T1.hi, 'fun2': temp_T2.lo, 'operator': 'AND'}, ], n_safety.hi, label='R2')
rule3 = FuzzyRule([{'fun1': temp_T1.lo, 'fun2': temp_T2.lo, 'operator': 'AND'}, ], n_safety.md, label='R3')
rule4 = FuzzyRule([{'fun1': temp_T1.hi, 'fun2': temp_T2.hi, 'operator': 'AND'}, ], n_safety.md, label='R4')

rule5 = FuzzyRule([{'fun1': temp_T1.md, 'fun2': temp_T2.hi, 'operator': 'AND'}, ], n_safety.lo, label='R5')
rule6 = FuzzyRule([{'fun1': temp_T1.md, 'fun2': temp_T2.lo, 'operator': 'AND'}, ], n_safety.hi, label='R6')
rule7 = FuzzyRule([{'fun1': temp_T1.lo, 'fun2': temp_T2.md, 'operator': 'AND'}, ], n_safety.lo, label='R7')
rule8 = FuzzyRule([{'fun1': temp_T1.hi, 'fun2': temp_T2.md, 'operator': 'AND'}, ], n_safety.hi, label='R8')

rule9 = FuzzyRule([{'fun1': temp_T1.md, 'fun2': temp_T2.md, 'operator': 'AND'}, ], n_safety.md, label='R9')

# Define fuzzy control system
rules = [rule1, rule2, rule3, rule4, rule5, rule6, rule7, rule8, rule9]
sim = FuzzySystem([temp_T1, temp_T2], n_safety, rules)

# %% Compute for given inputs
inputs = np.array([
    370.0,  # T1
    580.0,  # T2
])
n_safety_value, aggregate, n_safety_activation = sim.compute(inputs, normalize=True)
sim.view()

# Simulate at higher resolution the control space in 2D
n_levels = 100
grid_in = Design(lb[:2], ub[:2], n_levels, "fullfact").unscale()

# Loop through the system to collect the control surface
sim.reset()
z, a, _ = sim.compute(grid_in, normalize=True)

x = grid_in[:, 0].reshape((n_levels, n_levels))
y = grid_in[:, 1].reshape((n_levels, n_levels))
z = z.reshape((n_levels, n_levels))

# Figure 2

# Plot the result in 2D
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot()

surf = ax.contourf(x, y, z, cmap=plt.cm.jet, )
ax.set_xlabel('$T_1$')
ax.set_ylabel('$T_2$')

cbar = plt.cm.ScalarMappable(cmap=plt.cm.jet)
cbar.set_array(z)

boundaries = np.linspace(1, 6, 51)
cbar_h = fig.colorbar(cbar, boundaries=boundaries)
cbar_h.set_label('safety factor', rotation=90, labelpad=3)

# fig.savefig('images/sim2D.png', format='png', dpi=200)
plt.show()

# Define PDFs for input specifications

# For T1, T2
mu = np.array([370, 580])
Sigma = np.array([
    [50, 25],
    [75, 100],
])
Requirement = GaussianFunc(mu, Sigma, 'R1')
Requirement.random(10000)
Requirement.view(xlabel='Thermal specification $T_1,T_2 \
    \sim \mathcal{N}(\mu,\Sigma)$', savefile='R')
Requirement.reset()

# for n_decided
mu = np.array([4.0, ])
Sigma = np.array([[0.3 ** 2, ], ])
decided_safety = GaussianFunc(mu, Sigma, 'T1')
decided_safety.random(10000)
decided_safety.view(xlabel='Decided $n_\mathrm{safety}^\mathrm{decided}$', savefile='Th')
decided_safety.reset()

# define input specifications
s1 = InputSpec(Requirement.mu[0], 'S1', universe=(325, 550), cov_index=0, description='nacelle temperature',
               symbol='T1', distribution=Requirement, variable_type='FLOAT')
s2 = InputSpec(Requirement.mu[1], 'S2', universe=(325, 550), cov_index=1, description='gas surface temperature',
               symbol='T2', distribution=Requirement, variable_type='FLOAT')
s3 = InputSpec(decided_safety.mu, 'S3', universe=(1.0, 6.0), description='decided safety', 
               symbol='n_decided', variable_type='FLOAT')
input_specs = [s1, s2, s3]

# define the behaviour models
behaviour = Distribution(aggregate, lb=lb[-1], ub=ub[-1], label='B1')
print("Sample mean for n_safety = %f" % (behaviour.random(10000).mean()))  # should be close to n_safety_value
print(
    "Sample standard deviation for n_safety = %f" % (behaviour.random(10000).std()))  # should be close to
# n_safety_value
behaviour.view(xlabel='Behaviour $n_\mathrm{safety}$', savefile='B')
behaviour.reset()


# this is the n_safety model
class B1(Behaviour):
    def __call__(self, T1, T2):
        # Compute for given inputs
        sim.reset()
        _, aggregate, _ = sim.compute(np.array([[T1, T2], ]), normalize=True)
        behaviour = Distribution(aggregate, lb=lb[-1], ub=ub[-1])
        self.target_threshold = behaviour.random()


b1 = B1('B1')
behaviours = [b1, ]

# Defining a MarginNode object
e1 = MarginNode('E1', direction='must_exceed')
margin_nodes = [e1, ]


# Defining a MarginNetwork object
class MAN(MarginNetwork):

    def randomize(self):
        Requirement.random()
        decided_safety.random()
        s1.random()
        s2.random()
        s3.random()

    def forward(self):
        Requirement.random()
        decided_safety.random()

        # retrieve MAN components
        s1 = self.input_specs[0]  # stochastic
        s2 = self.input_specs[1]  # stochastic
        s3 = self.input_specs[2]  # stochastic
        b1 = self.behaviours[0]
        e1 = self.margin_nodes[0]

        # Execute behaviour models
        b1(s1.value, s2.value)
        e1(b1.target_threshold, s3.value)


man = MAN([],input_specs,[],behaviours,[],margin_nodes,[],'MAN_1')

for n in range(10000):
    man.forward()

e1.excess.view(xlabel='Excess $\Delta = n_\mathrm{safety} - \
    n_\mathrm{safety}^\mathrm{decided}$')
e1.excess.view_cdf(xlabel='Excess $\Delta = n_\mathrm{safety} - \
    n_\mathrm{safety}^\mathrm{decided}$')

# for N in range(10,1000,50):
#     ThermalNode_1(Requirement.random(N),None)
#     # behaviour.random(N) # dummy behaviour

#     Requirement.view(xlabel='Thermal specification $T_1,T_2 \
#         \sim \mathcal{N}(\mu,\Sigma)$',savefile='R_%i' %(N))
#     threshold.view(xlabel='Decided $n_\mathrm{safety}^\mathrm{decided}$',savefile='Th_%i' %(N))
#     # behaviour.view(xlabel='Response $n_\mathrm{safety}$',savefile='B_%i' %(N))

#     ThermalNode_1.excess.view(xlabel='Excess $\Delta = n_\mathrm{safety} - \
#         n_\mathrm{safety}^\mathrm{decided}$', savefile='E_%i' %(N))
#     ThermalNode_1.excess.view_cdf(xlabel='Excess $\Delta = n_\mathrm{safety} - \
#         n_\mathrm{safety}^\mathrm{decided}$', savefile='E_cdf_%i' %(N))
