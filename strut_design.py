import sys
from mvm.utilities import check_folder
import os
import pandas as pd

from man_defs import get_man

# folder = os.path.join('data','strut','C1'); lean = 0.0; height = 15.0
# folder = os.path.join('data','strut','C2'); lean = 30.0; height = 15.0
# folder = os.path.join('data','strut','C3'); lean = 0.0; height = 20.0
folder = os.path.join('data','strut','C4'); lean = 30.0; height = 20.0

base_folder = os.path.join('data','strut')
check_folder(base_folder)
check_folder(folder)

# off-the shelf parts
widths = list(range(60,120+5,10))
# create material dictionary
n_materials = 15

material_dict_user = {
    'start'   : {
        'sigma_y' : 92, # MPa
        'rho' : 11.95e-06, # kg/mm3
        'cost' : 0.1 # USD/kg
        },
    'end'  : {
        'sigma_y' : 828, # MPa
        'rho' : 4.43e-06, # kg/mm3
        'cost' : 1.10 # USD/kg
        },
}
# generate 10 materials by linearly interpolating
df = pd.DataFrame(columns=material_dict_user['start'].keys(), index=range(n_materials), dtype=float)
df.iloc[0] = material_dict_user['start']
df.iloc[-1] = material_dict_user['end']
df.interpolate(method='linear',axis=0,inplace=True)
material_dict = df.transpose().to_dict()

man = get_man(height=height,lean=lean,n_materials=material_dict,widths=widths,
    train_surrogate=False,man_folder=base_folder,overwrite=True,name='strut_s',num_threads=1)

# load the MAN
man.load('strut_s',folder=base_folder)
man.init_decisions()
man.allocate_margins()
man.forward()

# Perform Monte-Carlo simulation
n_epochs = 1000
for n in range(n_epochs):
    sys.stdout.write('Progress: %d%%   \r' % ((n / n_epochs) * 100))
    sys.stdout.flush()

    man.randomize()
    man.init_decisions()
    man.allocate_margins()
    man.forward()
    man.compute_impact()
    man.compute_absorption()

man.save('strut_s',folder=folder)

# man.load('strut_s',folder=folder)

# # View distribution of excess
# man.margin_nodes[0].excess.view(xlabel='E1')
# man.margin_nodes[1].excess.view(xlabel='E2')

# # View distribution of Impact on Performance
# man.impact_matrix.view(0, 0, xlabel='E1,P1')
# man.impact_matrix.view(1, 0, xlabel='E2,P1')

# man.deterioration_vector.view(0, xlabel='S1')
# man.deterioration_vector.view(1, xlabel='S2')

# man.absorption_matrix.view(0, 0, xlabel='E1,S1')
# man.absorption_matrix.view(1, 0, xlabel='E2,S1')

# man.absorption_matrix.view(0, 1, xlabel='E1,S2')
# man.absorption_matrix.view(1, 1, xlabel='E2,S2')

# impact_matrix = np.nanmean(man.impact_matrix.values,axis=(2,))  # average along performance parameters (assumes equal weighting)
# absorption_matrix = np.nanmean(man.absorption_matrix.values,axis=(2,))  # average along input specs (assumes equal weighting)

# rows = ['E%i'%i for i in range(absorption_matrix.shape[0])]
# cols = ['S%i'%i for i in range(absorption_matrix.shape[1])]
# plt.imshow(absorption_matrix, cmap='hot', interpolation='nearest')
# plt.xticks(range(0, absorption_matrix.shape[1]),cols)
# plt.yticks(range(0, absorption_matrix.shape[0]),rows)

# plt.show()

# rows = ['E%i'%i for i in range(impact_matrix.shape[0])]
# cols = ['S%i'%i for i in range(impact_matrix.shape[1])]
# plt.imshow(impact_matrix, cmap='hot', interpolation='nearest')
# plt.xticks(range(0, impact_matrix.shape[1]),cols)
# plt.yticks(range(0, impact_matrix.shape[0]),rows)

# plt.show()

# # display the margin value plot
# man.compute_mvp('scatter')