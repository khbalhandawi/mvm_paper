"""
Examples
--------
Example showing the usage if discrete density 
function values to define a multivariate Guassian
"""

# 1D example (uniform)
import numpy as np
from mvm import GaussianFunc, UniformFunc, Distribution

mean = 10.0
sd = 5.0

dist = UniformFunc(center=np.array([mean]), interval=np.array([[sd ** 2, ]]))

print(dist.random(10000).mean())  # should be close to 10.0
print(dist.random(10000).std())  # should be close to 5.0

dist.view()

# View distribution
import matplotlib.pyplot as plt

plt.hist(dist.random(10000).squeeze(), bins=100, density=True)
# plt.plot(np.linspace(-6,6),pdf)
plt.show()

# 1D example
import numpy as np
from mvm import GaussianFunc, Distribution

mean = 10.0
sd = 5.0
x = np.linspace(-100, 100, 100)  # 2D grid

function = GaussianFunc(np.array([mean]), np.array([[sd ** 2, ]]))
p = function.compute_density(x[:, None])  # get density values
dist = Distribution(p)

print(dist.random(10000).mean())  # should be close to 10.0
print(dist.random(10000).std())  # should be close to 5.0

# View distribution
import matplotlib.pyplot as plt

plt.hist(dist.random(10000).squeeze(), bins=100, density=True)
# plt.plot(np.linspace(-6,6),pdf)
plt.show()

# Arbitrary example
import numpy as np
from mvm import GaussianFunc, Distribution

x = np.linspace(-100, 100, 100)  # 2D grid

# define a piecewise linear function
def density(x):
    density = np.empty(0)
    for value in x:
        if value <= 0:
            p = -0.1 * value
        elif value > 0 and value <= 50:
            p = 0.5 * value
        else:
            p = 50
        
        density = np.append(density,p)
    
    return density

p = density(x)
plt.plot(x,p)

# define the arbitrary distribution by providing the density values
dist = Distribution(p,lb=-100,ub=100)

print(dist.random(10000).mean())
print(dist.random(10000).std())

# View distribution
import matplotlib.pyplot as plt

plt.hist(dist.random(10000).squeeze(), bins=100, density=True)
# plt.plot(np.linspace(-6,6),pdf)
plt.show()

# 2D example
import numpy as np
from mvm import Design, GaussianFunc, Distribution

mean = 10.0
sd = 5.0

x = Design(np.array([-100, -100]), np.array([100, 100]), 512, "fullfact").unscale()  # 2D grid
function = GaussianFunc(np.array([mean, mean]), np.array([[sd ** 2, 0], [0, sd ** 2]]))
p = function.compute_density(x)  # get density values

dist = Distribution(p.reshape((512, 512)))
print(dist.random(1000000).mean(axis=1))  # should be close to 10.0
print(dist.random(1000000).std(axis=1))  # should be close to 5.0

# View distribution
import matplotlib.pyplot as plt

plt.scatter(*dist.random(1000))
plt.show()
