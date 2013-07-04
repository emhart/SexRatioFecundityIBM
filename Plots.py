import matplotlib.pyplot as plt
import numpy as np
from ibmsimulation import ibm_help as ih

b = -100
c = -2
const = 1

x = np.linspace(0, 1, 20)

out = ih.holling(3,.3,x)
plt.plot(x,out)
print out
plt.show()

