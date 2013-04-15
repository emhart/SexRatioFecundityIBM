import matplotlib.pyplot as plt
import numpy as np
from ibmsimulation import ibm_help as ih

b = -10
c = -2
const = 1

x = np.linspace(0, 10, 20)

out = ih.gompertz(const,b,c,x)
plt.plot(x,out)
plt.show()

