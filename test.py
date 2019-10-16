import numpy as np
import matplotlib.pyplot as plt

def nanomagnet_field(theta,ratio):
    norm_theta = theta%(2*np.pi)
    norm_ratio = (ratio*2*np.pi)
    if (norm_theta < norm_ratio):
        new_theta = norm_theta/(2*ratio)
    elif (norm_theta >= norm_ratio):
        new_theta = np.pi+(norm_theta-norm_ratio)/(2-2*ratio)
    else:
        raise ValueError('input to nanomagnet_field is wrong.')
    return np.sin(new_theta)

x = np.linspace(0,6*np.pi,1000)
y = []

for i in x:
    y.append(nanomagnet_field(i,0.2))

plt.plot(x,y)
plt.show()