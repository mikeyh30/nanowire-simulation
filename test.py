import numpy as np
import matplotlib.pyplot as plt

def nanomagnet_field(theta,ratio):
    norm_theta = theta%(2*np.pi)
    norm_ratio = (ratio*2*np.pi)

    new_theta = np.where(norm_theta < norm_ratio,  norm_theta/(2*ratio),
                         np.pi+(norm_theta-norm_ratio)/(2-2*ratio)
                        )

    return np.sin(new_theta)

x = np.linspace(0,6*np.pi,1000)

plt.plot(x,nanomagnet_field(x,0.2))
plt.show()