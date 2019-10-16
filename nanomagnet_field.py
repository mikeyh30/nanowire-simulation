import numpy as np
import matplotlib.pyplot as plt

def staggered_sinusoid(theta,ratio):
    norm_theta = theta%(2*np.pi)
    norm_ratio = (ratio*2*np.pi)

    new_theta = np.where(norm_theta < norm_ratio,  norm_theta/(2*ratio),
                         np.pi+(norm_theta-norm_ratio)/(2-2*ratio)
                        )

    return np.sin(new_theta)

def staggered_cosinusoid(theta,ratio):
    norm_theta = theta%(2*np.pi)
    norm_ratio = (ratio*2*np.pi)

    new_theta = np.where(norm_theta < norm_ratio,  norm_theta/(2*ratio),
                         np.pi+(norm_theta-norm_ratio)/(2-2*ratio)
                        )
    new_theta = new_theta + np.pi/2
    return np.sin(new_theta)

if __name__ == "__main__":
    x = np.linspace(0,6*np.pi,1000)
    plt.plot(x,staggered_sinusoid(x,0.5))
    plt.show()