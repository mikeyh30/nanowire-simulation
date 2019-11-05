import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

def rick_sinusiod(theta):
    df = pd.read_csv('./2pi_1D_slice.csv')
    norm_theta = (theta/(2*np.pi))%1
    Xrange = df.shape[0]
    new_theta = np.floor(Xrange*norm_theta)
    return (df['v'].loc[new_theta], df['u'].loc[new_theta])

if __name__ == "__main__":
    x = np.linspace(0,6*np.pi,1000)
    # plt.plot(x,staggered_sinusoid(x,0.5))
    plt.plot(x,rick_sinusiod(x)[0])
    plt.show()