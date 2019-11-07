import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def staggered_sinusoid(theta,ratio):
    norm_theta = theta%(2*np.pi)
    norm_ratio = (ratio*2*np.pi)

    new_theta = np.where(norm_theta < norm_ratio,  norm_theta/(2*ratio),
                         np.pi+(norm_theta-norm_ratio)/(2-2*ratio)
                        )

    return (np.sin(new_theta), np.sin(new_theta+np.pi/2))

def rick_sinusiod(theta):
    df = pd.read_csv('.data/rick-simulation-profiles/2pi_1D_slice.csv')
    norm_theta = (theta/(2*np.pi))%1
    Xrange = df.shape[0]
    # Yrange = (np.max([df['v'].max(),np.abs(df['v'].min())]),
    #          np.max([df['u'].max(),np.abs(df['u'].min())])
    #          )
    new_theta = np.floor(Xrange*norm_theta)
    # return (df['v'].loc[new_theta]/Yrange[0], df['u'].loc[new_theta]/Yrange[1])
    return (df['v'].loc[new_theta], df['u'].loc[new_theta])

if __name__ == "__main__":
    x = np.linspace(0,6*np.pi,1000)
    # plt.plot(x,staggered_sinusoid(x,0.5)[0])
    plt.plot(x,rick_sinusiod(x)[0],x,rick_sinusiod(x)[1])
    plt.show()