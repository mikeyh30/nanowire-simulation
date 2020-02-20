import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks


def staggered_sinusoid(theta, ratio):
    norm_theta = theta % (2 * np.pi)
    norm_ratio = ratio * 2 * np.pi

    new_theta = np.where(
        norm_theta < norm_ratio,
        norm_theta / (2 * ratio),
        np.pi + (norm_theta - norm_ratio) / (2 - 2 * ratio),
    )

    return (np.sin(new_theta), np.sin(new_theta + np.pi / 2))


def rick_sinusiod(theta):
    df = pd.read_csv("./data/rick-simulation-profiles/2pi_1D_slice.csv")
    norm_theta = (theta / (2 * np.pi)) % 1
    Xrange = df.shape[0]
    # Yrange = (np.max([df['v'].max(),np.abs(df['v'].min())]),
    #          np.max([df['u'].max(),np.abs(df['u'].min())])
    #          )
    new_theta = np.floor(Xrange * norm_theta)
    # return (df['v'].loc[new_theta]/Yrange[0], df['u'].loc[new_theta]/Yrange[1])
    return (df["v"].loc[new_theta], df["u"].loc[new_theta])


def staggered_fourier(theta, ratio):
    norm_theta = theta % (2 * np.pi)
    norm_ratio = ratio * 2 * np.pi

    def sin(theta):
        return np.sin(theta)

    return (sin(norm_theta), sin(norm_theta + np.pi / 2))


def sines(theta):
    norm_theta = theta % (2 * np.pi)
    norm_ratio = ratio * 2 * np.pi

    def fun(theta):
        return 0.5 * np.sin(theta) + 0.5 * np.sin(2 * theta)

    return (fun(norm_theta), fun(norm_theta + np.pi / 2))


def curve(theta, peaks, fft, freq, mult):
    height = 0
    peaks = peaks[len(peaks) // 2 :]
    for peak in peaks:
        temp = mult * fft.real[peak] * np.sin(2 * np.pi * freq[peak] * theta)
        height += temp
    return height


def rick_fourier(theta):
    samples = 2048
    x = np.arange(samples)
    freq = np.fft.fftfreq(samples)
    fftsin = np.fft.fft(rick_sinusiod(x)[0])
    fftcos = np.fft.fft(rick_sinusiod(x)[1])
    peakssin, _ = find_peaks(np.abs(fftsin.real), prominence=0.8)
    peakscos, _ = find_peaks(np.abs(fftcos.real), prominence=0.8)
    sin = curve(theta, peakssin, fftsin, freq, 0.005)
    cos = curve(theta, peakscos, fftcos, freq, 0.0013)
    return sin, cos


if __name__ == "__main__":
    samples = 2048
    x = np.arange(samples)
    x1 = np.linspace(0, 6 * np.pi, samples)
    # x = np.arange(1024)
    ratio = 0.2
    # plt.plot(x,staggered_sinusoid(x,0.5)[0])
    # plt.plot(x,rick_sinusiod(x)[0],x,rick_sinusiod(x)[1])
    freq = np.fft.fftfreq(samples)
    fft = np.fft.fft(rick_sinusiod(x)[0])

    fig, (axs1, axs2) = plt.subplots(1, 2)
    axs1.plot(freq, fft.real, freq, fft.imag)
    axs2.plot(x1, rick_fourier(x1)[1])
    # plt.plot(freq,np.fft.fft(np.sin(x)).imag,
    #          freq,np.fft.fft(np.sin(x)).real)
    plt.show()
