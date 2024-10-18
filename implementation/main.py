import numpy as np
import matplotlib.pyplot as plt
# import cmath
# This is a sample Python script.

# Press Umschalt+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def bandgap(sim_time: int, no_cell:int = 1):
    m = np.repeat([3.98/no_cell, 1.69/no_cell, 1.69/no_cell, 3.98/no_cell], no_cell)
    k = np.repeat([70.9e9, 5.28e9, 5.28e9, 70.9e9], no_cell)
    no_ele = len(m)
    w = np.zeros(no_ele)
    a = np.zeros((no_ele,no_ele), dtype=complex)
    wave_number = np.zeros(sim_time)
    frequency = np.zeros((sim_time,no_ele))
    c = k / m
    x = np.linspace(0, np.pi/no_ele, sim_time)
    for t in range(sim_time):
        for i in range(no_ele):
            l = i - 1
            r = i + 1
            if i == 0:
                w[i] = (k[no_ele - 1] + k[i]) / m[i]
                l = no_ele-1
            else:
                w[i] = (k[i - 1] + k[i]) / m[i]

            if i == no_ele-1:
                r = 0
            a[i, i] = w[i]
            a[i, r] = -c[i] * np.exp(x[t]*1j)
            a[i, l] = -(w[i] - c[i]) * np.exp(-1j * x[t])
        #
        wave_number[t] = no_ele * x[t]
    #
        frequency[t, :] = np.sqrt(np.real(np.linalg.eig(a)[0])) / (2000 * np.pi)
    fig, ax = plt.subplots()
    ax.set_ylim(0, 16)
    ax.set_xlim(0, np.pi)
    ax.plot(wave_number, frequency) #  label=['qw','ee','ew','ee'])
    ax.set_xlabel('WAVE NUMBER')
    ax.set_ylabel('FREQUENCY in kHz')
    ax.legend(['qw','ee','ew','ee'])
    plt.show()
    return wave_number, frequency


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    _, _ = bandgap(no_cell=1, sim_time=200)
    # print(wave_no)
    # print(freq)
